using PDEBase
using ModelingToolkit
using ModelingToolkitBase: complete, unknowns, get_eqs
using Symbolics
using SymbolicUtils
using Test

const DAE_HELPERS = (
    :BrownFullBasicInitUnsafeError, :brown_init_offenders, :differential_unknowns,
    :dae_operating_point, :time_derivative_order,
)

@testset "DAE support names are public and documented" begin
    for name in DAE_HELPERS
        @test isdefined(PDEBase, name)
        @test haskey(Base.Docs.meta(PDEBase), Base.Docs.Binding(PDEBase, name))
        isdefined(Base, :ispublic) && @test Base.ispublic(PDEBase, name)
    end
end

@testset "time_derivative_order" begin
    @independent_variables τ
    @independent_variables s
    @variables a(τ)
    D = Differential(τ)
    t = Symbolics.unwrap(τ)
    @test PDEBase.time_derivative_order(a, t) == 0
    @test PDEBase.time_derivative_order(D(a), t) == 1
    @test PDEBase.time_derivative_order(D(D(a)) + a, t) == 2
    @test PDEBase.time_derivative_order(2 * D(a)^2, t) == 1
    @test PDEBase.time_derivative_order(Differential(s)(a), t) == 0
    @test PDEBase.time_derivative_order(1.0, t) == 0
end

# Scalar unknowns: `a` is differential, `b` algebraic.
@testset "scalar system" begin
    @independent_variables τ
    @variables a(τ) b(τ)
    D = Differential(τ)
    eqs = [D(a) ~ -a + b, b ~ 2a]
    build(init) = complete(
        System(eqs, τ, [a, b], []; initialization_eqs = init, name = :sys)
    )

    sys = build([a ~ 1.0])
    unks = Set{Any}(Symbolics.unwrap.(unknowns(sys)))
    diffvars = PDEBase.differential_unknowns(sys, unks)
    @test diffvars == Dict{Any, Any}(Symbolics.unwrap(a) => true)

    @test isempty(PDEBase.brown_init_offenders(sys))
    @test isempty(PDEBase.brown_init_offenders(build(Equation[])))

    algebraic = PDEBase.brown_init_offenders(build([b ~ 1.0]))
    @test length(algebraic) == 1
    @test occursin("algebraic rather than differential", last(only(algebraic)))

    coupled = PDEBase.brown_init_offenders(build([a ~ b]))
    @test length(coupled) == 1
    @test occursin("relates", last(only(coupled)))

    derivative = PDEBase.brown_init_offenders(build([D(a) ~ 0.0]))
    @test length(derivative) == 1
    @test occursin("time derivative", last(only(derivative)))

    err = PDEBase.BrownFullBasicInitUnsafeError(derivative)
    @test err isa Exception
    msg = sprint(showerror, err)
    @test occursin("BrownFullBasicInit", msg)
    @test occursin("discretize(pdesys, disc)", msg)
    @test occursin(string(first(first(derivative))), msg)

    op = PDEBase.dae_operating_point(sys, Dict(a => 1.0), Symbolics.unwrap(τ))
    @test op[Symbolics.unwrap(a)] == 1.0
    @test op[Symbolics.unwrap(b)] == 0.0
    @test op[Differential(Symbolics.unwrap(τ))(Symbolics.unwrap(a))] == 0.0
end

# An array unknown whose interior is differential and whose two end elements are algebraic.
@testset "array system" begin
    @independent_variables τ
    @variables u(τ)[1:5]
    D = Differential(τ)
    eqs = [
        D(u[2:4]) ~ -u[2:4],
        u[1:1] ~ [0.0],
        u[5:5] ~ [0.0],
    ]
    build(init) = complete(
        System(eqs, τ, [u], []; initialization_eqs = init, name = :sys)
    )
    sys = build(Equation[])
    unks = Set{Any}(Symbolics.unwrap.(unknowns(sys)))
    uu = Symbolics.unwrap(u)
    @test uu in unks

    diffvars = PDEBase.differential_unknowns(sys, unks)
    @test length(diffvars) == 1 && isequal(only(keys(diffvars)), uu)
    @test diffvars[uu] == BitVector([false, true, true, true, false])

    @test isempty(PDEBase.brown_init_offenders(sys))

    interior = PDEBase.brown_init_offenders(build([u[2:4] ~ [1.0, 2.0, 3.0]]))
    @test isempty(interior)
    element = PDEBase.brown_init_offenders(build([u[3] ~ 1.0]))
    @test isempty(element)

    boundary = PDEBase.brown_init_offenders(build([u[1] ~ 1.0]))
    @test length(boundary) == 1
    @test occursin("algebraic rather than differential", last(only(boundary)))

    straddling = PDEBase.brown_init_offenders(build([u[1:3] ~ [1.0, 2.0, 3.0]]))
    @test length(straddling) == 1
    @test occursin("algebraic rather than differential", last(only(straddling)))

    coupled = PDEBase.brown_init_offenders(build([u[2] ~ u[3]]))
    @test length(coupled) == 1
    @test occursin("relates", last(only(coupled)))

    op = PDEBase.dae_operating_point(sys, Dict(u => [1.0, 2.0, 3.0, 4.0, 5.0]), Symbolics.unwrap(τ))
    @test op[uu] == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test op[Differential(Symbolics.unwrap(τ))(uu)] == zeros(5)
end

@testset "array system with colon slices" begin
    @independent_variables τ
    @variables u(τ)[1:3, 1:5]
    D = Differential(τ)
    sys = complete(
        System([D(u[:, 2:4]) ~ -u[:, 2:4]], τ, [u], []; name = :sys)
    )
    uu = Symbolics.unwrap(u)
    unks = Set{Any}(Symbolics.unwrap.(unknowns(sys)))

    diffvars = PDEBase.differential_unknowns(sys, unks)
    expected = falses(3, 5)
    expected[:, 2:4] .= true
    @test diffvars[uu] == expected

    differential = PDEBase._is_differential(Symbolics.unwrap(u[:, 3:4]), diffvars, unks)
    @test differential === true
    @test PDEBase._is_differential(Symbolics.unwrap(u[:, 1:2]), diffvars, unks) === false
end
