using PDEBase
using ModelingToolkit
using ModelingToolkitBase: complete, guesses, initial_conditions, initialization_equations,
    inputs, mtkcompile, parameters, unknowns
using SciMLBase
using Symbolics
using SymbolicUtils: unwrap_const
using Test

struct UnknownsDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
end
struct ScalarUnknownsSpace <: PDEBase.AbstractDiscreteSpace
    discvars::Any
end
struct ArrayUnknownsSpace <: PDEBase.AbstractDiscreteSpace
    discvars::Any
    arrays::Any
end
struct InputUnknownsSpace <: PDEBase.AbstractDiscreteSpace
    discvars::Any
    state_unknowns::Any
    system_inputs::Any
end
struct UnknownsMetadata <: SciMLBase.AbstractDiscretizationMetadata{true}
    pdesys::Any
end

PDEBase.get_time(disc::UnknownsDiscretization) = disc.time
PDEBase.get_discvars(s::ScalarUnknownsSpace) = s.discvars
PDEBase.get_discvars(s::ArrayUnknownsSpace) = s.discvars
PDEBase.get_discvars(s::InputUnknownsSpace) = s.discvars
PDEBase.get_system_unknowns(s::ArrayUnknownsSpace) = s.arrays
PDEBase.get_system_unknowns(s::InputUnknownsSpace) = s.state_unknowns
PDEBase.get_system_inputs(s::InputUnknownsSpace) = s.system_inputs

@testset "get_system_unknowns chooses the unknowns of the generated system" begin
    @parameters t x
    @variables u(..)
    Dt = Differential(t)
    Dx = Differential(x)
    pdesys = PDESystem(
        [Dt(u(t, x)) ~ Dx(Dx(u(t, x)))],
        [u(0, x) ~ 0, u(t, 0) ~ 0, u(t, 1) ~ 0],
        [t ∈ (0, 1), x ∈ (0, 1)],
        [t, x],
        [u(t, x)];
        name = :pdebase_system_unknowns_test
    )

    @variables w(t)[1:3]
    elements = collect(w)
    state = PDEBase.EquationState(
        [Dt(w[1]) ~ -w[1], Dt(w[2]) ~ w[1] - w[2], Dt(w[3]) ~ w[2] - w[3]], Equation[]
    )
    u0 = [w[i] => Float64(i) for i in 1:3]
    disc = UnknownsDiscretization(t)
    metadata = UnknownsMetadata(pdesys)

    scalar_space = ScalarUnknownsSpace(Dict(:w => elements))
    @test isequal(PDEBase.get_system_unknowns(scalar_space), vec(elements))
    scalar_sys, _ = PDEBase.generate_system(
        state, scalar_space, u0, (0.0, 1.0), metadata, disc
    )
    @test length(unknowns(scalar_sys)) == 3
    @test isempty(inputs(scalar_sys))

    array_space = ArrayUnknownsSpace(Dict(:w => elements), [w])
    @test isequal(PDEBase.get_system_unknowns(array_space), [w])
    array_sys, _ = PDEBase.generate_system(
        state, array_space, u0, (0.0, 1.0), metadata, disc
    )
    @test length(unknowns(array_sys)) == 1
    @test isequal(only(unknowns(array_sys)), unwrap(w))

    # The state vector is the same whichever representation the unknowns take.
    scalar_prob = ODEProblem(complete(scalar_sys), u0, (0.0, 1.0))
    array_prob = ODEProblem(complete(array_sys), u0, (0.0, 1.0))
    @test scalar_prob.u0 == array_prob.u0 == [1.0, 2.0, 3.0]
    @test scalar_prob.f(scalar_prob.u0, scalar_prob.p, 0.0) ==
        array_prob.f(array_prob.u0, array_prob.p, 0.0) == [-1.0, -1.0, -1.0]
end

@testset "get_system_inputs declares external discrete fields" begin
    @parameters t x
    @variables u(..) a(t) b(t) f1(t) f2(t)
    Dt = Differential(t)
    Dx = Differential(x)
    pdesys = PDESystem(
        [Dt(u(t, x)) ~ Dx(Dx(u(t, x)))],
        [u(0, x) ~ 0, u(t, 0) ~ 0, u(t, 1) ~ 0],
        [t ∈ (0, 1), x ∈ (0, 1)],
        [t, x],
        [u(t, x)];
        name = :pdebase_system_inputs_test
    )

    state = PDEBase.EquationState(
        [Dt(a) ~ -a + f1, Dt(b) ~ a - b + f2], Equation[]
    )
    space = InputUnknownsSpace(
        Dict(:u => [a, b, f1, f2]), [a, b, f1, f2], [f1, f2]
    )
    u0 = [a => 1.0, b => 2.0, f1 => 3.0, f2 => 4.0]
    sys, _ = PDEBase.generate_system(
        state, space, u0, (0.0, 1.0), UnknownsMetadata(pdesys),
        UnknownsDiscretization(t)
    )

    @test isequal(unknowns(sys), unwrap.([a, b, f1, f2]))
    @test isequal(space.state_unknowns, [a, b, f1, f2])
    @test isequal(inputs(sys), unwrap.([f1, f2]))
    @test unwrap_const(initial_conditions(sys)[unwrap(f1)]) == 3.0
    @test unwrap_const(initial_conditions(sys)[unwrap(f2)]) == 4.0
    @test length(initialization_equations(sys)) == 2
    @test all(!haskey(guesses(sys), unwrap(f)) for f in (f1, f2))

    compiled = mtkcompile(sys; inputs = inputs(sys))
    @test isequal(inputs(compiled), unwrap.([f1, f2]))
    @test all(f -> any(isequal(f), parameters(compiled)), unwrap.([f1, f2]))
    @test all(f -> !any(isequal(f), unknowns(compiled)), unwrap.([f1, f2]))

    prob = ODEProblem(compiled, nothing, (0.0, 1.0))
    state_values = Dict(unknowns(compiled) .=> prob.u0)
    rates = Dict(unknowns(compiled) .=> prob.f(prob.u0, prob.p, 0.0))
    @test state_values[unwrap(a)] == 1.0
    @test state_values[unwrap(b)] == 2.0
    @test rates[unwrap(a)] == 2.0
    @test rates[unwrap(b)] == 3.0

    bad_space = InputUnknownsSpace(space.discvars, [a, b], [f1, f2])
    @test_throws ArgumentError PDEBase.generate_system(
        state, bad_space, u0, (0.0, 1.0), UnknownsMetadata(pdesys),
        UnknownsDiscretization(t)
    )
end
