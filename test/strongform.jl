using PDEBase
using PDEBase.StrongForm: StateSlots, check_whitelist, split_time
using DomainSets: ClosedInterval
using Symbolics
using SymbolicUtils: Operator, term
using ModelingToolkit: @named, @parameters, PDESystem
using IntervalSets: var".."
using Test

struct Doubler <: Operator end
(D::Doubler)(ex) = Num(term(D, Symbolics.unwrap(ex); type = Real))

@variables t x y
@variables u(..) w(..) k(..)
@register_symbolic opaque(x)
Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
U, Wf = u(t, x, y), w(t, x, y)

@testset "StateSlots" begin
    slots = StateSlots(2, 3, 1)
    @test length(slots.U) == length(slots.Ud) == 2
    @test length(slots.G) == 2 * 3
    @test length(slots.H) == 2 * 3 * 3
    @test length(slots.X) == length(slots.N) == 3
    @test length(slots.P) == 1
    @test slots.T isa Num
    @test isempty(StateSlots(1, 2, 0).P)
    @test string(slots.G[4]) == "G₄"
end

@testset "Operator whitelist" begin
    ivs, dvs = [t, x, y], [U, Wf]
    kc = 1 + x^2 * y
    accepted = [
        Dx(kc * Dx(U)) + Dy(kc * Dy(U)),
        Dx((1 + U^2) * Dx(U)),
        -U * Dx(U) - Dx(U^2) + 2U,
        Dx(ifelse(x < 0.5, 1.0, 3.0) * Dx(U)) + abs(x - 0.5) * U,
        Dx(Dy(U)) + Dy(Dy(U)),
        Dt(Wf) - Dx(Dx(U) + U * Dx(Wf)) + exp(-t) * sin(x),
        Dx(opaque(x) * Dx(U)),
        Dx(k(x, y) * Dx(U)),
        ifelse((x < 0.5) & !(y > 0.2), U, max(U, 0)),
    ]
    for ex in accepted
        @test check_whitelist(ex, ivs, dvs) === nothing
    end
    @test check_whitelist(Dt(U) ~ Dx(Dx(U)), ivs, dvs) === nothing

    @variables z
    integral = Integral(x in ClosedInterval(0, 1))
    @test_throws "operator Integral(x, 0 .. 1) is not supported by the strong-form analysis" check_whitelist(
        Dt(U) - integral(U), ivs, dvs
    )
    @test_throws "operator Differential(z, 1) is not supported by the strong-form analysis" check_whitelist(
        Differential(z)(U), ivs, dvs
    )
    @test_throws r"operator (.*\.)?Doubler\(\) is not supported by the strong-form analysis" check_whitelist(
        Doubler()(U) + U, ivs, dvs
    )
    # A registered function without a derivative rule is rejected only when the
    # state reaches it, since the tangent kernels need its partial derivatives.
    @test_throws "operator opaque is not supported by the strong-form analysis" check_whitelist(
        opaque(U) * Dx(U), ivs, dvs
    )
    @test check_whitelist(opaque(x) * Dx(U), ivs, dvs) === nothing
    # The same holds for a symbolic function that is not a dependent variable.
    @test_throws "operator k is not supported by the strong-form analysis" check_whitelist(
        Dx(k(U) * Dx(U)), ivs, dvs
    )
    @test check_whitelist(k(x, y) * Dx(U) + U * Wf, ivs, dvs) === nothing
    # Indexed parameter arrays carry no state and keep their argument shapes.
    @variables a[1:2]
    @test check_whitelist(a[1] * U + a[2] * Dx(U), ivs, dvs) === nothing
    @test check_whitelist(sin(a[1] * U) * Dx(U), ivs, dvs) === nothing
end

@testset "split_time" begin
    # Heat equation.
    mass, spatial = split_time(Dt(U) ~ Dx(Dx(U)) + Dy(Dy(U)), [U], t)
    @test isequal(mass, [Num(1)])
    @test isequal(sum(spatial), -Dx(Dx(U)) - Dy(Dy(U)))

    # Two fields with a coefficient on one mass term and a coupled flux.
    mass, spatial = split_time(2x * Dt(U) - Dt(Wf) / 3 ~ Dx(Dx(U) + U * Dx(Wf)), [U, Wf], t)
    @test isequal(mass, [2x, Num(-1 // 3)])
    @test length(spatial) == 2  # derivatives are distributed over the sum

    # Wave as a first-order system.
    @test isequal(first(split_time(Dt(U) ~ Wf, [U, Wf], t)), [Num(1), Num(0)])
    mass, spatial = split_time(Dt(Wf) ~ Dx(Dx(U)), [U, Wf], t)
    @test isequal(mass, [Num(0), Num(1)])
    @test isequal(only(spatial), -Dx(Dx(U)))

    # Algebraic row: no mass at all.
    mass, spatial = split_time(0 ~ Dx(Dx(Wf)) + exp(-t) * sin(x), [U, Wf], t)
    @test isequal(mass, [Num(0), Num(0)])
    @test length(spatial) == 2

    # Steady equation, with and without a time variable.
    V = u(x, y)
    for tv in (nothing, t)
        mass, spatial = split_time(Dx(Dx(V)) + Dy(Dy(V)) ~ 1, [V], tv)
        @test isequal(mass, [Num(0)])
        @test length(spatial) == 3
    end

    # Quotients with an additive numerator, which `expand` leaves whole.
    mass, spatial = split_time((Dt(U) + U) / x, [U], t)
    @test isequal(mass, [1 / x])
    @test isequal(only(spatial), U / x)
    @test isequal(first(split_time((Dt(U) + Dt(Wf)) / x ~ 0, [U, Wf], t)), [1 / x, 1 / x])

    # Residual form is accepted as well as an equation.
    @test isequal(first(split_time(Dt(U) - Dx(Dx(U)), [U], t)), [Num(1)])

    @test_throws r"second order in time; write as a first-order system" split_time(
        Dt(Dt(U)) ~ Dx(Dx(U)), [U], t
    )
    @test_throws r"is not ∂ₜu × coefficient" split_time(Dt(U)^2 ~ Dx(Dx(U)), [U], t)
    @test_throws r"is not ∂ₜu × coefficient" split_time(Dx(Dt(U)) ~ Dx(Dx(U)), [U], t)
    @test_throws r"is not ∂ₜu × coefficient" split_time(Dt(U) * Dt(Wf) ~ 0, [U, Wf], t)
    @test_throws r"state-dependent mass .* needs residual-form DAE" split_time(
        (1 + U^2) * Dt(U) ~ Dx(Dx(U)), [U], t
    )
end

@testset "split_time on flattened array variables" begin
    @parameters τ ξ
    @variables (q(..))[1:2]
    Dτ, Dξ = Differential(τ), Differential(ξ)
    eqs = [Dτ(q(τ, ξ)[1]) ~ Dξ(Dξ(q(τ, ξ)[1])), 0 ~ Dξ(Dξ(q(τ, ξ)[2])) - q(τ, ξ)[1]]
    bcs = [q(0, ξ)[1] ~ 0, q(τ, 0)[1] ~ 0, q(τ, 1)[1] ~ 0, q(τ, 0)[2] ~ 0, q(τ, 1)[2] ~ 0]
    @named sys = PDESystem(eqs, bcs, [τ ∈ (0 .. 1), ξ ∈ (0 .. 1)], [τ, ξ], [q(τ, ξ)])
    sys, replaced = PDEBase.make_pdesys_compatible(sys)
    dvs = sort(PDEBase.depvars(PDEBase.VariableMap(sys, nothing; replaced_vars = replaced)); by = string)
    flat = sys.eqs
    check_whitelist.(flat, Ref([τ, ξ]), Ref(dvs))
    @test isequal(first(split_time(flat[1], dvs, τ)), [Num(1), Num(0)])
    @test isequal(first(split_time(flat[2], dvs, τ)), [Num(0), Num(0)])
end

@testset "StrongForm developer API is documented" begin
    api = (
        :StateSlots, :IBPPolicy, :ProductRule, :Weighted, :Reject, :EssentialBC, :NaturalBC,
        :PeriodicBC, :StrongFormAnalysis, :check_whitelist, :split_time, :split_divergence,
        :distribute_derivatives, :NONSMOOTH_FUNCTIONS, :DISCONTINUOUS_FUNCTIONS, :SlotContext,
        :slot_volume, :check_evaluable, :analyse, :BoxFace, :slot_condition, :restrict_to_face,
        :match_flux, :route_natural, :essential_condition, :periodic_condition,
        :boundary_conditions, :pair_equations,
    )
    for name in api
        @test haskey(Base.Docs.meta(PDEBase.StrongForm), Base.Docs.Binding(PDEBase.StrongForm, name))
        isdefined(Base, :ispublic) && @test Base.ispublic(PDEBase.StrongForm, name)
    end
    @test !any(name -> Base.isexported(PDEBase, name), (:StrongForm, api...))
    @test [name for name in api if Base.isexported(PDEBase.StrongForm, name)] == [:analyse]
end
