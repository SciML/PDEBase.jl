using PDEBase
using SciMLBase
using ModelingToolkit
using Test

struct PlainDiscretization <: SciMLBase.AbstractDiscretization end
struct PreflightDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
end
const PREFLIGHT_SYSTEM = Ref{Any}(nothing)

PDEBase.get_time(disc::PreflightDiscretization) = disc.time
function PDEBase.interface_errors(sys::PDESystem, ::PreflightDiscretization)
    PREFLIGHT_SYSTEM[] = sys
    throw(ArgumentError("preflight sentinel"))
end

struct TestSpace <: PDEBase.AbstractDiscreteSpace end
struct TestMapping <: PDEBase.AbstractVarEqMapping end
struct TestDerivativeData <: PDEBase.AbstractDifferentialDiscretizer end
struct TestState <: PDEBase.AbstractDiscretizationState
    equations::Vector{Any}
end
struct TestMetadata <: SciMLBase.AbstractDiscretizationMetadata{false}
    metadata::Base.RefValue{Any}
end

@testset "Default discretization interface" begin
    @parameters t x
    @variables u(..)
    Dt = Differential(t)
    Dx = Differential(x)
    eq = Dt(u(t, x)) ~ Dx(Dx(u(t, x)))
    pdesys = PDESystem(
        [eq],
        [u(0, x) ~ 0, u(t, 0) ~ 0, u(t, 1) ~ 0],
        [t ∈ (0, 1), x ∈ (0, 1)],
        [t, x],
        [u(t, x)];
        name = :pdebase_interface_test
    )
    v = VariableMap(pdesys)
    plain = PlainDiscretization()
    space = TestSpace()
    mapping = TestMapping()
    derivative_data = TestDerivativeData()
    state = TestState(Any[])
    metadata = TestMetadata(Ref{Any}(nothing))

    @test PDEBase.interface_errors(pdesys, v, plain) === nothing
    @test PDEBase.interface_errors(pdesys, plain) === nothing
    @test PDEBase.check_boundarymap(Dict(), v, plain) === nothing
    @test PDEBase.should_transform(pdesys, plain, Dict()) === false
    @test PDEBase.transform_pde_system!(v, Dict(), pdesys, plain) === nothing
    @test isempty(PDEBase.construct_disc_state(plain))
    @test PDEBase.construct_discrete_space(v, plain) === nothing
    @test PDEBase.construct_var_equation_mapping([], Dict(), space, plain) === nothing
    @test PDEBase.construct_differential_discretizer(pdesys, space, plain, Dict()) === nothing
    @test PDEBase.discretize_equation!(
        state, eq, mapping, nothing, Dict(), Any[], space, derivative_data, Dict(), plain
    ) === nothing
    @test PDEBase.generate_ic_defaults([], space, plain) === nothing
    @test PDEBase.generate_metadata(space, plain, pdesys, Dict(), nothing, []) === nothing
    @test PDEBase.generate_system(state, space, [], nothing, metadata, plain; checks = false) === nothing
    @test PDEBase.get_time(plain) === nothing
    @test isempty(PDEBase.get_discvars(space))
    @test PDEBase.get_eqvar(mapping, eq) === nothing
    @test PDEBase.add_metadata!(metadata, :symbolic_system) === :symbolic_system
    @test metadata.metadata[] === :symbolic_system
end

@testset "Original-system preflight" begin
    @parameters t x
    @variables state(..)
    @variables forcing(..) [input = true]
    Dt = Differential(t)
    domains = [t ∈ (0, 1), x ∈ (0, 1)]
    pdesys = PDESystem(
        [Dt(state(t, x)) ~ forcing(t, x)],
        [state(0, x) ~ 0, forcing(0, x) ~ 1],
        domains,
        [t, x],
        [state(t, x), forcing(t, x)];
        name = :preflight_test
    )

    PREFLIGHT_SYSTEM[] = nothing
    @test_throws "preflight sentinel" SciMLBase.symbolic_discretize(
        pdesys, PreflightDiscretization(t)
    )
    @test PREFLIGHT_SYSTEM[] === pdesys

    differentiated_input = PDESystem(
        [Dt(state(t, x)) ~ Dt(forcing(t, x))],
        [state(0, x) ~ 0],
        domains,
        [t, x],
        [state(t, x), forcing(t, x)];
        name = :invalid_input_derivative
    )
    PREFLIGHT_SYSTEM[] = nothing
    @test_throws "cannot be differentiated" SciMLBase.symbolic_discretize(
        differentiated_input, PreflightDiscretization(t)
    )
    @test PREFLIGHT_SYSTEM[] === nothing

    timeless_input = PDESystem(
        [Dt(state(t, x)) ~ forcing(x)],
        [state(0, x) ~ 0],
        domains,
        [t, x],
        [state(t, x), forcing(x)];
        name = :timeless_input
    )
    @test_throws "must depend on time" SciMLBase.symbolic_discretize(
        timeless_input, PreflightDiscretization(t)
    )

    @variables array_input(..)[1:2] [input = true]
    array_input_system = PDESystem(
        [Dt(state(t, x)) ~ sum(array_input(t, x))],
        [state(0, x) ~ 0],
        domains,
        [t, x],
        [state(t, x), array_input(t, x)];
        name = :array_input
    )
    @test_throws "real and scalar-valued" SciMLBase.symbolic_discretize(
        array_input_system, PreflightDiscretization(t)
    )
end

struct MockDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
end
struct MockSpace <: PDEBase.AbstractDiscreteSpace
    eqvar::Any
end
struct MockMapping <: PDEBase.AbstractVarEqMapping
    eqvar::Any
end
struct MockDerivativeData <: PDEBase.AbstractDifferentialDiscretizer end
struct MockState <: PDEBase.AbstractDiscretizationState
    equations::Vector{Any}
end
struct MockMetadata <: SciMLBase.AbstractDiscretizationMetadata{true}
    metadata::Base.RefValue{Any}
end

PDEBase.get_time(disc::MockDiscretization) = disc.time
PDEBase.construct_disc_state(::MockDiscretization) = MockState(Any[])
PDEBase.construct_discrete_space(v::VariableMap, ::MockDiscretization) = MockSpace(first(depvars(v)))
PDEBase.construct_var_equation_mapping(_, _, space::MockSpace, ::MockDiscretization) =
    MockMapping(space.eqvar)
PDEBase.construct_differential_discretizer(_, _, ::MockDiscretization, _) = MockDerivativeData()
PDEBase.get_eqvar(mapping::MockMapping, _) = mapping.eqvar
function PDEBase.discretize_equation!(
        state::MockState, pde::Equation, ::MockMapping, _, _, _, ::MockSpace,
        ::MockDerivativeData, _, ::MockDiscretization
    )
    push!(state.equations, pde)
    return nothing
end
PDEBase.generate_ic_defaults(_, ::MockSpace, ::MockDiscretization) = :u0
PDEBase.generate_metadata(_, ::MockDiscretization, _, _, _, _) = MockMetadata(Ref{Any}(nothing))
function PDEBase.generate_system(
        state::MockState, space::MockSpace, u0, tspan, metadata::MockMetadata,
        ::MockDiscretization; checks = true
    )
    return state, space, u0, tspan, metadata, checks
end
PDEBase.get_discvars(space::MockSpace) = (space.eqvar,)

@testset "Generic discretization hook dispatch" begin
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
        name = :pdebase_generic_interface_test
    )

    result = SciMLBase.symbolic_discretize(pdesys, MockDiscretization(t); checks = false)
    state, space, u0, tspan, metadata, checks = result

    @test length(state.equations) == 1
    @test space isa MockSpace
    @test PDEBase.get_discvars(space) == (space.eqvar,)
    @test u0 === :u0
    @test tspan == (0, 1)
    @test metadata isa MockMetadata
    @test checks === false
    @test PDEBase.get_time(MockDiscretization(t)) === t
end
