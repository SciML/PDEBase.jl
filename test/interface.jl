using PDEBase
using SciMLBase
using ModelingToolkit
using Test

struct PlainDiscretization <: SciMLBase.AbstractDiscretization end
struct TestSpace <: PDEBase.AbstractDiscreteSpace end
struct TestMapping <: PDEBase.AbstractVarEqMapping end
struct TestDerivativeData <: PDEBase.AbstractDifferentialDiscretizer end
struct TestState <: PDEBase.AbstractDiscretizationState
    equations::Vector{Any}
end
struct TestMetadata <: SciMLBase.AbstractDiscretizationMetadata{false}
    metadata::Base.RefValue{Any}
end

struct PreflightDiscretization <: PDEBase.AbstractEquationSystemDiscretization end
struct PreflightError <: Exception end

PDEBase.interface_errors(::PDESystem, ::PreflightDiscretization) = throw(PreflightError())
PDEBase.get_time(::PreflightDiscretization) = error("preflight did not run first")

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
    @test isempty(PDEBase.get_system_inputs(space))
    @test PDEBase.get_eqvar(mapping, eq) === nothing
    @test PDEBase.add_metadata!(metadata, :symbolic_system) === :symbolic_system
    @test metadata.metadata[] === :symbolic_system
end

@testset "Original-system preflight runs first" begin
    @parameters t x
    @variables u(..)
    pdesys = PDESystem(
        [Differential(t)(u(t, x)) ~ 0], Equation[],
        [t ∈ (0, 1), x ∈ (0, 1)], [t, x], [u(t, x)];
        name = :pdebase_preflight_test
    )

    @test_throws PreflightError SciMLBase.symbolic_discretize(
        pdesys, PreflightDiscretization()
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
