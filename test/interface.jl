using PDEBase
using SciMLBase
using ModelingToolkit
using SymbolicUtils: operation
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

struct MockDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
end
struct MockBoxDiscretization <: PDEBase.AbstractBoxDiscretization
    time::Any
end
struct LoweringDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
    lowered::Any
    received::Base.RefValue{Any}
end
const MockDisc = Union{MockDiscretization, MockBoxDiscretization, LoweringDiscretization}
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

PDEBase.get_time(disc::MockDisc) = disc.time
PDEBase.construct_disc_state(::MockDisc) = MockState(Any[])
PDEBase.construct_discrete_space(v::VariableMap, ::MockDisc) = MockSpace(first(depvars(v)))
PDEBase.construct_var_equation_mapping(_, _, space::MockSpace, ::MockDisc) =
    MockMapping(space.eqvar)
PDEBase.construct_differential_discretizer(_, _, ::MockDisc, _) = MockDerivativeData()
PDEBase.get_eqvar(mapping::MockMapping, _) = mapping.eqvar
function PDEBase.discretize_equation!(
        state::MockState, pde::Equation, ::MockMapping, _, _, _, ::MockSpace,
        ::MockDerivativeData, _, ::MockDisc
    )
    push!(state.equations, pde)
    return nothing
end
PDEBase.generate_ic_defaults(_, ::MockSpace, ::MockDisc) = :u0
PDEBase.generate_metadata(_, ::MockDisc, _, _, _, _) = MockMetadata(Ref{Any}(nothing))
function PDEBase.generate_system(
        state::MockState, space::MockSpace, u0, tspan, metadata::MockMetadata,
        ::MockDisc; checks = true
    )
    return state, space, u0, tspan, metadata, checks
end
PDEBase.get_discvars(space::MockSpace) = (space.eqvar,)
function PDEBase.lower_boundaries(pdesys::PDESystem, disc::LoweringDiscretization)
    disc.received[] = pdesys
    return disc.lowered
end
function PDEBase.split_ics_bcs(boundarymap, v::VariableMap, disc::LoweringDiscretization)
    t = PDEBase.get_time(disc)
    ics = mapreduce(u -> boundarymap[u][t], vcat, keys(boundarymap))
    bcmap = Dict(
        u => filter(kv -> !isequal(first(kv), t), boundarymap[u])
            for u in keys(boundarymap)
    )
    return ics, bcmap
end

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

function heat_pdesys(; coefficient = 1)
    @parameters t x
    @variables u(..)
    Dt = Differential(t)
    Dx = Differential(x)
    return PDESystem(
            [Dt(u(t, x)) ~ coefficient * Dx(Dx(u(t, x)))],
            [u(0, x) ~ 0, u(t, 0) ~ 0, u(t, 1) ~ 0],
            [t ∈ (0, 1), x ∈ (0, 1)],
            [t, x],
            [u(t, x)];
            name = :pdebase_hierarchy_test
        ), t
end

@testset "Box and mesh discretization hierarchy" begin
    @test PDEBase.AbstractBoxDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    @test PDEBase.AbstractMeshDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    @test PDEBase.AbstractEquationSystemDiscretization <: SciMLBase.AbstractDiscretization
    @test !(PDEBase.AbstractBoxDiscretization <: PDEBase.AbstractMeshDiscretization)
    @test !(PDEBase.AbstractMeshDiscretization <: PDEBase.AbstractBoxDiscretization)
    @test !(PDEBase.AbstractBoxDiscretization <: PDEBase.AbstractOptimizationSystemDiscretization)
    @test !(PDEBase.AbstractMeshDiscretization <: PDEBase.AbstractOptimizationSystemDiscretization)
    @test AbstractBoxDiscretization === PDEBase.AbstractBoxDiscretization
    @test AbstractMeshDiscretization === PDEBase.AbstractMeshDiscretization

    pdesys, t = heat_pdesys()
    generic = SciMLBase.symbolic_discretize(pdesys, MockDiscretization(t); checks = false)
    pdesys, t = heat_pdesys()
    box = SciMLBase.symbolic_discretize(pdesys, MockBoxDiscretization(t); checks = false)
    @test isequal(box[1].equations, generic[1].equations)
    @test isequal(box[2].eqvar, generic[2].eqvar)
    @test box[3] === generic[3]
    @test box[4] == generic[4]
    @test box[6] === generic[6]
end

@testset "lower_boundaries hook" begin
    pdesys, t = heat_pdesys()
    @test PDEBase.lower_boundaries(pdesys, PlainDiscretization()) === pdesys
    @test PDEBase.lower_boundaries(pdesys, MockDiscretization(t)) === pdesys
    @test PDEBase.lower_boundaries(pdesys, MockBoxDiscretization(t)) === pdesys

    original, t = heat_pdesys()
    lowered, _ = heat_pdesys(coefficient = 2)
    disc = LoweringDiscretization(t, lowered, Ref{Any}(nothing))
    state, = SciMLBase.symbolic_discretize(original, disc; checks = false)
    @test disc.received[] === original

    direct_system, t = heat_pdesys(coefficient = 2)
    direct, = SciMLBase.symbolic_discretize(direct_system, MockDiscretization(t); checks = false)
    unlowered_system, t = heat_pdesys()
    unlowered, = SciMLBase.symbolic_discretize(unlowered_system, MockDiscretization(t); checks = false)
    @test isequal(state.equations, direct.equations)
    @test !isequal(state.equations, unlowered.equations)
end

@testset "split_ics_bcs hook" begin
    pdesys, t = heat_pdesys()
    disc = MockDiscretization(t)
    v = VariableMap(pdesys, disc)
    bcorders = Dict(map(x -> x => d_orders(x, PDEBase.get_bcs(pdesys)), all_ivs(v)))
    boundarymap = PDEBase.parse_bcs(PDEBase.get_bcs(pdesys), v, bcorders)

    expected_ics = mapreduce(u -> boundarymap[u][t], vcat, operation.(depvars(v)))
    expected_bcmap = Dict(
        map(collect(keys(boundarymap))) do u
            u => Dict(map(x -> x => boundarymap[u][x], indvars(v)))
        end
    )

    ics, bcmap = PDEBase.split_ics_bcs(boundarymap, v, disc)
    @test length(ics) == 1
    @test isequal(ics, expected_ics)
    @test isequal(collect(keys(bcmap)), collect(keys(expected_bcmap)))
    for u in keys(expected_bcmap)
        @test isequal(collect(keys(bcmap[u])), collect(keys(expected_bcmap[u])))
        for x in keys(expected_bcmap[u])
            @test isequal(bcmap[u][x], expected_bcmap[u][x])
        end
    end

    override_ics, override_bcmap = PDEBase.split_ics_bcs(
        boundarymap, v, LoweringDiscretization(t, nothing, Ref{Any}(nothing))
    )
    @test isequal(override_ics, expected_ics)
    for u in keys(expected_bcmap), x in keys(expected_bcmap[u])
        @test isequal(override_bcmap[u][x], expected_bcmap[u][x])
    end

    steady_ics, steady_bcmap = PDEBase.split_ics_bcs(boundarymap, v, PlainDiscretization())
    @test isempty(steady_ics)
    @test isequal(collect(keys(steady_bcmap)), collect(keys(expected_bcmap)))
end
