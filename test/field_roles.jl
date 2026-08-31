using ModelingToolkit
using PDEBase
using SciMLBase
using SymbolicIndexingInterface
using Test

struct RoleDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
    kwargs::NamedTuple
end

struct RoleSpace <: PDEBase.AbstractDiscreteSpace
    discvars::Any
end

struct RoleMetadata <: SciMLBase.AbstractDiscretizationMetadata{true}
    pdesys::Any
    metadata::Base.RefValue{Any}
end

struct ProblemRoleDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    system::Any
    time::Any
    kwargs::NamedTuple
end

PDEBase.get_time(disc::RoleDiscretization) = disc.time
PDEBase.get_discvars(space::RoleSpace) = space.discvars
PDEBase.add_metadata!(metadata::RoleMetadata, sys) = (metadata.metadata[] = sys)
function SciMLBase.symbolic_discretize(
        ::PDESystem, disc::ProblemRoleDiscretization; checks = true
    )
    return disc.system, (0.0, 1.0)
end

contains_symbol(xs, x) = any(isequal(PDEBase.safe_unwrap(x)), PDEBase.safe_unwrap.(xs))

@testset "Active PDE field roles" begin
    @parameters t x symbolic_default
    @variables state(..)
    @variables forcing(..) [input = true]
    @variables measured(..) [output = true]
    Dt = Differential(t)
    pdesys = PDESystem(
        [Dt(state(t, x)) ~ forcing(t, x) + measured(t, x)],
        [state(0, x) ~ 0],
        [t ∈ (0, 1), x ∈ (0, 1)],
        [t, x],
        [state(t, x), forcing(t, x), measured(t, x)];
        name = :field_role_lowering
    )

    @variables state_1(t) state_2(t)
    @variables forcing_1(t) [input = true]
    @variables forcing_2(t) [input = true]
    @variables forcing_3(t) [input = true]
    @variables measured_1(t) [output = true]
    @variables measured_2(t) [output = true]
    @variables measured_3(t) [output = true]

    discvars = Dict(
        operation(PDEBase.safe_unwrap(state(t, x))) => [state_1, state_2],
        operation(PDEBase.safe_unwrap(forcing(t, x))) => [forcing_1, forcing_2, forcing_3],
        operation(PDEBase.safe_unwrap(measured(t, x))) =>
            [measured_1, measured_2, measured_3],
    )
    equations = [
        Dt(state_1) ~ forcing_3 + measured_1,
        Dt(state_2) ~ forcing_1,
        measured_1 ~ state_1,
    ]
    boundary_equations = [measured_2 ~ state_2]
    state = PDEBase.EquationState(equations, boundary_equations)
    space = RoleSpace(discvars)
    metadata = RoleMetadata(pdesys, Ref{Any}(nothing))
    u0 = [
        state_1 => 1.0,
        state_2 => 2.0,
        forcing_1 => Float32(1.25),
        forcing_2 => 100.0,
        forcing_3 => symbolic_default,
        measured_1 => 3.0,
        measured_2 => 4.0,
        measured_3 => 5.0,
    ]

    sys, tspan = PDEBase.generate_system(
        state, space, u0, (0.0, 1.0), metadata,
        RoleDiscretization(t, (;)); checks = false
    )

    @test tspan == (0.0, 1.0)
    @test isequal(ModelingToolkit.inputs(sys), PDEBase.safe_unwrap.([forcing_1, forcing_3]))
    @test isequal(ModelingToolkit.outputs(sys), PDEBase.safe_unwrap.([measured_1, measured_2]))
    @test all(x -> contains_symbol(ModelingToolkit.unknowns(sys), x), [forcing_1, forcing_3])
    @test !contains_symbol(ModelingToolkit.unknowns(sys), forcing_2)
    @test contains_symbol(ModelingToolkit.unknowns(sys), measured_3)
    @test PDEBase.SymbolicUtils.iscall(PDEBase.safe_unwrap(sys.forcing_1))

    defaults = ModelingToolkit.initial_conditions(sys)
    @test PDEBase._discrete_ic_value(defaults[forcing_1]) isa Float32
    @test isequal(PDEBase._discrete_ic_value(defaults[forcing_1]), Float32(1.25))
    @test isequal(defaults[forcing_3], symbolic_default)
    @test !haskey(defaults, forcing_2)
    @test isequal(PDEBase._discrete_ic_value(ModelingToolkit.guesses(sys)[measured_3]), 5.0)

    init_eqs = ModelingToolkit.initialization_equations(sys)
    @test all(init_eqs) do eq
        !any(x -> contains_symbol(ModelingToolkit.inputs(sys), x), Symbolics.get_variables(eq))
    end
    @test all(
        x -> !contains_symbol(keys(ModelingToolkit.guesses(sys)), x),
        ModelingToolkit.inputs(sys)
    )

    compiled = ModelingToolkit.mtkcompile(
        sys; inputs = ModelingToolkit.inputs(sys), outputs = ModelingToolkit.outputs(sys)
    )
    @test isequal(
        ModelingToolkit.inputs(compiled), PDEBase.safe_unwrap.([forcing_1, forcing_3])
    )
    @test all(
        x -> !contains_symbol(ModelingToolkit.unknowns(compiled), x),
        ModelingToolkit.inputs(compiled)
    )
    @test all(
        x -> contains_symbol(ModelingToolkit.parameters(compiled), x),
        ModelingToolkit.inputs(compiled)
    )
    @test isequal(
        ModelingToolkit.outputs(compiled), PDEBase.safe_unwrap.([measured_1, measured_2])
    )

    missing_sys, = PDEBase.generate_system(
        state, space, filter(p -> !isequal(first(p), forcing_3), u0), (0.0, 1.0), metadata,
        RoleDiscretization(t, (;)); checks = false
    )
    @test_throws "forcing_3" PDEBase._require_input_defaults(missing_sys)
    @test PDEBase._require_input_defaults(sys) === nothing
    @test_throws "forcing_3" SciMLBase.discretize(
        pdesys, ProblemRoleDiscretization(missing_sys, t, (;)); checks = false
    )
    numeric_u0 = map(u0) do pair
        isequal(first(pair), forcing_3) ? forcing_3 => 3.0 : pair
    end
    numeric_sys, = PDEBase.generate_system(
        state, space, numeric_u0, (0.0, 1.0), metadata,
        RoleDiscretization(t, (;)); checks = false
    )
    prob = SciMLBase.discretize(
        pdesys, ProblemRoleDiscretization(numeric_sys, t, (;)); checks = false
    )
    @test prob isa ODEProblem
    @test parameter_index(prob, PDEBase.safe_unwrap(forcing_1)) !== nothing
    @test parameter_index(prob, PDEBase.safe_unwrap(forcing_2)) === nothing
    get_forcing = getp(prob, PDEBase.safe_unwrap(forcing_1))
    set_forcing! = setp(prob, PDEBase.safe_unwrap(forcing_1))
    @test isequal(get_forcing(prob), Float32(1.25))
    set_forcing!(prob, 2.5f0)
    @test isequal(get_forcing(prob), Float32(2.5))
end

@testset "Slice-form role activity" begin
    @parameters t x
    @variables state(..)
    @variables measured(..) [output = true]
    Dt = Differential(t)
    pdesys = PDESystem(
        [Dt(state(t, x)) ~ measured(t, x)],
        [state(0, x) ~ 0],
        [t ∈ (0, 1), x ∈ (0, 1)],
        [t, x],
        [state(t, x), measured(t, x)];
        name = :slice_role_activity
    )

    @variables state_cells(t)[1:2]
    @variables measured_cells(t)[1:3] [output = true]
    discvars = Dict(
        operation(PDEBase.safe_unwrap(state(t, x))) => state_cells,
        operation(PDEBase.safe_unwrap(measured(t, x))) => measured_cells,
    )
    equation_state = PDEBase.EquationState(
        [Dt(measured_cells[1:2]) ~ state_cells], Equation[]
    )
    sys, = PDEBase.generate_system(
        equation_state,
        RoleSpace(discvars),
        [state_cells[1] => 1.0, state_cells[2] => 2.0],
        (0.0, 1.0),
        RoleMetadata(pdesys, Ref{Any}(nothing)),
        RoleDiscretization(t, (;));
        checks = false
    )

    @test isequal(
        ModelingToolkit.outputs(sys), PDEBase.safe_unwrap.(collect(measured_cells[1:2]))
    )
    @test contains_symbol(ModelingToolkit.unknowns(sys), measured_cells[3])
end

@testset "Derivative initial conditions with field roles" begin
    @parameters t x
    @variables state(..)
    @variables forcing(..) [input = true]
    Dt = Differential(t)
    pdesys = PDESystem(
        [Dt(Dt(state(t, x))) ~ forcing(t, x)],
        [state(0, x) ~ 1, Dt(state(0, x)) ~ 2, forcing(0, x) ~ 3],
        [t ∈ (0, 1), x ∈ (0, 1)],
        [t, x],
        [state(t, x), forcing(t, x)];
        name = :derivative_initial_conditions_with_roles
    )
    @variables state_1(t)
    @variables forcing_1(t) [input = true]
    equation_state = PDEBase.EquationState(
        [Dt(Dt(state_1)) ~ forcing_1], Equation[]
    )
    sys, = PDEBase.generate_system(
        equation_state,
        RoleSpace(
            Dict(
                operation(PDEBase.safe_unwrap(state(t, x))) => [state_1],
                operation(PDEBase.safe_unwrap(forcing(t, x))) => [forcing_1],
            )
        ),
        [state_1 => 1.0, Dt(state_1) => 2.0, forcing_1 => 3.0],
        (0.0, 1.0),
        RoleMetadata(pdesys, Ref{Any}(nothing)),
        RoleDiscretization(t, (;));
        checks = false
    )

    init_eqs = ModelingToolkit.initialization_equations(sys)
    @test any(isequal(PDEBase.safe_unwrap(state_1 ~ 1.0)), PDEBase.safe_unwrap.(init_eqs))
    @test any(
        isequal(PDEBase.safe_unwrap(Dt(state_1) ~ 2.0)), PDEBase.safe_unwrap.(init_eqs)
    )
    @test isequal(PDEBase._discrete_ic_value(ModelingToolkit.guesses(sys)[state_1]), 1.0)
    @test isequal(
        PDEBase._discrete_ic_value(ModelingToolkit.guesses(sys)[Dt(state_1)]), 2.0
    )
end

@testset "PDE role ordering and precedence" begin
    @parameters t x
    @variables state(..)
    @variables forcing_a(..) [input = true]
    @variables forcing_b(..) [input = true, output = true]
    Dt = Differential(t)
    pdesys = PDESystem(
        [Dt(state(t, x)) ~ forcing_a(t, x) + forcing_b(t, x)],
        [state(0, x) ~ 0],
        [t ∈ (0, 1), x ∈ (0, 1)],
        [t, x],
        [state(t, x), forcing_b(t, x), forcing_a(t, x)];
        name = :field_role_ordering
    )

    @variables state_1(t)
    @variables a_1(t) [input = true]
    @variables a_2(t) [input = true]
    @variables b_1(t) [input = true, output = true]
    discvars = Dict(
        operation(PDEBase.safe_unwrap(state(t, x))) => [state_1],
        operation(PDEBase.safe_unwrap(forcing_a(t, x))) => reshape([a_1, a_2], 1, 2),
        operation(PDEBase.safe_unwrap(forcing_b(t, x))) => [b_1],
    )
    state = PDEBase.EquationState(
        [Dt(state_1) ~ a_1 + a_2 + b_1], Equation[]
    )
    metadata = RoleMetadata(pdesys, Ref{Any}(nothing))
    u0 = [state_1 => 0.0, a_1 => 1.0, a_2 => 2.0, b_1 => 3.0]

    first_sys, = PDEBase.generate_system(
        state, RoleSpace(discvars), u0, (0.0, 1.0), metadata,
        RoleDiscretization(t, (;)); checks = false
    )
    second_sys, = PDEBase.generate_system(
        state, RoleSpace(discvars), u0, (0.0, 1.0), metadata,
        RoleDiscretization(t, (;)); checks = false
    )

    @test isequal(ModelingToolkit.inputs(first_sys), PDEBase.safe_unwrap.([b_1, a_1, a_2]))
    @test isempty(ModelingToolkit.outputs(first_sys))
    @test isequal(ModelingToolkit.inputs(first_sys), ModelingToolkit.inputs(second_sys))
end

@testset "No-role lowering regression" begin
    @parameters t x
    @variables state(..)
    Dt = Differential(t)
    pdesys = PDESystem(
        [Dt(state(t, x)) ~ -state(t, x)],
        [state(0, x) ~ 1],
        [t ∈ (0, 1), x ∈ (0, 1)],
        [t, x],
        [state(t, x)];
        name = :no_field_roles
    )
    @variables state_1(t) state_2(t)
    equation_state = PDEBase.EquationState(
        [Dt(state_1) ~ -state_1, Dt(state_2) ~ -state_2], Equation[]
    )
    u0 = [state_1 => 1.0, state_2 => 2.0]
    sys, = PDEBase.generate_system(
        equation_state,
        RoleSpace(
            Dict(operation(PDEBase.safe_unwrap(state(t, x))) => [state_1, state_2])
        ),
        u0,
        (0.0, 1.0),
        RoleMetadata(pdesys, Ref{Any}(nothing)),
        RoleDiscretization(t, (;));
        checks = false
    )
    legacy_init_eqs, legacy_guesses = PDEBase._discrete_initialization(
        equation_state.eqs, t, u0
    )
    legacy_sys = ModelingToolkit.System(
        equation_state.eqs, t, PDEBase.safe_unwrap.([state_1, state_2]), Num[];
        initial_conditions = Dict{Any, Any}(pdesys.initial_conditions),
        initialization_eqs = legacy_init_eqs,
        guesses = legacy_guesses,
        name = :legacy_no_field_roles,
        checks = false
    )

    @test isempty(ModelingToolkit.inputs(sys))
    @test isempty(ModelingToolkit.outputs(sys))
    @test isequal(ModelingToolkit.unknowns(sys), ModelingToolkit.unknowns(legacy_sys))
    @test isequal(ModelingToolkit.equations(sys), ModelingToolkit.equations(legacy_sys))
    @test isequal(
        ModelingToolkit.initial_conditions(sys), ModelingToolkit.initial_conditions(legacy_sys)
    )
    @test isequal(
        ModelingToolkit.initialization_equations(sys),
        ModelingToolkit.initialization_equations(legacy_sys)
    )
    @test isequal(ModelingToolkit.guesses(sys), ModelingToolkit.guesses(legacy_sys))
    compiled = ModelingToolkit.mtkcompile(sys)
    legacy_compiled = ModelingToolkit.mtkcompile(legacy_sys)
    @test isequal(
        ModelingToolkit.unknowns(compiled), ModelingToolkit.unknowns(legacy_compiled)
    )
    @test isequal(
        ModelingToolkit.parameters(compiled), ModelingToolkit.parameters(legacy_compiled)
    )
end
