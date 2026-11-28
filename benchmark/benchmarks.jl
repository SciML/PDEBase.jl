using PDEBase, BenchmarkTools
using ModelingToolkit, Symbolics, SciMLBase

const SUITE = BenchmarkGroup()

# =============================================================================
# VariableMap construction — symbolic variable ↔ equation mapping
# =============================================================================

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
    name = :heat_1d
)

SUITE["varmap"] = BenchmarkGroup()

SUITE["varmap"]["construct"] = @benchmarkable VariableMap($pdesys)

# Small maps, where the sort key of the dependent calls is computed at most a couple of times:
# a steady mesh system with a Dirichlet trace and a transient one with a trace and an initial
# condition, both on a (y, x) mesh pairing.
struct VarMapBenchMesh end
PDEBase.meshdim(::VarMapBenchMesh) = 2
PDEBase.facetset_names(::VarMapBenchMesh) = (:inlet, :wall)
PDEBase.cellset_names(::VarMapBenchMesh) = (:all,)
PDEBase.contains_point(::VarMapBenchMesh, p) = all(c -> 0 <= c <= 1, p)

struct VarMapBenchDiscretization <: PDEBase.AbstractMeshDiscretization
    time::Any
end
PDEBase.get_time(disc::VarMapBenchDiscretization) = disc.time

@parameters y k
mesh = MeshDomain(VarMapBenchMesh())
inlet = boundarypiece(mesh, :inlet)
Dy = Differential(y)

steady_mesh = PDESystem(
    [0 ~ k * (Dx(Dx(u(x, y))) + Dy(Dy(u(x, y))))], [u(inlet) ~ 1], [(y, x) ∈ mesh],
    [x, y], [u(x, y)], [k]; name = :steady_mesh_trace
)
transient_mesh = PDESystem(
    [Dt(u(t, x, y)) ~ k * (Dx(Dx(u(t, x, y))) + Dy(Dy(u(t, x, y))))],
    [u(t, inlet) ~ sin(t), u(0, x, y) ~ 0], [(y, x) ∈ mesh, t ∈ (0, 1)],
    [t, x, y], [u(t, x, y)], [k]; name = :transient_mesh_trace_ic
)

SUITE["varmap"]["construct_steady_mesh_trace"] = @benchmarkable VariableMap(
    $steady_mesh, VarMapBenchDiscretization(nothing)
)
SUITE["varmap"]["construct_transient_mesh_trace_ic"] = @benchmarkable VariableMap(
    $transient_mesh, VarMapBenchDiscretization($t)
)

v = VariableMap(pdesys)

SUITE["varmap"]["x2i"] = @benchmarkable x2i(
    $v, $(Symbolics.unwrap(u(t, x))), $(Symbolics.unwrap(x))
)
SUITE["varmap"]["depvar_ops"] = @benchmarkable getproperty($v, :depvar_ops)
SUITE["varmap"]["count_differentials"] = @benchmarkable count_differentials(
    $(Symbolics.unwrap(Dx(Dx(u(t, x))))), $(Symbolics.unwrap(x))
)
SUITE["varmap"]["differential_order"] = @benchmarkable differential_order(
    $eq, $(Symbolics.unwrap(x))
)
SUITE["varmap"]["find_derivative"] = @benchmarkable find_derivative(
    $(Symbolics.unwrap(Dx(u(t, x)) + 1)),
    $(operation(Symbolics.unwrap(u(t, x))))
)

# =============================================================================
# Interface errors / boundary checks on the default discretization path
# =============================================================================

struct PlainDiscretization <: SciMLBase.AbstractDiscretization end

SUITE["interface"] = BenchmarkGroup()

SUITE["interface"]["interface_errors"] = @benchmarkable PDEBase.interface_errors(
    $pdesys, $v, PlainDiscretization()
)
SUITE["interface"]["check_boundarymap"] = @benchmarkable PDEBase.check_boundarymap(
    Dict(), $v, PlainDiscretization()
)
SUITE["interface"]["should_transform"] = @benchmarkable PDEBase.should_transform(
    $pdesys, PlainDiscretization(), Dict()
)

# =============================================================================
# Fresh-process regression guards for the cost of loading PDEBase and of its first calls.
# Each sample starts a new Julia process, so what is measured is the precompiled state of the
# package rather than a warm session.
#   "cold_import"  whole-process wall time of loading PDEBase and exiting, without building
#                  or analysing a model.
#   "load"         whole-process wall time of also loading ModelingToolkit and Symbolics and
#                  building the system.
#   "first_calls"  whole-process wall time of also running the first parsed_boundary_map and
#                  StrongForm.analyse calls. To collect those two first calls in isolation,
#                  timed inside the child after package loading and model construction, run
#                  `julia --project=<environment> benchmark/ttfx_child.jl first_calls` and read
#                  its `first_calls parsed_boundary_map=S analyse=S` line.
#   "import_deps_loaded"  `using PDEBase` in a fresh process whose direct dependencies are
#                  already loaded, which is the cost a session that has those dependencies
#                  loaded pays for loading PDEBase. The child loads the dependencies and
#                  reports `ready` in the untimed setup, and the sample covers only the round
#                  trip that triggers `using PDEBase`.
# =============================================================================

const TTFX_CHILD = joinpath(@__DIR__, "ttfx_child.jl")
ttfx_cmd(stage) = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $TTFX_CHILD $stage`

SUITE["ttfx"] = BenchmarkGroup()

# The children print records that only the tests and a direct invocation read.
for stage in ("cold_import", "load", "first_calls")
    SUITE["ttfx"][stage] = @benchmarkable run(pipeline($(ttfx_cmd(stage)); stdout = devnull)) samples = 5 evals = 1 seconds = 600
end

SUITE["ttfx"]["import_deps_loaded"] = @benchmarkable begin
    println(child, "go")
    name, seconds = split(readline(child))
    name == "PDEBase" || error("the import stage of ttfx_child.jl reported $name")
    parse(Float64, seconds)
end setup = (
    child = open($(ttfx_cmd("import")), "r+"); readline(child) == "ready" || error("the import stage of ttfx_child.jl did not report ready")
) teardown = (wait(child); success(child) || error("the import stage of ttfx_child.jl failed")) samples = 5 evals = 1 seconds = 600
