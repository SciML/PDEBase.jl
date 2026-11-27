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
# Fresh-process costs. Each sample starts a new Julia process, so the precompiled workload
# (rather than a warm session) decides the result. "load" stops after loading the packages
# and building the system, and "first_calls" also runs the first expand_vector_operators,
# parsed_boundary_map and StrongForm.analyse calls. Both stages report whole-process wall
# time. To measure the three first calls in isolation, time them inside the child after
# package loading and model construction.
#
# "import_deps_loaded" times `using PDEBase` in a fresh process whose direct dependencies are
# already loaded, which is the cost every session that has the dependencies loaded pays for
# the precompiled workload. The child loads the dependencies and reports `ready` in the untimed
# setup, and the sample covers only the round trip that triggers `using PDEBase`.
# =============================================================================

const TTFX_CHILD = joinpath(@__DIR__, "ttfx_child.jl")
ttfx_cmd(stage) = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $TTFX_CHILD $stage`

SUITE["ttfx"] = BenchmarkGroup()

for stage in ("load", "first_calls")
    SUITE["ttfx"][stage] = @benchmarkable run($(ttfx_cmd(stage))) samples = 5 evals = 1 seconds = 600
end

SUITE["ttfx"]["import_deps_loaded"] = @benchmarkable begin
    println(child, "go")
    readline(child)
end setup = (child = open($(ttfx_cmd("import")), "r+"); readline(child)) teardown = wait(child) samples = 5 evals = 1 seconds = 600
