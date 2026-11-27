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
