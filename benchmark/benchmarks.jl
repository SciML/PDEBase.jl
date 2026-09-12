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
