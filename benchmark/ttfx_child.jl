# Fresh-process time-to-first-call workload, run as a subprocess by benchmark/benchmarks.jl.
# `load` stops after building the system, `first_calls` also runs the main entry points once.
using PDEBase, ModelingToolkit, Symbolics

struct TTFXDiscretization <: PDEBase.AbstractBoxDiscretization
    time::Any
end
PDEBase.get_time(disc::TTFXDiscretization) = disc.time

@parameters t x y
@variables u(..)
Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
U = u(t, x, y)
pdesys = PDESystem(
    [Dt(U) ~ Dx(Dx(U)) + Dy(Dy(U))],
    [u(0, x, y) ~ 0, u(t, 0, y) ~ 0, u(t, 1, y) ~ 0, u(t, x, 0) ~ 0, u(t, x, 1) ~ 0],
    [t ∈ (0, 1), x ∈ (0, 1), y ∈ (0, 1)], [t, x, y], [U]; name = :heat_2d
)

if only(ARGS) == "first_calls"
    PDEBase.expand_vector_operators(pdesys)
    parsed = PDEBase.parsed_boundary_map(pdesys, TTFXDiscretization(t))
    PDEBase.StrongForm.analyse(parsed.pdesys, parsed.v, parsed.boundarymap)
end
