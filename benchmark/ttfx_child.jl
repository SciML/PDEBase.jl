# Fresh-process workloads, run as a subprocess by benchmark/benchmarks.jl. `load` stops after
# building the system, `first_calls` also runs the main entry points once, and `import` loads
# the direct dependencies of PDEBase, then prints `ready`, waits for a line on stdin and times `using PDEBase`.
const STAGES = ("load", "first_calls", "import")
const STAGE = only(ARGS)
STAGE in STAGES ||
    throw(ArgumentError("unknown stage $(repr(STAGE)); the allowed stages are $(join(STAGES, ", "))"))

if STAGE == "import"
    pdebase = Base.PkgId(Base.UUID("a7812802-0625-4b9e-961c-d332478797e5"), "PDEBase")
    project = Base.parsed_toml(joinpath(dirname(Base.locate_package(pdebase)), "..", "Project.toml"))
    for (name, uuid) in project["deps"]
        Base.require(Base.PkgId(Base.UUID(uuid), name))
    end
    println("ready")
    flush(stdout)
    readline()
    seconds = @elapsed Base.require(pdebase)
    println(seconds)
    flush(stdout)
    exit()
end

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

if STAGE == "first_calls"
    PDEBase.expand_vector_operators(pdesys)
    parsed = PDEBase.parsed_boundary_map(pdesys, TTFXDiscretization(t))
    PDEBase.StrongForm.analyse(parsed.pdesys, parsed.v, parsed.boundarymap)
end
