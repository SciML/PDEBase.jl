# Fresh-process workloads, run as a subprocess by benchmark/benchmarks.jl. `load` stops after
# building the system. `first_calls` also runs the main entry points once and prints an
# `analysis nfields=N dim=D` record from the analysis result. `import` loads the direct
# dependencies of PDEBase, prints `ready`, waits for a `go` line on stdin and prints the seconds
# that `using PDEBase` took; any other input exits unsuccessfully without a duration.
const STAGES = ("load", "first_calls", "import")
const STAGE = only(ARGS)
STAGE in STAGES ||
    throw(ArgumentError("unknown stage $(repr(STAGE)); the allowed stages are $(join(STAGES, ", "))"))

if STAGE == "import"
    import Pkg, TOML
    pdebase_uuid = Base.UUID("a7812802-0625-4b9e-961c-d332478797e5")
    active = TOML.parsefile(Base.active_project())
    # PDEBase is the active project itself in a checkout, and a dependency of it otherwise.
    pdebase_dir = get(active, "uuid", "") == string(pdebase_uuid) ? dirname(Base.active_project()) :
        Pkg.dependencies()[pdebase_uuid].source
    pdebase_project = joinpath(pdebase_dir, "Project.toml")
    # The project of PDEBase goes first on the load path so that its direct dependencies resolve
    # by name, and the stack is restored before the timed `using PDEBase`.
    pushfirst!(LOAD_PATH, pdebase_project)
    try
        for (name, uuid) in TOML.parsefile(pdebase_project)["deps"]
            Core.eval(Main, :(import $(Symbol(name))))
            # A checkout with its own manifest could resolve a different version than the active
            # environment, which the timed import would then reuse.
            expected = Pkg.dependencies()[Base.UUID(uuid)].version
            loaded = pkgversion(Base.invokelatest(getfield, Main, Symbol(name)))
            if expected !== nothing && loaded !== nothing && loaded != expected
                println(stderr, "$name $loaded was preloaded, but the active environment resolves $expected")
                exit(1)
            end
        end
    finally
        popfirst!(LOAD_PATH)
    end
    println("ready")
    flush(stdout)
    if readline() != "go"
        println(stderr, "expected `go` on stdin")
        exit(1)
    end
    seconds = @elapsed Core.eval(Main, :(using PDEBase))
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
    analysis = PDEBase.StrongForm.analyse(parsed.pdesys, parsed.v, parsed.boundarymap)
    println("analysis nfields=", analysis.nfields, " dim=", analysis.dim)
end
