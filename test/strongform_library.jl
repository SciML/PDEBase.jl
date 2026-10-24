# Symbolic golden-outcome gate (spec §10.5) over every system of PDESystemLibrary 0.1.10:
# accepted or rejected-with-message per system, and for each accepted system the
# reconstruction identity, linearity, BC-kind split and absence of warnings.
using PDEBase, Symbolics, Test
using PDEBase: VariableMap, all_ivs, d_orders, parse_bcs, indvars, LowerBoundary,
    UpperBoundary, InterfaceBoundary, isupper
using PDEBase.StrongForm: StrongFormAnalysis, analyse
using Base.CoreLogging: Warn, with_logger
using Symbolics: Differential, unwrap
using SymbolicUtils: iscall, operation, arguments
import PDESystemLibrary

const THIRD = "third-order operator needs a mixed or C¹ formulation"
const FOURTH = "fourth-order operator needs a mixed or C¹ formulation"

# Expected outcomes. Every system not listed here is accepted. The 7 higher-order systems are
# rejected for a feature with no milestone (spec §10.5), so their acceptance is @test_broken.
const REJECTED = Dict(
    :adv3 => THIRD, :advdiff3 => THIRD, :gen1 => THIRD, :gen2 => THIRD, :gen3 => THIRD,
    :diff4 => FOURTH, :advdiffno3 => FOURTH,
)
# Periodic-only systems: the analysis accepts them with `PeriodicBC`s. Only the Ferrite backend
# rejects them before M4, so they are not @test_broken here.
const PERIODIC = Set(
    [:convsin, :convcos, :convtri, :convsquare, :convsquare2, :adv, :diff, :advdiff, :bruss]
)
const NONLINEAR = Set([:inviscid_burgers_monotonic, :burgers_2d, :bruss])

const NSTATES = 5
const RTOL = 1.0e-9

# Deterministic draws in (0, 1), a Weyl sequence, so that the gate needs no Random dependency.
mutable struct Draws
    k::Int
end
Base.rand(d::Draws) = (d.k += 1; mod(0.5 + d.k * 0.6180339887498949, 1.0))

struct LibraryTime <: PDEBase.AbstractEquationSystemDiscretization
    t::Any
end
PDEBase.get_time(d::LibraryTime) = d.t

"The `VariableMap` and `parse_bcs` boundary map of a library system, as `symbolic_discretize` builds them."
function library_setup(sys)
    it = findfirst(iv -> Symbol(iv) == :t, sys.ivs)
    v = VariableMap(sys, LibraryTime(it === nothing ? nothing : sys.ivs[it]))
    orders = Dict(x => d_orders(x, sys.bcs) for x in all_ivs(v))
    return v, parse_bcs(sys.bcs, v, orders)
end

params(v) = Any[unwrap(p isa Pair ? p.first : p) for p in v.ps]
fieldops(v) = Any[operation(unwrap(u)) for u in v.ū]
hasderivative(ex) = (ex = unwrap(ex); iscall(ex) && (operation(ex) isa Differential || any(hasderivative, arguments(ex))))
depends(ex, syms) = any(s -> any(isequal(unwrap(s)), Symbolics.get_variables(ex)), syms)

hessian_index(f, d, e, dim) = ((f - 1) * dim + d - 1) * dim + e

# Slot of the state atom `ex`: a field call under at most two Differentials.
function atom_slot(ex, v, slots)
    ds = Any[]
    while iscall(ex) && operation(ex) isa Differential
        append!(ds, fill(operation(ex).x, operation(ex).order))
        ex = unwrap(only(arguments(ex)))
    end
    f = findfirst(op -> isequal(op, operation(ex)), fieldops(v))
    sivs, dim = indvars(v), length(indvars(v))
    isempty(ds) && return slots.U[f]
    length(ds) == 1 && isequal(ds[1], v.time) && return slots.Ud[f]
    dirs = [findfirst(x -> isequal(x, d), sivs) for d in ds]
    any(isnothing, dirs) && error("mixed time-space derivative $ex has no slot")
    length(dirs) == 1 && return slots.G[(f - 1) * dim + dirs[1]]
    length(dirs) == 2 && return slots.H[hessian_index(f, dirs[1], dirs[2], dim)]
    return error("derivative of order $(length(dirs)) in $ex has no slot")
end

function collect_atoms!(acc, ex, ops)
    ex = unwrap(ex)
    iscall(ex) || return acc
    op = operation(ex)
    if op isa Differential || any(o -> isequal(o, op), ops)
        push!(acc, ex)
    else
        foreach(a -> collect_atoms!(acc, a, ops), arguments(ex))
    end
    return acc
end

"""
The residual `lhs - rhs` of equation `eq` in slot form, built independently of `analyse`:
derivatives expanded, then every state atom, coordinate, time and parameter replaced by its slot.
"""
function reference_slots(eq, v, slots)
    e = Symbolics.expand_derivatives(eq.lhs - eq.rhs)
    rules = Dict{Any, Any}(a => atom_slot(a, v, slots) for a in collect_atoms!(Set{Any}(), e, fieldops(v)))
    for (d, x) in enumerate(indvars(v))
        rules[x] = slots.X[d]
    end
    v.time === nothing || (rules[v.time] = slots.T)
    for (j, p) in enumerate(params(v))
        rules[p] = slots.P[j]
    end
    return substitute(e, rules; fold = Val(false))
end

"Random slot values: coordinates and time inside the domain, `H` symmetric per field."
function random_state(rng, slots, v)
    vals = Dict{Num, Float64}()
    draw() = 0.5 + rand(rng)
    inside(x) = ((lo, hi) = v.intervals[x]; lo + (hi - lo) * (0.1 + 0.8 * rand(rng)))
    for (d, x) in enumerate(indvars(v))
        vals[slots.X[d]] = inside(x)
    end
    vals[slots.T] = v.time === nothing ? draw() : inside(v.time)
    for s in vcat(slots.U, slots.Ud, slots.G, slots.N, slots.P)
        vals[s] = draw()
    end
    dim = length(slots.X)
    for f in eachindex(slots.U), d in 1:dim, e in d:dim
        vals[slots.H[hessian_index(f, d, e, dim)]] = vals[slots.H[hessian_index(f, e, d, dim)]] = draw()
    end
    return vals
end

evaluate(ex, vals) = Float64(Symbolics.value(substitute(Num(ex), vals; fold = Val(true))))

# ∂ᵢ of a slot expression holding at most first derivatives: the total derivative through
# the state slots, ∂/∂Xᵢ + Σ_f ∂/∂U_f ⋅ G_fi + Σ_fd ∂/∂G_fd ⋅ H_fdi.
function total_derivative(ex, slots, i)
    dim = length(slots.X)
    D = Symbolics.derivative(ex, slots.X[i])
    for f in eachindex(slots.U)
        D += Symbolics.derivative(ex, slots.U[f]) * slots.G[(f - 1) * dim + i]
        for d in 1:dim
            D += Symbolics.derivative(ex, slots.G[(f - 1) * dim + d]) * slots.H[hessian_index(f, d, i, dim)]
        end
    end
    return D
end

statevars(slots) = vcat(slots.U, slots.G, slots.Ud, slots.H)

# Affine in the state slots iff f(a) + f(b) = f(a + b) + f(0) at random state values.
function isaffine(exs, slots, rng, base)
    state = statevars(slots)
    for _ in 1:3
        a = Dict(s => rand(rng) for s in state)
        b = Dict(s => rand(rng) for s in state)
        at(S) = evaluate.(exs, Ref(merge(base, S)))
        fa, fb = at(a), at(b)
        r = fa .+ fb .- at(Dict(s => a[s] + b[s] for s in state)) .- at(Dict(s => 0.0 for s in state))
        all(abs.(r) .<= RTOL .* (1 .+ abs.(fa) .+ abs.(fb))) || return false
    end
    return true
end

# Spatial conditions by syntactic kind: (iv, side, :essential | :natural), and the periodic count.
function syntactic_bcs(v, bmap)
    kinds = Tuple{String, Symbol, Symbol}[]
    nperiodic = 0
    for (_, bykey) in bmap, (x, bs) in bykey, b in bs
        v.time !== nothing && isequal(x, v.time) && continue
        if b isa InterfaceBoundary
            isupper(b) || (nperiodic += 1)
        elseif b isa Union{LowerBoundary, UpperBoundary}
            kind = hasderivative(b.eq.lhs - b.eq.rhs) ? :natural : :essential
            push!(kinds, (string(b.x), b isa LowerBoundary ? :lo : :hi, kind))
        end
    end
    return sort(kinds), nperiodic
end

function analysed_bcs(a::StrongFormAnalysis)
    kinds = vcat(
        [(string(bc.piece.iv), bc.piece.side, :essential) for bc in a.essential],
        [(string(bc.piece.iv), bc.piece.side, :natural) for bc in a.natural]
    )
    return sort(kinds), length(a.periodic)
end

# Run `analyse`, returning the analysis or the exception, and every warning it logged.
function analyse_logged(sys, v, bmap)
    logger = Test.TestLogger(; min_level = Warn, respect_maxlog = false)
    result = with_logger(logger) do
        try
            analyse(sys, v, bmap)
        catch err
            err
        end
    end
    return result, logger.logs
end

const SYSTEMS = PDESystemLibrary.all_systems

@testset "corpus" begin
    @test length(SYSTEMS) == 427
    @test allunique(getproperty.(SYSTEMS, :name))
    @test issubset(keys(REJECTED), getproperty.(SYSTEMS, :name))
    @test issubset(PERIODIC ∪ NONLINEAR, getproperty.(SYSTEMS, :name))
    # L6: trans_sin is defined but never added to all_systems
    @test_broken any(s -> s.name === :trans_sin, SYSTEMS)
end

# Checks on an accepted analysis `a` of library system `sys`.
function check_accepted(sys, v, bmap, a, logs)
    slots = a.slots
    @test isempty(logs)
    @test isempty(a.periodic) == !(sys.name in PERIODIC)
    @test (a.nfields, a.dim) == (length(v.ū), length(indvars(v)))
    @test length(a.strong) == length(sys.eqs)
    @test !any(G -> depends(G, vcat(slots.H, slots.Ud)), reduce(vcat, a.flux))
    rng = Draws(0)
    bases = [random_state(rng, slots, v) for _ in 1:NSTATES]
    for (k, eq) in enumerate(sys.eqs)
        reference = reference_slots(eq, v, slots)
        divergence = sum(total_derivative(G, slots, i) for (i, G) in enumerate(a.flux[k]); init = Num(0))
        reconstructed = a.mass[k] + divergence + a.remainder[k]
        for vals in bases
            e = evaluate(a.strong[k], vals)
            @test isapprox(e, evaluate(reference, vals); rtol = RTOL, atol = RTOL)
            @test isapprox(e, evaluate(reconstructed, vals); rtol = RTOL, atol = RTOL)
        end
    end
    linear = isaffine(vcat(a.mass, reduce(vcat, a.flux), a.remainder), slots, rng, bases[1])
    @test linear == !(sys.name in NONLINEAR)
    @test analysed_bcs(a) == syntactic_bcs(v, bmap)
    return nothing
end

outcomes = Dict{Symbol, String}()
@testset "golden outcome: $(sys.name)" for sys in SYSTEMS
    v, bmap = library_setup(sys)
    result, logs = analyse_logged(sys, v, bmap)
    accepted = result isa StrongFormAnalysis
    outcomes[sys.name] = accepted ? "accepted" : sprint(showerror, result)
    if haskey(REJECTED, sys.name)
        @test_broken accepted
        @test result isa ErrorException && occursin(REJECTED[sys.name], result.msg)
    else
        @test accepted
        accepted && check_accepted(sys, v, bmap, result, logs)
    end
end

naccepted = count(==("accepted"), values(outcomes))
println("StrongForm golden outcomes: $naccepted accepted, $(length(outcomes) - naccepted) rejected")
for name in sort!(collect(keys(outcomes)))
    outcomes[name] == "accepted" || println("  rejected ", name, ": ", first(outcomes[name], 200))
end
