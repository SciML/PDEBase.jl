# Symbolic golden-outcome gate (spec §10.5) over every system of PDESystemLibrary 0.1.10:
# accepted or rejected-with-message per system, and for each accepted system the
# reconstruction identity, linearity, BC-kind split and absence of warnings.
using PDEBase, Symbolics, Test
using PDEBase: VariableMap, all_ivs, d_orders, parse_bcs, indvars, LowerBoundary,
    UpperBoundary, InterfaceBoundary, isupper
using PDEBase.StrongForm: StrongFormAnalysis, Weighted, analyse, declared_fields
using Base.CoreLogging: Warn, with_logger
using Symbolics: Differential, unwrap
using SymbolicUtils: iscall, operation, arguments
import PDESystemLibrary

include(joinpath(@__DIR__, "shared", "strongform_identity.jl"))

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
# Systems whose default ProductRule analysis warns about a non-symmetric correction term. sph,
# 4/r² ∂ᵣ(r² ∂ᵣu), is the spec's own example (§5.1), whose remedy is `Weighted(r²)`.
const PRODUCT_RULE_WARNS = Set([:sph])

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

"""
A library system preprocessed in `symbolic_discretize`'s order, with its `VariableMap` and
`parse_bcs` boundary map. `handle_complex` returns a copy with fresh equation vectors, so
`cardinalize_eqs!` never mutates the shared library system.
"""
function library_setup(sys)
    sys, _ = PDEBase.handle_complex(sys)
    PDEBase.cardinalize_eqs!(sys)
    sys, replaced = PDEBase.make_pdesys_compatible(sys)
    it = findfirst(iv -> Symbol(iv) == :t, sys.ivs)
    v = VariableMap(sys, LibraryTime(it === nothing ? nothing : sys.ivs[it]); replaced_vars = replaced)
    orders = Dict(x => d_orders(x, sys.bcs) for x in all_ivs(v))
    return sys, v, parse_bcs(sys.bcs, v, orders)
end

# Parameter slots: a Pair stands for its key, an array parameter for its elements column-major.
function params(v)
    ps = Any[]
    for p in v.ps
        p = unwrap(p isa Pair ? p.first : p)
        Symbolics.symtype(p) <: AbstractArray ? append!(ps, unwrap.(vec(collect(Symbolics.scalarize(p))))) : push!(ps, p)
    end
    return ps
end
# Field f is the f-th dependent variable in the order the system declares them.
fieldops(sys, v) = Any[operation(unwrap(u)) for u in declared_fields(sys, v)]
hasderivative(ex) = (ex = unwrap(ex); iscall(ex) && (operation(ex) isa Differential || any(hasderivative, arguments(ex))))
depends(ex, syms) = any(s -> any(isequal(unwrap(s)), Symbolics.get_variables(ex)), syms)

hessian_index(f, d, e, dim) = ((f - 1) * dim + d - 1) * dim + e

# Slot of the state atom `ex`: a field call under at most two Differentials.
function atom_slot(ex, ops, v, slots)
    ds = Any[]
    while iscall(ex) && operation(ex) isa Differential
        append!(ds, fill(operation(ex).x, operation(ex).order))
        ex = unwrap(only(arguments(ex)))
    end
    f = findfirst(op -> isequal(op, operation(ex)), ops)
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
function reference_slots(eq, sys, v, slots)
    e = Symbolics.expand_derivatives(eq.lhs - eq.rhs)
    ops = fieldops(sys, v)
    rules = Dict{Any, Any}(a => atom_slot(a, ops, v, slots) for a in collect_atoms!(Set{Any}(), e, ops))
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
function analyse_logged(sys, v, bmap; ibp = Dict())
    logger = Test.TestLogger(; min_level = Warn, respect_maxlog = false)
    result = with_logger(logger) do
        try
            analyse(sys, v, bmap; ibp)
        catch err
            err
        end
    end
    return result, logger.logs
end

const SYSTEMS = PDESystemLibrary.all_systems
# Snapshot before any preprocessing, to show that the gate never mutates the library.
const ORIGINAL = [(copy(sys.eqs), copy(sys.bcs)) for sys in SYSTEMS]

@testset "corpus" begin
    @test length(SYSTEMS) == 427
    @test allunique(getproperty.(SYSTEMS, :name))
    @test issubset(keys(REJECTED), getproperty.(SYSTEMS, :name))
    @test issubset(PERIODIC ∪ NONLINEAR, getproperty.(SYSTEMS, :name))
    # L6: trans_sin is defined but never added to all_systems
    @test_broken any(s -> s.name === :trans_sin, SYSTEMS)
    # The Weighted(r²) remedy: the weighted conormal vanishes on r = 0, so Dr(u(t, 0)) ~ 0 is a
    # zero facet flux, and the Dirichlet condition at r = 1 stays essential.
    sph, v, bmap = library_setup(only(filter(s -> s.name === :sph, SYSTEMS)))
    weighted, logs = analyse_logged(sph, v, bmap; ibp = Dict(1 => Weighted(first(indvars(v))^2)))
    @test weighted isa StrongFormAnalysis && isempty(logs)
    if weighted isa StrongFormAnalysis
        @test analysed_bcs(weighted) == syntactic_bcs(v, bmap)
        @test isequal(Symbolics.simplify(only(weighted.natural).Φ), 0)
        @test only(weighted.natural).piece.side === :lo && only(weighted.essential).piece.side === :hi
    end
end

# Checks on an accepted analysis `a` of the preprocessed library system `sys` named `name`.
function check_accepted(name, sys, v, bmap, a, logs)
    slots = a.slots
    if name in PRODUCT_RULE_WARNS
        @test length(logs) == 1 && occursin("product-rule correction terms", only(logs).message)
    else
        @test isempty(logs)
    end
    @test isempty(a.periodic) == !(name in PERIODIC)
    @test (a.nfields, a.dim) == (length(v.ū), length(indvars(v)))
    @test length(a.strong) == length(sys.eqs)
    @test !any(G -> depends(G, vcat(slots.H, slots.Ud)), reduce(vcat, a.flux))
    rng = Draws(0)
    bases = [random_state(rng, slots, v) for _ in 1:NSTATES]
    for (k, eq) in enumerate(sys.eqs)
        reference = reference_slots(eq, sys, v, slots)
        divergence = sum(total_derivative(G, i, slots, a.dim) for (i, G) in enumerate(a.flux[k]); init = Num(0))
        reconstructed = a.mass[k] + divergence + a.remainder[k]
        for vals in bases
            e = evaluate(a.strong[k], vals)
            @test isapprox(e, evaluate(reference, vals); rtol = RTOL, atol = RTOL)
            @test isapprox(e, evaluate(reconstructed, vals); rtol = RTOL, atol = RTOL)
        end
    end
    linear = isaffine(vcat(a.mass, reduce(vcat, a.flux), a.remainder), slots, rng, bases[1])
    @test linear == !(name in NONLINEAR)
    @test analysed_bcs(a) == syntactic_bcs(v, bmap)
    return nothing
end

outcomes = Dict{Symbol, String}()
@testset "golden outcome: $(sys.name)" for sys in SYSTEMS
    prepared, v, bmap = library_setup(sys)
    result, logs = analyse_logged(prepared, v, bmap)
    accepted = result isa StrongFormAnalysis
    outcomes[sys.name] = accepted ? "accepted" : sprint(showerror, result)
    if haskey(REJECTED, sys.name)
        @test_broken accepted
        @test result isa ErrorException && occursin(REJECTED[sys.name], result.msg)
    else
        @test accepted
        accepted && check_accepted(sys.name, prepared, v, bmap, result, logs)
    end
end

@test all(isequal((sys.eqs, sys.bcs), o) for (sys, o) in zip(SYSTEMS, ORIGINAL))

naccepted = count(==("accepted"), values(outcomes))
println("StrongForm golden outcomes: $naccepted accepted, $(length(outcomes) - naccepted) rejected")
for name in sort!(collect(keys(outcomes)))
    outcomes[name] == "accepted" || println("  rejected ", name, ": ", first(outcomes[name], 200))
end
