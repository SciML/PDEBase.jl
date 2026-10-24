"""
    SlotContext(v::VariableMap; fields = depvars(v), params = v.ps)

How the symbols of a `PDESystem` map to slots: field `f` is `fields[f]`, coordinate `X[d]` is
`indvars(v)[d]`, `T` is the time and `P` holds the parameters in the order of `params`, each
array parameter contributing its scalar elements in column-major order (`a[1], a[2], …`).
`Pair` entries `p => default` stand for `p`.
"""
struct SlotContext
    ops::Vector{Any}
    args::Any
    sivs::Vector{Any}
    t::Any
    params::Vector{Any}
    intervals::Any
end
function SlotContext(v::VariableMap; fields = depvars(v), params = v.ps)
    ps = Any[]
    for p in params
        p = unwrap(p isa Pair ? p.first : p)
        symtype(p) <: AbstractArray ? append!(ps, unwrap.(vec(collect(Symbolics.scalarize(p))))) :
            push!(ps, p)
    end
    ops = Any[operation(unwrap(u)) for u in fields]
    return SlotContext(ops, v.args, collect(Any, indvars(v)), v.time, ps, v.intervals)
end

_field(ctx::SlotContext, ex) = iscall(ex) ? findfirst(op -> isequal(op, operation(ex)), ctx.ops) : nothing

"Field index of `ex` if it is a field called at its declared (interior) arguments."
function _interior_field(ctx::SlotContext, ex)
    f = _field(ctx, ex)
    f === nothing && return nothing
    return isequal(collect(arguments(ex)), collect(ctx.args[ctx.ops[f]])) ? f : nothing
end

function _slot_leaf(ctx::SlotContext, ex, slots)
    d = findfirst(x -> isequal(x, ex), ctx.sivs)
    d === nothing || return unwrap(slots.X[d])
    ctx.t !== nothing && isequal(ex, ctx.t) && return unwrap(slots.T)
    j = findfirst(p -> isequal(p, ex), ctx.params)
    j === nothing || return unwrap(slots.P[j])
    return ex
end

hessian_index(f, d, e, dim) = ((f - 1) * dim + d - 1) * dim + e

"Slot of a derivative of an interior field: `Ud`, `G` or `H`; `nothing` for anything else."
function _derivative_slot(ctx::SlotContext, ex, slots)
    op, inner, dim = operation(ex), only(arguments(ex)), length(ctx.sivs)
    e = findfirst(isequal(op.x), ctx.sivs)
    f = _interior_field(ctx, inner)
    if f !== nothing
        ctx.t !== nothing && isequal(op.x, ctx.t) && op.order == 1 && return unwrap(slots.Ud[f])
        e === nothing && return nothing
        op.order == 1 && return unwrap(slots.G[(f - 1) * dim + e])
        op.order == 2 && return unwrap(slots.H[hessian_index(f, e, e, dim)])
        return nothing
    end
    (e !== nothing && op.order == 1 && iscall(inner) && operation(inner) isa Differential) ||
        return nothing
    d = findfirst(isequal(operation(inner).x), ctx.sivs)
    f = _interior_field(ctx, only(arguments(inner)))
    (d === nothing || f === nothing || operation(inner).order != 1) && return nothing
    return unwrap(slots.H[hessian_index(f, d, e, dim)])
end

function _slot_interior(ctx::SlotContext, ex, slots)
    iscall(ex) || return _slot_leaf(ctx, ex, slots)
    j = findfirst(p -> isequal(p, ex), ctx.params)
    j === nothing || return unwrap(slots.P[j])
    op = operation(ex)
    if op isa Differential
        slot = _derivative_slot(ctx, ex, slots)
        slot === nothing || return slot
    else
        f = _interior_field(ctx, ex)
        f === nothing || return unwrap(slots.U[f])
    end
    args = map(a -> _slot_interior(ctx, a, slots), arguments(ex))
    return maketerm(typeof(ex), op, args, metadata(ex))
end

function _leftovers!(found, ex, slotset)
    if iscall(ex)
        op = operation(ex)
        (op isa Differential || op isa SymbolicUtils.BasicSymbolic) && return push!(found, ex)
        foreach(a -> _leftovers!(found, a, slotset), arguments(ex))
    elseif SymbolicUtils.issym(ex) && !(ex in slotset)
        push!(found, ex)
    end
    return found
end

"""
    slot_volume(ctx::SlotContext, ex, slots::StateSlots) -> Num

`ex` from a volume equation in slot form: derivatives are expanded, then interior fields become
`U`, their derivatives `Ud`, `G` and `H`, coordinates `X`, time `T` and parameters `P`. Errors
on any symbol left over, such as an unregistered coefficient or its derivative. Only for volume
equations: a boundary equation must not be derivative-expanded while it holds traces.
"""
function slot_volume(ctx::SlotContext, ex, slots::StateSlots)
    ex = unwrap(Symbolics.expand_derivatives(Symbolics.wrap(ex)))
    out = _slot_interior(ctx, ex, slots)
    slotset = Set(unwrap.(vcat(slots.U, slots.G, slots.Ud, slots.H, slots.X, [slots.T], slots.P)))
    leftover = unique(_leftovers!(Any[], out, slotset))
    isempty(leftover) ||
        error("unresolved symbols $(Num.(leftover)) in $(Num(out)); register the derivative")
    return Num(out)
end
