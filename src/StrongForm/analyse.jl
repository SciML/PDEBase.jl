same_residual(a, b) =
    isequal(Symbolics.expand(Num(a - b)), Num(0)) || isequal(Symbolics.expand(Num(a + b)), Num(0))

same_policy(a, b) = typeof(a) == typeof(b) && (!(a isa Weighted) || isequal(a.w, b.w))

# Equation keys are matched by residual up to sign, since preprocessing such as
# `cardinalize_eqs!` rewrites `lhs ~ rhs` as `lhs - rhs ~ 0`, after flattening array variables
# the way `make_pdesys_compatible` does. A key matching nothing, or two keys giving one
# equation different policies, is an error.
function ibp_policies(ibp, eqs, v::VariableMap)
    flatten = chain_flatten_array_variables(collect(values(replaced_vars(v))))
    policies = Union{Nothing, IBPPolicy}[nothing for _ in eqs]
    residuals = [eq.lhs - eq.rhs for eq in eqs]
    for (key, policy) in ibp
        key isa Union{Integer, Equation} ||
            error("ibp keys must be Equations or integer indices, got $key")
        matched = key isa Integer ? filter(==(key), eachindex(eqs)) :
            findall(r -> same_residual(r, flatten(unwrap(key.lhs - key.rhs))), residuals)
        isempty(matched) && error("ibp key $key matches no equation of the system")
        for i in matched
            policies[i] === nothing || same_policy(policies[i], policy) ||
                error("equation $i is matched by ibp keys with different policies, $(policies[i]) and $policy")
            policies[i] = policy
        end
    end
    return IBPPolicy[something(p, ProductRule()) for p in policies]
end

"""
    analyse(pdesys, v::VariableMap, boundarymap; ibp = Dict()) -> StrongFormAnalysis

Strong-form analysis of `pdesys`, which must already be in the form `symbolic_discretize`
hands to the discretizer hooks (after `make_pdesys_compatible`), with `v` its `VariableMap`
and `boundarymap` from `parse_bcs`. Each equation `e = lhs - rhs` passes the operator whitelist
([`check_whitelist`](@ref)), is split in time ([`split_time`](@ref)) and into divergence form
([`split_divergence`](@ref)), and is written in [`StateSlots`](@ref) ordered as by
[`SlotContext`](@ref).

`ibp` maps an equation, or its index in the flattened equation list
(`_flatten_eqs(get_eqs(pdesys))`, the same order as the outputs), to an [`IBPPolicy`](@ref);
other keys are an error. Equations are matched by
their residual `lhs - rhs` up to sign, after flattening array variables through
`replaced_vars(v)`, so a key written before preprocessing still matches. A key that matches no
equation, or keys giving one equation different policies, are errors. The default is
[`ProductRule`](@ref), which warns once per call when it generates a correction term.
Under [`Weighted`](@ref)`(w)` the whole analysis is of the weighted equation `w⋅e`.
Equations are not paired with fields here ([`pair_equations`](@ref) is the discretizer's
helper). Spatial boundary conditions in `boundarymap` are classified by
[`boundary_conditions`](@ref) into essential, natural (matched against the slot-form fluxes)
and periodic conditions, keyed by `PDEBase.BoxFace` for a box discretization's map and by
`BoundaryPiece` for the `PieceBoundary` entries of a mesh discretization's map; initial
conditions are skipped.
"""
function analyse(pdesys, v::VariableMap, boundarymap; ibp = Dict())
    fields = declared_fields(pdesys, v)
    ctx = SlotContext(pdesys, v; fields)
    nfields, dim = length(fields), length(ctx.sivs)
    slots = StateSlots(nfields, dim, length(ctx.params))
    allivs = ctx.t === nothing ? ctx.sivs : vcat(ctx.sivs, ctx.t)
    strong, mass, remainder = Num[], Num[], Num[]
    flux, conormal = Vector{Num}[], Vector{Num}[]
    corrected = String[]
    eqs = _flatten_eqs(get_eqs(pdesys))
    for (index, (eq, policy)) in enumerate(zip(eqs, ibp_policies(ibp, eqs, v)))
        e = unwrap(eq.lhs - eq.rhs)
        policy isa Weighted && (e = unwrap(policy.w * e))
        check_whitelist(e, allivs, fields)
        m, spatial = split_time(e, fields, ctx.t)
        foreach(mf -> check_expandable(mf, ctx.sivs, fields), m)
        G, R, hascorrection = split_divergence(spatial, ctx.sivs, policy; depvars = fields)
        hascorrection && push!(corrected, "$index ($(nameof(typeof(policy))))")
        slot(ex) = slot_volume(ctx, ex, slots)
        push!(strong, slot_volume(ctx, e, slots; derivatives = true))
        push!(mass, sum(slot(m[f]) * slots.Ud[f] for f in 1:nfields))
        push!(flux, slot.(G))
        push!(remainder, slot(sum(R; init = Num(0))))
        Φ = sum(flux[end][d] * slots.N[d] for d in 1:dim)
        push!(conormal, [Symbolics.derivative(Φ, g) for g in slots.G])
        foreach(check_evaluable, vcat(mass[end], flux[end], remainder[end], conormal[end]))
    end
    essential, natural, periodic = boundary_conditions(boundarymap, ctx, slots, flux; mass, bcs = get_bcs(pdesys))
    isempty(corrected) || @warn "Integration by parts generated product-rule correction terms for equations $(join(corrected, ", ")), so their stiffness matrices are non-symmetric. Write them in divergence form, or pass `ibp = Dict(eq => Weighted(w))` or `Reject()`."
    return StrongFormAnalysis{_bc_key(essential)}(
        slots, nfields, dim, strong, mass, flux, remainder, conormal, essential, natural, periodic
    )
end

_bc_key(::AbstractVector{EssentialBC{K}}) where {K} = K
