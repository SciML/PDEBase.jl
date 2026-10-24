same_residual(a, b) =
    isequal(Symbolics.expand(Num(a - b)), Num(0)) || isequal(Symbolics.expand(Num(a + b)), Num(0))

# Equation keys are matched by residual up to sign, since preprocessing such as
# `cardinalize_eqs!` rewrites `lhs ~ rhs` as `lhs - rhs ~ 0`. A key matching nothing is an error.
function ibp_policies(ibp, eqs)
    policies = IBPPolicy[ProductRule() for _ in eqs]
    residuals = [eq.lhs - eq.rhs for eq in eqs]
    for (key, policy) in ibp
        matched = key isa Integer ? filter(==(key), eachindex(eqs)) :
            findall(r -> same_residual(r, key.lhs - key.rhs), residuals)
        isempty(matched) && error("ibp key $key matches no equation of the system")
        policies[matched] .= Ref(policy)
    end
    return policies
end

function check_no_boundary_conditions(boundarymap, ctx::SlotContext)
    for (_, bykey) in boundarymap, (x, bs) in bykey
        ctx.t !== nothing && isequal(x, ctx.t) && continue
        isempty(bs) ||
            error("boundary-condition analysis is not implemented yet; analyse needs a boundary map without spatial conditions")
    end
    return nothing
end

"""
    analyse(pdesys, v::VariableMap, boundarymap; ibp = Dict()) -> StrongFormAnalysis

Strong-form analysis of `pdesys`, which must already be in the form `symbolic_discretize`
hands to the discretizer hooks (after `make_pdesys_compatible`), with `v` its `VariableMap`
and `boundarymap` from `parse_bcs`. Each equation `e = lhs - rhs` passes the operator whitelist
([`check_whitelist`](@ref)), is split in time ([`split_time`](@ref)) and into divergence form
([`split_divergence`](@ref)), and is written in [`StateSlots`](@ref) ordered as by
[`SlotContext`](@ref).

`ibp` maps an equation, or its index, to an [`IBPPolicy`](@ref). Equations are matched by
their residual `lhs - rhs` up to sign, so a key written before preprocessing still matches, and
a key that matches no equation is an error. The default is
[`ProductRule`](@ref), which warns once per call when it generates a correction term.
Under [`Weighted`](@ref)`(w)` the whole analysis is of the weighted equation `w⋅e`.
Equations are not paired with fields here. Boundary conditions are not analysed yet: a
boundary map holding any spatial condition is an error.
"""
function analyse(pdesys, v::VariableMap, boundarymap; ibp = Dict())
    ctx = SlotContext(v)
    fields = depvars(v)
    nfields, dim = length(fields), length(ctx.sivs)
    slots = StateSlots(nfields, dim, length(ctx.params))
    allivs = ctx.t === nothing ? ctx.sivs : vcat(ctx.sivs, ctx.t)
    strong, mass, remainder = Num[], Num[], Num[]
    flux, conormal = Vector{Num}[], Vector{Num}[]
    corrected = Int[]
    eqs = _flatten_eqs(get_eqs(pdesys))
    for (index, (eq, policy)) in enumerate(zip(eqs, ibp_policies(ibp, eqs)))
        e = unwrap(eq.lhs - eq.rhs)
        policy isa Weighted && (e = unwrap(policy.w * e))
        check_whitelist(e, allivs, fields)
        m, spatial = split_time(e, fields, ctx.t)
        G, R, hascorrection = split_divergence(spatial, ctx.sivs, policy)
        hascorrection && push!(corrected, index)
        slot(ex) = slot_volume(ctx, ex, slots)
        push!(strong, slot(e))
        push!(mass, sum(slot(m[f]) * slots.Ud[f] for f in 1:nfields))
        push!(flux, slot.(G))
        push!(remainder, slot(sum(R; init = Num(0))))
        Φ = sum(flux[end][d] * slots.N[d] for d in 1:dim)
        push!(conormal, [Symbolics.derivative(Φ, g) for g in slots.G])
    end
    check_no_boundary_conditions(boundarymap, ctx)
    isempty(corrected) || @warn "The ProductRule integration-by-parts policy generated correction terms for equations $corrected, so their stiffness matrices are non-symmetric. Write them in divergence form, or pass `ibp = Dict(eq => Weighted(w))` or `Reject()`."
    return StrongFormAnalysis{Any}(
        slots, nfields, dim, strong, mass, flux, remainder, conormal,
        EssentialBC{Any}[], NaturalBC{Any}[], PeriodicBC{Any}[]
    )
end
