# Equation-to-field pairing (spec §5.4): equation `e` is tested with the test function of
# `pair_equations(...)[e]`.

"""
    pair_equations(mass, slots; pairing = Dict{Int, Int}()) -> Vector{Int}

The field each equation is paired with. `mass[e]` is equation `e`'s time-derivative part
`Σ_f m_ef Ud[f]` in slots. An explicit `pairing[e] = f` wins; otherwise the equation is paired with
the one field whose time derivative it contains. An algebraic row (no time derivative) or a row
holding several fields' time derivatives must be named in `pairing`. Errors if a field is paired
twice or `pairing` names an equation or field that does not exist.
"""
function pair_equations(mass, slots; pairing = Dict{Int, Int}())
    for (e, f) in pairing
        e in eachindex(mass) || error("`pairing` names equation $e, but there are $(length(mass)) equations")
        f in eachindex(slots.U) ||
            error("`pairing` pairs equation $e with field $f, but there are $(length(slots.U)) fields")
    end
    fields = map(eachindex(mass)) do e
        haskey(pairing, e) && return pairing[e]
        owned = [f for f in eachindex(slots.Ud) if _depends(mass[e], [slots.Ud[f]])]
        isempty(owned) &&
            error("equation $e has no time derivative; name the field it is tested with in `pairing`")
        length(owned) == 1 ||
            error("equation $e contains the time derivatives of fields $owned; name the field it is tested with in `pairing`")
        return only(owned)
    end
    for f in unique(fields)
        eqs = findall(==(f), fields)
        length(eqs) == 1 ||
            error("field $f is paired with equations $(eqs[1]) and $(eqs[2]); name the pairing explicitly")
    end
    return fields
end
