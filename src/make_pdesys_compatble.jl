function make_pdesys_compatible(pdesys::PDESystem)
    dvs = pdesys.dvs
    eqs = pdesys.eqs
    bcs = pdesys.bcs

    replaced_vars = Dict()
    ch = chain_flatten_array_variables(dvs)
    safe_ch(x) = safe_unwrap(x) |> ch
    baddvs = filter(dvs) do u
        isequal(operation(safe_unwrap(u)), getindex)
    end
    _replacedvars = map(baddvs) do u
        u => Num(safe_ch(u))
    end |> Dict
    merge!(replaced_vars, _replacedvars)

    eqs = apply_lhs_rhs(ch, eqs)
    bcs = apply_lhs_rhs(ch, bcs)
    dvs = map(safe_ch, dvs)

    return PDESystem(eqs, bcs, pdesys.ivs, dvs, pdesys.domain, pdesys.ps)
end
