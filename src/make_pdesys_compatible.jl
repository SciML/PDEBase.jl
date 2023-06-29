
function chain_flatten_array_variables(dvs)
	rs = []
	for dv in dvs
		dv = safe_unwrap(dv)
		if isequal(operation(dv), getindex)
			name = operation(arguments(dv)[1])
			args = arguments(arguments(dv)[1])
			idxs = arguments(dv)[2:end]
			fullname = Symbol(string(name) * "_" * string(idxs))
			newop = (@variables $fullname(..))[1]
			push!(rs, @rule getindex($(name)(~~a), idxs...) => newop(~a...))
		end
	end
	return isempty(rs) ? identity : Prewalk(Chain(rs))
end

function apply_lhs_rhs(f, eqs)
	map(eqs) do eq
		f(eq.lhs) ~ f(eq.rhs)
	end
end

function make_pdesys_compatible(pdesys::PDESystem)
    eqs = pdesys.eqs
    bcs = pdesys.bcs
    dvs = pdesys.dvs
    if any(u -> u isa Symbolics.Arr, dvs)
        dvs = reduce(vcat, collect.(dvs))
    end

    ch = chain_flatten_array_variables(dvs)
    safe_ch(x) = safe_unwrap(x) |> ch
    baddvs = filter(dvs) do u
        isequal(operation(safe_unwrap(u)), getindex)
    end
    replaced_vars = map(baddvs) do u
        safe_ch(u) => u
    end |> Dict
    eqs = apply_lhs_rhs(ch, eqs)
    bcs = apply_lhs_rhs(ch, bcs)
    dvs = map(safe_ch, dvs)

    return PDESystem(eqs, bcs, pdesys.domain, pdesys.ivs, dvs, pdesys.ps,
                     defaults = pdesys.defaults, systems = pdesys.systems,
                     connector_type = pdesys.connector_type, metadata = pdesys.metadata,
                     analytic = pdesys.analytic, analytic_func = pdesys.analytic_func,
                     gui_metadata = pdesys.gui_metadata,
                     name = pdesys.name), replaced_vars

end
