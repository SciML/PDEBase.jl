
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
    eqs = get_eqs(pdesys)
    bcs = get_bcs(pdesys)
    dvs = get_dvs(pdesys)
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

    return PDESystem(eqs, bcs, get_domain(pdesys), get_ivs(pdesys), dvs, get_ps(pdesys),
                     defaults = get_defaults(pdesys), systems = get_systems(pdesys),
                     connector_type = get_connector_type(pdesys), metadata = get_metadata(pdesys),
                     analytic = getfield(pdesys, :analytic), analytic_func = getfield(pdesys, :analytic_func),
                     gui_metadata = get_gui_metadata(pdesys),
                     name = getfield(pdesys, :name)), replaced_vars

end

function split_complex_eq(eq, redvmaps, imdvmaps)
    eq = split_complex(eq)
    if eq isa Vector
        eq1 = eq[1]
        eq2 = eq[2]
        reeq1 = substitute(eq1.lhs, redvmaps) ~ substitute(eq1.rhs, redvmaps)
        imeq2 = substitute(eq2.lhs, imdvmaps) ~ substitute(eq2.rhs, imdvmaps)
        reeq2 = substitute(eq2.lhs, redvmaps) ~ substitute(eq2.rhs, redvmaps)
        imeq1 = substitute(eq1.lhs, imdvmaps) ~ substitute(eq1.rhs, imdvmaps)
        return [reeq1.lhs - imeq2.lhs ~ reeq1.rhs - imeq2.rhs , reeq2.lhs + imeq1.lhs ~ reeq2.rhs + imeq1.rhs]
    else
        eq1 = substitute(eq.lhs, redvmaps) ~ substitute(eq.rhs, redvmaps)
        eq2 = substitute(eq.lhs, imdvmaps) ~ substitute(eq.rhs, imdvmaps)
        return [eq1, eq2]
    end
end

struct ComplexEq
    reeq1
    imeq1
    reeq2
    imeq2
end

function split_complex_bc(eq, redvmaps, imdvmaps)
    eq = eq isa Pair ? eq : split_complex(eq) 
    if eq isa Vector 
        eq1 = eq[1]
        eq2 = eq[2]   
        reeq1 = substitute(eq1.lhs, redvmaps) ~ substitute(eq1.rhs, redvmaps)
        imeq2 = substitute(eq2.lhs, imdvmaps) ~ substitute(eq2.rhs, imdvmaps)
        reeq2 = substitute(eq2.lhs, redvmaps) ~ substitute(eq2.rhs, redvmaps)
        imeq1 = substitute(eq1.lhs, imdvmaps) ~ substitute(eq1.rhs, imdvmaps)
        return [reeq1.lhs - imeq2.lhs ~ reeq1.rhs - imeq2.rhs , reeq2.lhs + imeq1.lhs ~ reeq2.rhs + imeq1.rhs]
    elseif eq isa Pair
        rhs = split_complex(unwrap(eq.second))
        eq1 = substitute(eq.first, redvmaps) ~ rhs[1]
        eq2 = substitute(eq.first, imdvmaps) ~ rhs[2]
    else
        eq1 = substitute(eq1.lhs, redvmaps) ~ substitute(eq1.rhs, redvmaps)
        eq2 = substitute(eq2.lhs, imdvmaps) ~ substitute(eq2.rhs, imdvmaps)
    end
    return [eq1, eq2]
end

function handle_complex(pdesys)
    eqs = get_eqs(pdesys)
    if any(eq -> (eq isa Vector) || hascomplex(eq), eqs)
        dvmaps = map(get_dvs(pdesys)) do dv
            args = arguments(safe_unwrap(dv))
            dv = operation(safe_unwrap(dv))
            resym = Symbol("Re"*string(dv))
            imsym = Symbol("Im"*string(dv))
            redv = first(@variables $resym(..))
            imdv = first(@variables $imsym(..))
            redv = operation(unwrap(redv(args...)))
            imdv = operation(unwrap(imdv(args...)))
            (dv => redv, dv => imdv)
        end
        redvmaps = map(dvmaps) do dvmap
            dvmap[1]
        end
        imdvmaps = map(dvmaps) do dvmap
            dvmap[2]
        end
        dvmaps = Dict(map(dvmaps) do dvmap
            dvmap[1].first => (dvmap[1].second, dvmap[2].second)
        end)

        eqs = mapreduce(vcat, eqs) do eq
            split_complex_eq(eq, redvmaps, imdvmaps)
        end

        bcs = mapreduce(vcat, get_bcs(pdesys)) do eq
            split_complex_bc(eq, redvmaps, imdvmaps)
        end

        redvmaps = Dict(redvmaps)
        imdvmaps = Dict(imdvmaps)

        dvs = mapreduce(vcat, get_dvs(pdesys)) do dv
            dv = safe_unwrap(dv)
            redv = redvmaps[operation(dv)](arguments(dv)...)
            imdv = imdvmaps[operation(dv)](arguments(dv)...)
            [redv, imdv]
        end
        
        pdesys = PDESystem(eqs, bcs, get_domain(pdesys), get_ivs(pdesys), dvs, get_ps(pdesys), name = getfield(pdesys, :name))
        return pdesys, dvmaps
    else
        dvmaps = nothing
        return pdesys, dvmaps
    end
end

