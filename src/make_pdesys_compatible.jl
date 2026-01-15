# Helper function to replace operations in symbolic terms
# In SymbolicUtils v4, substitute with operation mappings doesn't work,
# so we use manual tree traversal instead
function _replace_ops(term, op_map)
    # Handle Num wrapper
    if term isa Num
        return Num(_replace_ops(unwrap(term), op_map))
    end
    # Handle symbolic terms
    term_unwrapped = safe_unwrap(term)
    if iscall(term_unwrapped)
        old_op = operation(term_unwrapped)
        args = arguments(term_unwrapped)
        new_args = [_replace_ops(a, op_map) for a in args]
        if haskey(op_map, old_op)
            new_op = op_map[old_op]
            return new_op(new_args...)
        else
            if all(isequal.(args, new_args))
                return term
            else
                return maketerm(typeof(term_unwrapped), old_op, new_args, metadata(term_unwrapped))
            end
        end
    end
    return term
end

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
    return map(eqs) do eq
        if eq isa Pair
            # Handle initial conditions specified as Pairs (u(0,x) => value)
            f(eq.first) => f(eq.second)
        else
            f(eq.lhs) ~ f(eq.rhs)
        end
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

    return PDESystem(
            eqs, bcs, get_domain(pdesys), get_ivs(pdesys), dvs, get_ps(pdesys),
            initial_conditions = pdesys.initial_conditions, systems = get_systems(pdesys),
            connector_type = get_connector_type(pdesys), metadata = get_metadata(pdesys),
            analytic = getfield(pdesys, :analytic), analytic_func = getfield(pdesys, :analytic_func),
            gui_metadata = get_gui_metadata(pdesys),
            name = getfield(pdesys, :name)
        ),
        replaced_vars
end

function split_complex_eq(eq, redvmaps, imdvmaps)
    eq = split_complex(eq)
    if eq isa Vector
        eq1 = eq[1]
        eq2 = eq[2]
        reeq1 = _replace_ops(eq1.lhs, redvmaps) ~ _replace_ops(eq1.rhs, redvmaps)
        imeq2 = _replace_ops(eq2.lhs, imdvmaps) ~ _replace_ops(eq2.rhs, imdvmaps)
        reeq2 = _replace_ops(eq2.lhs, redvmaps) ~ _replace_ops(eq2.rhs, redvmaps)
        imeq1 = _replace_ops(eq1.lhs, imdvmaps) ~ _replace_ops(eq1.rhs, imdvmaps)
        return [
            reeq1.lhs - imeq2.lhs ~ reeq1.rhs - imeq2.rhs,
            reeq2.lhs + imeq1.lhs ~ reeq2.rhs + imeq1.rhs,
        ]
    else
        eq1 = _replace_ops(eq.lhs, redvmaps) ~ _replace_ops(eq.rhs, redvmaps)
        eq2 = _replace_ops(eq.lhs, imdvmaps) ~ _replace_ops(eq.rhs, imdvmaps)
        return [eq1, eq2]
    end
end

struct ComplexEq
    reeq1::Any
    imeq1::Any
    reeq2::Any
    imeq2::Any
end

function split_complex_bc(eq, redvmaps, imdvmaps)
    # Handle pre-split complex BCs from Symbolics v7
    # These are already split into real/imag parts, just need variable renaming
    if eq isa PreSplitComplexBC
        # eq.real_eq becomes the Reψ equation, eq.imag_eq becomes the Imψ equation
        real_renamed = _replace_ops(eq.real_eq.lhs, redvmaps) ~ _replace_ops(eq.real_eq.rhs, redvmaps)
        imag_renamed = _replace_ops(eq.imag_eq.lhs, imdvmaps) ~ _replace_ops(eq.imag_eq.rhs, imdvmaps)
        return [real_renamed, imag_renamed]
    end

    # For Pair type (initial conditions), handle specially
    if eq isa Pair
        rhs = split_complex(unwrap(eq.second))
        eq1 = _replace_ops(eq.first, redvmaps) ~ rhs[1]
        eq2 = _replace_ops(eq.first, imdvmaps) ~ rhs[2]
        return [eq1, eq2]
    end

    # For Equation type, check if it actually has complex values
    if !hascomplex(eq)
        # Real BC: just duplicate with variable renaming
        # e.g., ψ(t, 0) ~ 0  becomes  Reψ(t, 0) ~ 0  and  Imψ(t, 0) ~ 0
        eq1 = _replace_ops(eq.lhs, redvmaps) ~ _replace_ops(eq.rhs, redvmaps)
        eq2 = _replace_ops(eq.lhs, imdvmaps) ~ _replace_ops(eq.rhs, imdvmaps)
        return [eq1, eq2]
    end

    # Complex BC: split into real and imaginary parts
    eq_split = split_complex(eq)
    if eq_split isa Vector
        eq1 = eq_split[1]
        eq2 = eq_split[2]
        reeq1 = _replace_ops(eq1.lhs, redvmaps) ~ _replace_ops(eq1.rhs, redvmaps)
        imeq2 = _replace_ops(eq2.lhs, imdvmaps) ~ _replace_ops(eq2.rhs, imdvmaps)
        reeq2 = _replace_ops(eq2.lhs, redvmaps) ~ _replace_ops(eq2.rhs, redvmaps)
        imeq1 = _replace_ops(eq1.lhs, imdvmaps) ~ _replace_ops(eq1.rhs, imdvmaps)
        return [
            reeq1.lhs - imeq2.lhs ~ reeq1.rhs - imeq2.rhs,
            reeq2.lhs + imeq1.lhs ~ reeq2.rhs + imeq1.rhs,
        ]
    else
        eq1 = _replace_ops(eq.lhs, redvmaps) ~ _replace_ops(eq.rhs, redvmaps)
        eq2 = _replace_ops(eq.lhs, imdvmaps) ~ _replace_ops(eq.rhs, imdvmaps)
        return [eq1, eq2]
    end
end

function handle_complex(pdesys)
    eqs = get_eqs(pdesys)
    bcs = get_bcs(pdesys)
    # In MTK v11, complex equations may already be nested Vector{Equation}
    # Flatten first before processing
    eqs_flat = _flatten_eqs(eqs)
    bcs_flat = _flatten_bcs(bcs)

    eqs_have_complex = any(eq -> hascomplex(eq), eqs_flat)
    bcs_have_complex = any(bc -> hascomplex(bc), bcs_flat)

    # Check both equations and BCs for complex values
    if eqs_have_complex || bcs_have_complex
        dvmaps = map(get_dvs(pdesys)) do dv
            args = arguments(safe_unwrap(dv))
            dv = operation(safe_unwrap(dv))
            resym = Symbol("Re" * string(dv))
            imsym = Symbol("Im" * string(dv))
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
        dvmaps = Dict(
            map(dvmaps) do dvmap
                dvmap[1].first => (dvmap[1].second, dvmap[2].second)
            end
        )

        # Convert to Dict before calling split functions (required for SymbolicUtils v4)
        redvmaps_dict = Dict(redvmaps)
        imdvmaps_dict = Dict(imdvmaps)

        if eqs_have_complex
            # Equations have complex values - split them into real/imaginary parts
            eqs = mapreduce(vcat, eqs_flat) do eq
                split_complex_eq(eq, redvmaps_dict, imdvmaps_dict)
            end
        else
            # Equations are already real (MTK v11 may have pre-split them)
            # Just rename variables without re-splitting
            # In MTK v11, we expect equations to come in pairs (real, imag)
            # Map them directly to Reψ and Imψ equations
            n_eqs = length(eqs_flat)
            if n_eqs % 2 == 0
                # Assume first half are "real part" equations, second half are "imag part"
                # Just replace ψ with Reψ in first half and Imψ in second half
                half = n_eqs ÷ 2
                eqs = vcat(
                    [_replace_ops(eq.lhs, redvmaps_dict) ~ _replace_ops(eq.rhs, redvmaps_dict) for eq in eqs_flat[1:half]],
                    [_replace_ops(eq.lhs, imdvmaps_dict) ~ _replace_ops(eq.rhs, imdvmaps_dict) for eq in eqs_flat[half+1:end]]
                )
            else
                # Odd number of equations - just rename all with real maps
                eqs = [_replace_ops(eq.lhs, redvmaps_dict) ~ _replace_ops(eq.rhs, redvmaps_dict) for eq in eqs_flat]
            end
        end

        bcs = mapreduce(vcat, bcs_flat) do eq
            split_complex_bc(eq, redvmaps_dict, imdvmaps_dict)
        end

        dvs = mapreduce(vcat, get_dvs(pdesys)) do dv
            dv = safe_unwrap(dv)
            redv = redvmaps_dict[operation(dv)](arguments(dv)...)
            imdv = imdvmaps_dict[operation(dv)](arguments(dv)...)
            [redv, imdv]
        end

        pdesys = PDESystem(
            eqs, bcs, get_domain(pdesys), get_ivs(pdesys), dvs,
            get_ps(pdesys), name = getfield(pdesys, :name),
            initial_conditions = pdesys.initial_conditions
        )
        return pdesys, dvmaps
    else
        dvmaps = nothing
        # Even if no complex equations need splitting, we still need to flatten
        # nested equations that may have been created by MTK v11's equation processing
        # (bcs_flat was already computed above)
        pdesys = PDESystem(
            eqs_flat, bcs_flat, get_domain(pdesys), get_ivs(pdesys), get_dvs(pdesys),
            get_ps(pdesys), name = getfield(pdesys, :name),
            initial_conditions = pdesys.initial_conditions
        )
        return pdesys, dvmaps
    end
end
