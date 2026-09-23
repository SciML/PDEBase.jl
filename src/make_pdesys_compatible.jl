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

_dependent_variable_term(dv) = safe_unwrap(dv)
_dependent_variable_term(dv::Complex{Num}) = only(arguments(unwrap(real(dv))))

function _rename_complex_typed_components(term, redvmaps, imdvmaps)
    term isa Num && return _rename_complex_typed_components(unwrap(term), redvmaps, imdvmaps)
    if term isa Complex{Num}
        return complex(
            _rename_complex_typed_components(real(term), redvmaps, imdvmaps),
            _rename_complex_typed_components(imag(term), redvmaps, imdvmaps)
        )
    end
    iscall(term) || return term

    op = operation(term)
    args = arguments(term)
    if (op === real || op === imag) && length(args) == 1 && iscall(only(args))
        dv = only(args)
        dvmap = op === real ? redvmaps : imdvmaps
        if haskey(dvmap, operation(dv))
            return dvmap[operation(dv)](arguments(dv)...)
        end
    end
    new_args = [_rename_complex_typed_components(arg, redvmaps, imdvmaps) for arg in args]
    return all(isequal.(args, new_args)) ? term :
        maketerm(typeof(term), op, new_args, metadata(term))
end

function chain_flatten_array_variables(dvs)
    rs = []
    for dv in dvs
        dv = _dependent_variable_term(dv)
        if isequal(operation(dv), getindex)
            name = operation(arguments(dv)[1])
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
        if eq isa AbstractVector
            apply_lhs_rhs(f, eq)
        elseif eq isa Pair
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
        isequal(operation(_dependent_variable_term(u)), getindex)
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

function _split_complex_components(term, preserve_differentials = false)
    term = safe_unwrap(term)
    if preserve_differentials && iscall(term) && operation(term) isa Differential
        op = operation(term)
        arg = first(arguments(term))
        constant_arg = safe_unwrap(unwrap_const(arg))
        if constant_arg isa Number
            return 0, 0
        end
        try
            Symbolics.value(constant_arg) isa Number && return 0, 0
        catch
        end
        re, im = _split_complex_components(arg, preserve_differentials)
        return op(re), op(im)
    end
    term = unwrap_const(term)
    term isa Number && return real(term), imag(term)
    iscall(term) || return term, 0

    op = operation(term)
    args = arguments(term)
    if op === (+)
        real_part, imag_part = _split_complex_components(first(args), preserve_differentials)
        for arg in Iterators.drop(args, 1)
            argre, argim = _split_complex_components(arg, preserve_differentials)
            real_part += argre
            imag_part += argim
        end
        return real_part, imag_part
    elseif op === (-)
        real_part, imag_part = _split_complex_components(first(args), preserve_differentials)
        if length(args) == 1
            return -real_part, -imag_part
        end
        for arg in Iterators.drop(args, 1)
            argre, argim = _split_complex_components(arg, preserve_differentials)
            real_part -= argre
            imag_part -= argim
        end
        return real_part, imag_part
    elseif op === (*)
        real_part, imag_part = _split_complex_components(first(args), preserve_differentials)
        for arg in Iterators.drop(args, 1)
            argre, argim = _split_complex_components(arg, preserve_differentials)
            real_part, imag_part = real_part * argre - imag_part * argim,
                real_part * argim + imag_part * argre
        end
        return real_part, imag_part
    elseif op === (/)
        are, aim = _split_complex_components(args[1], preserve_differentials)
        b_realpart, bim = _split_complex_components(args[2], preserve_differentials)
        denom = b_realpart^2 + bim^2
        return (are * b_realpart + aim * bim) / denom, (aim * b_realpart - are * bim) / denom
    end
    return real(term), imag(term)
end

function _is_false_constant(term)
    term = unwrap_const(safe_unwrap(term))
    term === false && return true
    return false
end

function _dependent_variable_instances!(found, term, dependent_operations)
    term = safe_unwrap(term)
    if iscall(term)
        any(op -> isequal(operation(term), op), dependent_operations) && push!(found, term)
        foreach(arg -> _dependent_variable_instances!(found, arg, dependent_operations), arguments(term))
    end
    return found
end

function _same_dependent_variable_instances(eq1, eq2, dependent_operations)
    instances(eq) = unique(vcat(
        _dependent_variable_instances!(Any[], eq.lhs, dependent_operations),
        _dependent_variable_instances!(Any[], eq.rhs, dependent_operations)
    ))
    first_instances, second_instances = instances(eq1), instances(eq2)
    return !isempty(first_instances) && length(first_instances) == length(second_instances) &&
        all(x -> any(y -> isequal(x, y), second_instances), first_instances)
end

_is_zero_constant(term) = isequal(unwrap_const(safe_unwrap(term)), 0)

function _is_presplit_bc_candidate(eq1, eq2, dependent_operations)
    first_instances = _dependent_variable_instances!(Any[], eq1.lhs, dependent_operations)
    append!(first_instances, _dependent_variable_instances!(Any[], eq1.rhs, dependent_operations))
    second_instances = _dependent_variable_instances!(Any[], eq2.lhs, dependent_operations)
    append!(second_instances, _dependent_variable_instances!(Any[], eq2.rhs, dependent_operations))
    _same_dependent_variable_instances(eq1, eq2, dependent_operations) && return true
    return (isempty(first_instances) && _is_zero_constant(eq1.lhs) && !isempty(second_instances)) ||
        (isempty(second_instances) && _is_zero_constant(eq2.lhs) && !isempty(first_instances))
end

function _ambiguous_presplit_bc(bcs, dependent_operations)
    for bc in bcs
        if bc isa AbstractVector
            if length(bc) == 2 && all(eq -> eq isa Equation, bc) &&
                    _is_presplit_bc_candidate(bc[1], bc[2], dependent_operations)
                return true
            end
            _ambiguous_presplit_bc(bc, dependent_operations) && return true
        end
    end
    return false
end

function split_complex_eq(eq, redvmaps, imdvmaps; expand_derivatives = true)
    lhs, rhs = if eq isa AbstractVector
        real_eq, imag_eq = eq
        real_eq.lhs + im * imag_eq.lhs, real_eq.rhs + im * imag_eq.rhs
    else
        eq.lhs, eq.rhs
    end
    complexmap = Dict(
        op => ((args...) -> redop(args...) + im * imdvmaps[op](args...))
            for (op, redop) in redvmaps
    )
    if expand_derivatives
        lhs = Symbolics.expand(Symbolics.expand_derivatives(_replace_ops(lhs, complexmap)))
        rhs = Symbolics.expand(Symbolics.expand_derivatives(_replace_ops(rhs, complexmap)))
    else
        lhs = Symbolics.expand(_replace_ops(lhs, complexmap))
        rhs = Symbolics.expand(_replace_ops(rhs, complexmap))
    end
    lhsre, lhsim = _split_complex_components(lhs, !expand_derivatives)
    rhsre, rhsim = _split_complex_components(rhs, !expand_derivatives)
    return [lhsre ~ rhsre, lhsim ~ rhsim]
end

function split_complex_bc(eq, redvmaps, imdvmaps)
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

    return split_complex_eq(eq, redvmaps, imdvmaps)
end

function handle_complex(pdesys)
    eqs = get_eqs(pdesys)
    bcs = get_bcs(pdesys)
    typed_dvs = all(dv -> dv isa Complex{Num}, get_dvs(pdesys)) && !isempty(get_dvs(pdesys))
    # In MTK v11, complex equations may already be nested Vector{Equation}
    # Flatten first before processing
    eqs_flat = _flatten_eqs(eqs)
    dependent_operations = map(dv -> operation(_dependent_variable_term(dv)), get_dvs(pdesys))
    eqs_have_complex = any(eq -> hascomplex(eq), eqs_flat) || any(eq -> eq isa AbstractVector, eqs)
    if !typed_dvs && _ambiguous_presplit_bc(bcs, dependent_operations)
        throw(ArgumentError(
            "Symbolics has pre-split a complex boundary condition before PDEBase can verify its meaning. " *
                "Declare the dependent variable as `::Complex` to preserve complex boundary expressions."
        ))
    end
    bcs_flat = _flatten_bcs(bcs)

    if any(bc -> bc isa Equation && (_is_false_constant(bc.lhs) || _is_false_constant(bc.rhs)), bcs_flat)
        throw(
            ArgumentError(
                "A boundary condition reduced to `false` before PDEBase could inspect its complex value. " *
                    "Symbolics may truncate complex constants in `~`; declare the dependent variable as `::Complex`."
            )
        )
    end

    bcs_have_complex = any(bc -> hascomplex(bc), bcs_flat)

    # Check both equations and BCs for complex values
    if eqs_have_complex || bcs_have_complex || typed_dvs
        dvmaps = map(get_dvs(pdesys)) do dv
            dv = _dependent_variable_term(dv)
            args = arguments(dv)
            dv = operation(dv)
            resym = Symbol("Re" * string(dv))
            imsym = Symbol("Im" * string(dv))
            redv = first(@variables $resym(..)::Real)
            imdv = first(@variables $imsym(..)::Real)
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

        if typed_dvs
            rename(term) = _rename_complex_typed_components(term, redvmaps_dict, imdvmaps_dict)
            eqs = [rename(eq.lhs) ~ rename(eq.rhs) for eq in eqs_flat]
            bcs = mapreduce(vcat, bcs_flat) do bc
                if bc isa Pair
                    dv = _dependent_variable_term(bc.first)
                    op = operation(dv)
                    args = arguments(dv)
                    [
                        redvmaps_dict[op](args...) ~ rename(real(bc.second)),
                        imdvmaps_dict[op](args...) ~ rename(imag(bc.second)),
                    ]
                else
                    [rename(bc.lhs) ~ rename(bc.rhs)]
                end
            end
        elseif eqs_have_complex || bcs_have_complex
            # Equations have complex values - split them into real/imaginary parts
            eqs = mapreduce(vcat, eqs) do eq
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
                    [_replace_ops(eq.lhs, imdvmaps_dict) ~ _replace_ops(eq.rhs, imdvmaps_dict) for eq in eqs_flat[(half + 1):end]]
                )
            else
                # Odd number of equations - just rename all with real maps
                eqs = [_replace_ops(eq.lhs, redvmaps_dict) ~ _replace_ops(eq.rhs, redvmaps_dict) for eq in eqs_flat]
            end
        end

        if !typed_dvs
            bcs = mapreduce(vcat, bcs_flat) do eq
                split_complex_bc(eq, redvmaps_dict, imdvmaps_dict)
            end
        end

        dvs = mapreduce(vcat, get_dvs(pdesys)) do dv
            dv = _dependent_variable_term(dv)
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
