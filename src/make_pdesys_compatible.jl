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
        if old_op === complex && length(args) == 2
            # `complex(re, im)` only accepts Real-typed arguments, so a
            # substitution that makes an argument complex cannot rebuild the
            # term; `complex(a, b)` is `a + b*im` arithmetic.
            return _replace_ops(args[1], op_map) + im * _replace_ops(args[2], op_map)
        end
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
_dependent_variable_term(dv::Symbolics.Arr) = unwrap(dv)
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

"""
    make_pdesys_compatible(pdesys::PDESystem) -> (pdesys, replaced_vars)

Flatten the array-valued dependent variables of `pdesys` into scalar variables, rewriting its
equations, boundary conditions and unknowns. Returns the new system and a `Dict` mapping each
flattened scalar variable to the original indexed variable, the `replaced_vars` that
`VariableMap(pdesys, disc; replaced_vars)` records. Systems without array-valued dependent
variables come back with an empty map. `symbolic_discretize` calls it before building the
`VariableMap`; discretizers that run their own pipeline, such as one calling
`StrongForm.analyse`, must too.
"""
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

_is_false_constant(term) = unwrap_const(safe_unwrap(term)) === false

function _dependent_variable_instances!(found, term, dependent_operations)
    term = safe_unwrap(term)
    if iscall(term)
        any(op -> isequal(operation(term), op), dependent_operations) && push!(found, term)
        foreach(arg -> _dependent_variable_instances!(found, arg, dependent_operations), arguments(term))
    end
    return found
end

function _same_dependent_variable_instances(eq1, eq2, dependent_operations)
    instances(eq) = unique(
        vcat(
            _dependent_variable_instances!(Any[], eq.lhs, dependent_operations),
            _dependent_variable_instances!(Any[], eq.rhs, dependent_operations)
        )
    )
    first_instances, second_instances = instances(eq1), instances(eq2)
    return !isempty(first_instances) && length(first_instances) == length(second_instances) &&
        all(x -> any(y -> isequal(x, y), second_instances), first_instances)
end

_is_zero_constant(term) = isequal(unwrap_const(safe_unwrap(term)), 0)

function _dv_application(term, dependent_operations)
    term = safe_unwrap(term)
    iscall(term) || return false
    op = operation(term)
    op isa Differential && return _dv_application(first(arguments(term)), dependent_operations)
    return any(o -> isequal(op, o), dependent_operations)
end

_dv_free(term, dependent_operations) =
    isempty(_dependent_variable_instances!(Any[], term, dependent_operations))

const _MARKER_HINT = isdefined(Symbolics, :split_complex_equation) ?
    " or build the condition with `Symbolics.split_complex_equation`" : ""

_presplit_bc_error() = ArgumentError(
    "Symbolics has pre-split a complex boundary condition before PDEBase can verify its meaning, " *
        "or a `0 ~` condition was nested inside a boundary group. Declare the dependent " *
        "variable as `::Complex`" * _MARKER_HINT * " to preserve complex boundary expressions; " *
        "if these are real boundary conditions, pass them un-nested."
)

# Symbolics rewrites a complex `lhs ~ rhs` over real-typed variables to
# `[real(lhs) ~ real(rhs), imag(lhs) ~ imag(rhs)]`; for a real-typed `lhs`
# the imaginary half carries a literal `0` left-hand side. A nested pair
# `[L ~ a, 0 ~ b]` where `L` is a single dependent-variable application,
# possibly under a `Differential`, and neither `a` nor `b` mentions a
# dependent variable can only have come from splitting `L(ψ) ~ a + i b` — as
# a real grouping `0 ~ b` constrains no unknown — so it is rebuilt as that
# equation. A `[0 ~ a, 0 ~ b]` pair is a real grouping when both members
# mention a dependent variable, and the split of a residual-form complex
# condition (`0 ~ ψ(t, 0) - exp(im*t)` pre-splits to `[0 ~ ψ - cos(t),
# 0 ~ -sin(t)]`) when exactly one member is free of dependent variables —
# as a real grouping such a member constrains no unknown — so it is rebuilt
# as `0 ~ a + i b` and split downstream. Returns `nothing` when neither
# reading applies.
function _rebuild_presplit_pair(eq1, eq2, dependent_operations)
    if _is_zero_constant(eq1.lhs)
        free1 = _dv_free(eq1.rhs, dependent_operations)
        free2 = _dv_free(eq2.rhs, dependent_operations)
        free1 && free2 && throw(_presplit_bc_error())
        return free1 ⊻ free2 ? Equation(0, eq1.rhs + im * eq2.rhs) : nothing
    elseif _dv_application(eq1.lhs, dependent_operations) &&
            _dv_free(eq1.rhs, dependent_operations) &&
            _dv_free(eq2.rhs, dependent_operations)
        return Equation(eq1.lhs, eq1.rhs + im * eq2.rhs)
    end
    return nothing
end

_flip_equation(eq) = Equation(eq.rhs, eq.lhs)

_constrains_unknown(eq, dependent_operations) =
    !(_dv_free(eq.lhs, dependent_operations) && _dv_free(eq.rhs, dependent_operations))

# A condition written with the field on the right, `exp(im*t) ~ ψ(t, 0)` or
# `ψ(t, 0) - exp(im*t) ~ 0`, pre-splits to the mirror image `[a ~ L, b ~ 0]`
# or `[a ~ 0, b ~ 0]` and is read the same way with its sides swapped; an
# unmatched mirror is an ordinary real grouping such as `[u ~ 0, Dx(u) ~ 0]`.
# A `[0 ~ a, L ~ b]` pair whose members both constrain an unknown is left for
# `_ambiguous_presplit_bc`. Any other nested pair containing a `0 ~` member is
# rejected: dependent-variable data behind the `0 ~` (the pre-split of e.g.
# `L ~ i*ψ(t, 0)`) couples the fields in a way this shape does not record.
function _resolve_presplit_bc_pairs(bcs, dependent_operations)
    return map(bcs) do bc
        bc isa AbstractVector || return bc
        length(bc) == 2 && all(eq -> eq isa Equation, bc) ||
            return _resolve_presplit_bc_pairs(bc, dependent_operations)
        eq1, eq2 = bc
        if _is_zero_constant(eq2.lhs)
            rebuilt = _rebuild_presplit_pair(eq1, eq2, dependent_operations)
            rebuilt === nothing && !_is_zero_constant(eq1.lhs) && throw(_presplit_bc_error())
            return something(rebuilt, bc)
        elseif _is_zero_constant(eq1.lhs)
            _constrains_unknown(eq1, dependent_operations) &&
                _constrains_unknown(eq2, dependent_operations) || throw(_presplit_bc_error())
            return bc
        elseif _is_zero_constant(eq2.rhs)
            rebuilt = _rebuild_presplit_pair(
                _flip_equation(eq1), _flip_equation(eq2), dependent_operations
            )
            return rebuilt === nothing ? bc : _flip_equation(rebuilt)
        end
        return bc
    end
end

# A same-point pair, or a `[0 ~ a, L ~ b]` pair (the pre-split of
# `i*L ~ a + i b`), can equally be user grouping, so it is only flagged when
# the system already shows complex values.
function _ambiguous_presplit_bc(bcs, dependent_operations, complex_evidence)
    complex_evidence || return false
    for bc in bcs
        if bc isa AbstractVector
            if length(bc) == 2 && all(eq -> eq isa Equation, bc) && (
                    _same_dependent_variable_instances(bc[1], bc[2], dependent_operations) ||
                        _is_zero_constant(bc[1].lhs) && !_is_zero_constant(bc[2].lhs)
                )
                return true
            end
            _ambiguous_presplit_bc(bc, dependent_operations, complex_evidence) && return true
        end
    end
    return false
end

# A `Symbolics.SplitComplexEquation` records that its two equations came from
# splitting one complex `~` and keeps the original equation, so restoring it
# routes the condition through the coupled split instead of the ambiguous-pair
# checks that a hand-written or plain `~` pair still receives.
@static if isdefined(Symbolics, :SplitComplexEquation)
    _restore_marked_complex_eq(eq::Symbolics.SplitComplexEquation) = eq.original
end
_restore_marked_complex_eq(bcs::AbstractVector) = map(_restore_marked_complex_eq, bcs)
_restore_marked_complex_eq(bc) = bc

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
    # For Pair type (initial conditions), treat the pair as an equation so the
    # right-hand side is split through the complex map like an equation's rhs.
    eq isa Pair && (eq = eq.first ~ eq.second)

    # Boundary conditions are evaluated at points, where expand_derivatives
    # reduces `Differential(x)(f(t, 1))` to 0 because the boundary argument
    # does not contain x; preserve differentials instead of expanding them.
    return split_complex_eq(eq, redvmaps, imdvmaps; expand_derivatives = false)
end

function handle_complex(pdesys)
    eqs = get_eqs(pdesys)
    bcs = get_bcs(pdesys)
    dvs = get_dvs(pdesys)
    # A real field promoted to Complex{Num} has no symbolic real(...) wrapper.
    complex_typed = map(dvs) do dv
        dv isa Complex{Num} && iscall(unwrap(real(dv))) &&
            operation(unwrap(real(dv))) === real
    end
    if any(complex_typed) && !all(complex_typed)
        throw(ArgumentError("Complex-typed and real dependent variables cannot be mixed in handle_complex"))
    end
    typed_dvs = !isempty(dvs) && all(complex_typed)
    dependent_operations = map(dv -> operation(_dependent_variable_term(dv)), get_dvs(pdesys))
    # A `SplitComplexEquation` marker in the boundary conditions is restored to
    # the complex equation it was split from, so it is split below with its
    # real and imaginary parts coupled. Typed fields do not need the marker:
    # their real(ψ)/imag(ψ) parts survive `~` and rename directly.
    if !typed_dvs
        bcs = _restore_marked_complex_eq(bcs)
        bcs = _resolve_presplit_bc_pairs(bcs, dependent_operations)
    end
    # In MTK v11, complex equations may already be nested Vector{Equation}
    # Flatten first before processing
    eqs_flat = _flatten_eqs(eqs)
    eqs_have_complex = any(eq -> hascomplex(eq), eqs_flat) || any(eq -> eq isa AbstractVector, eqs)
    bcs_flat = _flatten_bcs(bcs)
    bcs_have_complex = any(bc -> hascomplex(bc), bcs_flat)
    if !typed_dvs && _ambiguous_presplit_bc(bcs, dependent_operations, eqs_have_complex || bcs_have_complex)
        throw(_presplit_bc_error())
    end

    if !typed_dvs && any(bc -> bc isa Equation && (_is_false_constant(bc.lhs) || _is_false_constant(bc.rhs)), bcs_flat)
        throw(
            ArgumentError(
                "A boundary condition reduced to `false` before PDEBase could inspect its complex value. " *
                    "Symbolics may truncate complex constants in `~`; declare the dependent variable as `::Complex`."
            )
        )
    end

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
