struct EquationState <: AbstractEquationSystemDiscretization
    eqs::Vector{Equation}
    bceqs::Vector{Equation}
end

function EquationState()
    return EquationState(Equation[], Equation[])
end

function _split_discrete_ic_key(x)
    x = unwrap(x)
    iscall(x) || return x, (), false

    op = operation(x)
    args = arguments(x)
    if op === getindex
        idx = Tuple(unwrap_const(safe_unwrap(i)) for i in args[2:end])
        return first(args), idx, true
    elseif op isa Symbolics.Differential && length(args) == 1
        arr, idx, isarr = _split_discrete_ic_key(first(args))
        return isarr ? (op(arr), idx, true) : (x, (), false)
    end
    return x, (), false
end

_discrete_ic_value(v) = unwrap_const(safe_unwrap(v))

function _time_derivative_order(x, t)
    x = unwrap(x)
    iscall(x) || return x, 0
    op = operation(x)
    if op isa Symbolics.Differential && isequal(op.x, t)
        return first(arguments(x)), op.order
    end
    return x, 0
end

# Array (slice-form) equations differentiate whole slices, e.g. D(u[2:n-1]); expand
# them to their elements so that they can be matched against elementwise IC keys.
function _expand_differential_var(var)
    SymbolicUtils.symtype(var) <: AbstractArray || return (var,)
    return vec(unwrap.(collect(Symbolics.scalarize(Symbolics.wrap(var)))))
end

function _discrete_initialization(eqs, t, u0)
    isempty(u0) && return Equation[], Dict{Any, Any}()
    t = unwrap(t)
    differential_orders = Dict{Any, Real}()
    for eq in eqs, derivative in Symbolics.get_differential_vars(eq)
        var, order = _time_derivative_order(derivative, t)
        order > 0 || continue
        for v in _expand_differential_var(var)
            differential_orders[v] = max(order, get(differential_orders, v, 0))
        end
    end

    init_eqs = Equation[]
    guesses = Dict{Any, Any}()
    for (key, value) in u0
        key = unwrap(key)
        value = _discrete_ic_value(value)
        var, order = _time_derivative_order(key, t)
        if order < get(differential_orders, var, 0)
            push!(init_eqs, key ~ value)
        end

        array, index, isindexed = _split_discrete_ic_key(key)
        if isindexed
            values = get!(() -> fill!(Array{Any}(undef, size(array)), NaN), guesses, array)
            values[index...] = value
        else
            guesses[key] = value
        end
    end
    return init_eqs, guesses
end

# Normalize an equation to `0 ~ ...` form for NonlinearSystem construction. Array
# (slice-form) equations cannot equate an array with a scalar zero, so subtract via
# broadcast and equate with a zero array of matching size.
# symtype falls back to typeof for non-symbolic values, so this covers literal arrays,
# symbolic arrays and scalars of either kind.
_is_array_valued(x) = SymbolicUtils.symtype(safe_unwrap(x)) <: AbstractArray

function _normalize_nonlinear_eq(eq::Equation)
    if _is_array_valued(eq.lhs) || _is_array_valued(eq.rhs)
        diff = Broadcast.materialize(
            Broadcast.broadcasted(-, Symbolics.wrap(eq.rhs), Symbolics.wrap(eq.lhs))
        )
        return zeros(size(diff)) ~ diff
    end
    return 0 ~ eq.rhs - eq.lhs
end

function generate_system(
        disc_state::EquationState, s, u0, tspan, metadata,
        disc::AbstractEquationSystemDiscretization;
        checks = true
    )
    discvars = get_discvars(s)
    t = get_time(disc)
    name = getfield(metadata.pdesys, :name)
    pdesys = metadata.pdesys
    alleqs = vcat(disc_state.eqs, unique(disc_state.bceqs))
    alldepvarsdisc = vec(reduce(vcat, vec(unique(reduce(vcat, vec.(values(discvars)))))))

    sys_defaults = Dict{Any, Any}(pdesys.initial_conditions)
    init_eqs = Equation[]
    guesses = Dict{Any, Any}()
    if t !== nothing && !isempty(u0)
        init_eqs, guesses = _discrete_initialization(alleqs, t, u0)
    end

    ps_raw = get_ps(pdesys)
    ps_raw = ps_raw === nothing || ps_raw === SciMLBase.NullParameters() ? Num[] : ps_raw
    # get_ps may return Pairs (e.g. [v => 0.5]); extract symbols and merge values into defaults
    if !isempty(ps_raw) && first(ps_raw) isa Pair
        ps = Num[first(p) for p in ps_raw]
        merge!(sys_defaults, Dict{Any, Any}(first(p) => last(p) for p in ps_raw))
    else
        ps = ps_raw
    end
    return try
        if t === nothing
            # At the time of writing, NonlinearProblems require that the system of equations be in this form:
            # 0 ~ ...
            # Thus, before creating a NonlinearSystem we normalize the equations s.t. the lhs is zero.
            eqs = map(_normalize_nonlinear_eq, alleqs)
            sys = System(
                eqs, alldepvarsdisc, ps, initial_conditions = sys_defaults, name = name,
                metadata = [ProblemTypeCtx => metadata], checks = checks
            )
            return sys, nothing
        else
            # * In the end we have reduced the problem to a system of equations in terms of Dt that can be solved by an ODE solver.

            sys = System(
                alleqs, t, alldepvarsdisc, ps;
                initial_conditions = sys_defaults,
                initialization_eqs = init_eqs,
                guesses = guesses,
                name = name,
                metadata = [ProblemTypeCtx => metadata], checks = checks
            )
            return sys, tspan
        end
    catch e
        println("The system of equations is:")
        println("Number of equations: ", length(alleqs))
        for (i, eq) in enumerate(alleqs)
            println("Eq $i: $eq")
        end
        println()
        println("Discretization failed, please post an issue on https://github.com/SciML/MethodOfLines.jl with the failing code and system at low point count.")
        println()
        rethrow(e)
    end
end

function SciMLBase.discretize(
        pdesys::PDESystem,
        discretization::AbstractEquationSystemDiscretization;
        analytic = nothing, checks = true, kwargs...
    )
    sys, tspan = SciMLBase.symbolic_discretize(pdesys, discretization; checks = checks)
    return try
        simpsys = mtkcompile(sys)
        if tspan === nothing
            add_metadata!(getmetadata(sys, ProblemTypeCtx, nothing), sys)
            # MTK v11 requires symbolic map for initial guess
            unknowns_list = ModelingToolkit.unknowns(simpsys)
            u0_guess = Dict(u => 1.0 for u in unknowns_list)
            return prob = NonlinearProblem(
                simpsys, u0_guess;
                discretization.kwargs..., kwargs...
            )
        else
            mol_metadata = getmetadata(simpsys, ProblemTypeCtx, nothing)
            add_metadata!(mol_metadata, sys)
            prob = ODEProblem(
                simpsys, nothing, tspan;
                discretization.kwargs...,
                kwargs...
            )
            if analytic === nothing
                return prob
            else
                f = ODEFunction(
                    pdesys, discretization, analytic = analytic,
                    discretization.kwargs..., kwargs...
                )

                return ODEProblem(
                    f, prob.u0, prob.tspan, prob.p;
                    discretization.kwargs..., kwargs...
                )
            end
        end
    catch e
        error_analysis(sys, e)
    end
end

function error_analysis(sys::System, e)
    is_time_dependent(sys) || rethrow(e)
    eqs = get_eqs(sys)
    unknowns = get_unknowns(sys)
    t = get_iv(sys)
    println("The system of equations is:")
    println(eqs)
    return if e isa ModelingToolkit.ExtraVariablesSystemException
        rs = [Differential(t)(state) => state for state in unknowns]
        extraunknowns = [state for state in unknowns]
        extraeqs = [eq for eq in eqs]
        numderivs = 0
        for r in rs
            for eq in extraeqs
                if subsmatch(eq.lhs, r) | subsmatch(eq.rhs, r)
                    extraunknowns = vec(setdiff(extraunknowns, [r.second]))
                    extraeqs = vec(setdiff(extraeqs, [eq]))
                    numderivs += 1
                    break
                end
            end
        end
        println()
        println("There are $(length(unknowns)) variables and $(length(eqs)) equations.\n")
        println("There are $numderivs time derivatives.\n")
        println("The variables without time derivatives are:")
        println(extraunknowns)
        println()
        println("The equations without time derivatives are:")
        println(extraeqs)
        rethrow(e)
    else
        rethrow(e)
    end
end
