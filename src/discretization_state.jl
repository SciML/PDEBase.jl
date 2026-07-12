struct EquationState <: AbstractEquationSystemDiscretization
    eqs::Vector{Equation}
    bceqs::Vector{Equation}
end

function EquationState()
    return EquationState(Equation[], Equation[])
end

"""
    discrete_u0_to_atomic_map(u0)

Convert indexed discrete IC pairs `(u(t))[i] => val` into a map compatible with
ModelingToolkit v11 `AtomicArrayDict` (whole-array keys only, e.g. `u(t) => values`).
Used for `guesses` / non-indexed maps. Empty input returns an empty `Dict`.
"""
function discrete_u0_to_atomic_map(u0)
    isempty(u0) && return Dict{Any, Any}()
    groups = Dict{Any, Dict{Any, Any}}()
    scalars = Dict{Any, Any}()
    for (k, v) in u0
        k = unwrap(k)
        arr, isarr = split_indexed_var(k)
        if isarr
            buf = get!(() -> Dict{Any, Any}(), groups, arr)
            buf[get_stable_index(k)] = v
        else
            scalars[k] = v
        end
    end
    out = Dict{Any, Any}()
    for (arr, idxs) in groups
        vals = Array{Float64}(undef, size(arr))
        fill!(vals, NaN)
        for (idx, v) in idxs
            vals[idx] = Float64(try
                Symbolics.value(v)
            catch
                v
            end)
        end
        out[arr] = vals
    end
    return merge!(out, scalars)
end

_numeric_ic_value(v) = Float64(try
    Symbolics.value(v)
catch
    v
end)

"""
    variables_with_time_derivative(eqs, t)

Return discrete unknowns that appear under `Differential(t)` in `eqs`. These are
the differential states whose initial values should be hard initialization
equations; remaining (algebraic) discrete states should only be guesses.
"""
function variables_with_time_derivative(eqs, t)
    found = Any[]
    t = unwrap(t)
    function walk(ex)
        ex = unwrap(ex)
        if iscall(ex)
            op = operation(ex)
            args = arguments(ex)
            if op isa Differential && isequal(op.x, t)
                push!(found, args[1])
            end
            foreach(walk, args)
        end
        return nothing
    end
    for eq in eqs
        walk(eq.lhs)
        walk(eq.rhs)
    end
    return unique(found)
end

"""
    discrete_initialization(eqs, t, u0)

Build MTK initialization data from discrete IC pairs:

- **Hard** `initialization_eqs`: `var ~ val` only for variables that appear under
  `Differential(t)` (differential states). This matches the continuous IC
  structure and avoids over-determining algebraic residuals.
- **Guesses**: all discrete ICs as whole-array maps (algebraic nodes use the
  projected continuous IC as a starting guess for the init solve).

Returns `(initialization_eqs, guesses)`.
"""
function discrete_initialization(eqs, t, u0)
    isempty(u0) && return Equation[], Dict{Any, Any}()
    u0_dict = Dict{Any, Any}(unwrap(k) => v for (k, v) in u0)
    init_eqs = Equation[]
    for dv in variables_with_time_derivative(eqs, t)
        dv_u = unwrap(dv)
        if haskey(u0_dict, dv_u)
            push!(init_eqs, dv ~ _numeric_ic_value(u0_dict[dv_u]))
        else
            # try isequal match (identity of keys can differ by wrapping)
            for (k, v) in u0_dict
                if isequal(k, dv_u)
                    push!(init_eqs, dv ~ _numeric_ic_value(v))
                    break
                end
            end
        end
    end
    guesses = discrete_u0_to_atomic_map(u0)
    return init_eqs, guesses
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

    # Continuous/non-indexed defaults from the PDESystem (parameters, etc.)
    sys_defaults = Dict{Any, Any}(pdesys.initial_conditions)
    # Discrete grid ICs: hard init eqs for differential states + guesses for all
    # (including algebraic). Do not dump indexed keys into AtomicArrayDict ICs.
    init_eqs = Equation[]
    guesses = Dict{Any, Any}()
    if t !== nothing && !isempty(u0)
        init_eqs, guesses = discrete_initialization(alleqs, t, u0)
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
            eqs = map(eq -> 0 ~ eq.rhs - eq.lhs, alleqs)
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
            # Discrete ICs are already on the System as initialization_eqs (differential
            # states) and guesses (all grid values). Do not re-pass indexed u0 pairs as
            # the ODEProblem operating point: that forces every grid value as a hard IC
            # and over-determines algebraic residuals. Only pass parameters here.
            pdesys_ic = mol_metadata.pdesys.initial_conditions
            ps_raw = get_ps(mol_metadata.pdesys)
            param_vals = Dict{Any, Any}()
            if ps_raw !== nothing && ps_raw !== SciMLBase.NullParameters() && !isempty(ps_raw)
                # get_ps may return Pairs (e.g. [v => 0.5]); extract parameter values
                if first(ps_raw) isa Pair
                    for p in ps_raw
                        param_vals[first(p)] = last(p)
                    end
                else
                    # Fall back to looking up parameters in initial_conditions
                    ps_unwrapped = [safe_unwrap(p) for p in ps_raw]
                    for (k, v) in pairs(pdesys_ic)
                        k_unwrapped = safe_unwrap(k)
                        if any(p -> isequal(k_unwrapped, safe_unwrap(p)), ps_unwrapped)
                            v_numeric = try
                                Symbolics.value(v)
                            catch
                                safe_unwrap(v)
                            end
                            param_vals[k] = v_numeric
                        end
                    end
                end
            end
            # Enable MTK initialization (DefaultInit -> OverrideInit).
            # Callers may override build_initializeprob / op via kwargs.
            op = isempty(param_vals) ? nothing : param_vals
            prob = ODEProblem(
                simpsys, op, tspan; build_initializeprob = true,
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
