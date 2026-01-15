struct EquationState <: AbstractEquationSystemDiscretization
    eqs::Vector{Equation}
    bceqs::Vector{Equation}
end

function EquationState()
    return EquationState(Equation[], Equation[])
end

function generate_system(
        disc_state::EquationState, s, u0, tspan, metadata,
        disc::AbstractEquationSystemDiscretization
    )
    discvars = get_discvars(s)
    t = get_time(disc)
    name = getfield(metadata.pdesys, :name)
    pdesys = metadata.pdesys
    alleqs = vcat(disc_state.eqs, unique(disc_state.bceqs))
    alldepvarsdisc = vec(reduce(vcat, vec(unique(reduce(vcat, vec.(values(discvars)))))))

    # u0 is now stored in metadata (passed during metadata construction)
    # MTK v11's AtomicArrayDict doesn't allow indexed array symbolics as keys in initial_conditions
    # Only pass non-indexed initial conditions to System
    sys_defaults = Dict(pdesys.initial_conditions)

    ps = get_ps(pdesys)
    ps = ps === nothing || ps === SciMLBase.NullParameters() ? Num[] : ps
    # Finalize
    # if haskey(metadata.disc.kwargs, :checks)
    #     checks = metadata.disc.kwargs[:checks]
    # else
    checks = true
    # end
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
                alleqs, t, alldepvarsdisc, ps, initial_conditions = sys_defaults, name = name,
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
        analytic = nothing, kwargs...
    )
    sys, tspan = SciMLBase.symbolic_discretize(pdesys, discretization)
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
            # Get u0 from metadata (stored there for MTK v11 compatibility)
            u0 = hasproperty(mol_metadata, :u0) ? mol_metadata.u0 : []
            # Get parameter values from the original pdesys initial_conditions
            # MTK v11 needs parameter values passed explicitly when creating ODEProblem
            pdesys_ic = mol_metadata.pdesys.initial_conditions
            ps = get_ps(mol_metadata.pdesys)
            if ps !== nothing && ps !== SciMLBase.NullParameters() && !isempty(ps)
                # Extract only parameter values (not unknown initial conditions)
                # Use isequal for proper symbolic comparison
                ps_unwrapped = [safe_unwrap(p) for p in ps]
                param_vals = Dict{Any,Any}()
                for (k, v) in pairs(pdesys_ic)
                    k_unwrapped = safe_unwrap(k)
                    if any(p -> isequal(k_unwrapped, safe_unwrap(p)), ps_unwrapped)
                        # Extract numeric value from symbolic constant if needed
                        v_numeric = try
                            Symbolics.value(v)
                        catch
                            safe_unwrap(v)
                        end
                        param_vals[k] = v_numeric
                    end
                end
                if !isempty(param_vals)
                    # MTK v11 API: merge u0 and params into a single Dict
                    op = merge(Dict(u0), param_vals)
                    prob = ODEProblem(
                        simpsys, op, tspan; build_initializeprob = false,
                        discretization.kwargs...,
                        kwargs...
                    )
                else
                    prob = ODEProblem(
                        simpsys, u0, tspan; build_initializeprob = false,
                        discretization.kwargs...,
                        kwargs...
                    )
                end
            else
                prob = ODEProblem(
                    simpsys, u0, tspan; build_initializeprob = false,
                    discretization.kwargs...,
                    kwargs...
                )
            end
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
