struct EquationState <: AbstractEquationSystemDiscretization
    eqs::Vector{Equation}
    bceqs::Vector{Equation}
end

function EquationState()
    return EquationState(Equation[], Equation[])
end

function generate_system(disc_state::EquationState, s, u0, tspan, metadata,
        disc::AbstractEquationSystemDiscretization)
    discvars = get_discvars(s)
    t = get_time(disc)
    name = getfield(metadata.pdesys, :name)
    pdesys = metadata.pdesys
    alleqs = vcat(disc_state.eqs, unique(disc_state.bceqs))
    alldepvarsdisc = vec(reduce(vcat, vec(unique(reduce(vcat, vec.(values(discvars)))))))

    defaults = merge(ModelingToolkit.defaults(pdesys), Dict(u0))
    ps = get_ps(pdesys)
    ps = ps === nothing || ps === SciMLBase.NullParameters() ? Num[] : ps
    # Finalize
    # if haskey(metadata.disc.kwargs, :checks)
    #     checks = metadata.disc.kwargs[:checks]
    # else
    checks = true
    # end
    try
        if t === nothing
            # At the time of writing, NonlinearProblems require that the system of equations be in this form:
            # 0 ~ ...
            # Thus, before creating a NonlinearSystem we normalize the equations s.t. the lhs is zero.
            eqs = map(eq -> 0 ~ eq.rhs - eq.lhs, alleqs)
            sys = System(eqs, alldepvarsdisc, ps, defaults = defaults, name = name,
                metadata = [ProblemTypeCtx => metadata], checks = checks)
            return sys, nothing
        else
            # * In the end we have reduced the problem to a system of equations in terms of Dt that can be solved by an ODE solver.

            sys = System(alleqs, t, alldepvarsdisc, ps, defaults = defaults, name = name,
                metadata = [ProblemTypeCtx => metadata], checks = checks)
            return sys, tspan
        end
    catch e
        println("The system of equations is:")
        #println(alleqs)
        println()
        println("Discretization failed, please post an issue on https://github.com/SciML/MethodOfLines.jl with the failing code and system at low point count.")
        println()
        rethrow(e)
    end
end

function SciMLBase.discretize(pdesys::PDESystem,
        discretization::AbstractEquationSystemDiscretization;
        analytic = nothing, kwargs...)
    sys, tspan = SciMLBase.symbolic_discretize(pdesys, discretization)
    try
        simpsys = mtkcompile(sys)
        if tspan === nothing
            add_metadata!(getmetadata(sys, ProblemTypeCtx, nothing), sys)
            return prob = NonlinearProblem(simpsys, ones(length(get_eqs(simpsys)));
                discretization.kwargs..., kwargs...)
        else
            add_metadata!(getmetadata(simpsys, ProblemTypeCtx, nothing), sys)
            prob = ODEProblem(simpsys, Pair[], tspan; build_initializeprob = false,
                discretization.kwargs...,
                kwargs...)
            if analytic === nothing
                return prob
            else
                f = ODEFunction(pdesys, discretization, analytic = analytic,
                    discretization.kwargs..., kwargs...)

                return ODEProblem(f, prob.u0, prob.tspan, prob.p;
                    discretization.kwargs..., kwargs...)
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
    if e isa ModelingToolkit.ExtraVariablesSystemException
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
