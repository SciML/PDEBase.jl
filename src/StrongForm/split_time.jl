"Additive terms of an expanded `ex`; `expand` leaves `(a + b) / d` whole, so quotients are split too."
function addterms(ex)
    iscall(ex) || return [ex]
    operation(ex) === (+) && return reduce(vcat, addterms.(arguments(ex)))
    if operation(ex) === (/)
        num, den = arguments(ex)
        return [term / den for term in addterms(num)]
    end
    return [ex]
end
factors(ex) = iscall(ex) && operation(ex) === (*) ? collect(arguments(ex)) :
    iscall(ex) && operation(ex) === (/) ? vcat(factors(arguments(ex)[1]), [1 / arguments(ex)[2]]) : [ex]

"Total order in `t` of the deepest derivative chain in `ex`; reads `Differential.order` since Symbolics folds `Dt(Dt(u))`."
function time_order(ex, t)
    iscall(ex) || return 0
    op = operation(ex)
    own = op isa Differential && isequal(op.x, t) ? op.order : 0
    return own + maximum((time_order(a, t) for a in arguments(ex)); init = 0)
end

"""
    split_time(eq, depvars, t) -> (mass, spatial)

Stage a of the strong-form analysis. Expands the residual `e = lhs - rhs` of `eq` (an
`Equation` or a residual expression) and splits it into `mass[f]`, the coefficient of `∂ₜu_f`
for each `u_f` in `depvars`, and `spatial`, the vector of remaining additive terms (with
derivatives distributed over sums, see [`distribute_derivatives`](@ref)), so that
`e = Σ_f mass[f]⋅∂ₜu_f + sum(spatial)`. `t` is `nothing` for a steady problem.

Errors on a term that is second order in time, not of the form coefficient × `∂ₜu_f`
(nonlinear in `∂ₜu` or `∂ₜ` inside a spatial derivative), or whose coefficient depends on
the state.
"""
split_time(eq::Equation, depvars, t) = split_time(eq.lhs - eq.rhs, depvars, t)
function split_time(ex, depvars, t)
    mass = Num[0 for _ in depvars]
    spatial = Num[]
    for term in addterms(unwrap(Symbolics.expand(distribute_derivatives(Symbolics.expand(ex)))))
        order = t === nothing ? 0 : time_order(term, t)
        if order == 0
            push!(spatial, term)
            continue
        end
        order > 1 && error("$term is second order in time; write as a first-order system")
        fs = factors(term)
        f = findfirst(u -> any(isequal(unwrap(Differential(t)(u))), fs), depvars)
        f === nothing && error("time derivative in $term is not ∂ₜu × coefficient")
        dudt = unwrap(Differential(t)(depvars[f]))
        m = prod(fs[findall(!isequal(dudt), fs)]; init = 1)
        time_order(m, t) == 0 && count(isequal(dudt), fs) == 1 ||
            error("time derivative in $term is not ∂ₜu × coefficient")
        is_state(m, depvars) && error("state-dependent mass in $term needs residual-form DAE")
        mass[f] += m
    end
    return mass, spatial
end
