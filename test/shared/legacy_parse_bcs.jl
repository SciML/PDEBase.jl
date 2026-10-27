# A verbatim copy of `parse_bcs(bcs, v, orders)` and of the dependent-variable collection
# of the boundary constructors as they were before `parse_bcs` dispatched on the
# discretization, kept as the oracle for the equivalence tests in test/parse_bcs.jl.

using PDEBase: VariableMap, generate_boundary_matching_rules, split_terms, flatten_vardict,
    subsmatch, LowerBoundary, UpperBoundary, InterfaceBoundary, HigherOrderInterfaceBoundary,
    get_depvars, depvar, all_ivs
using SymbolicUtils: operation, arguments

function legacy_boundary_vars(eq, t, v)
    depvars_lhs = get_depvars(eq.lhs, v.depvar_ops)
    depvars_rhs = get_depvars(eq.rhs, v.depvar_ops)
    depvars = collect(depvars_lhs ∪ depvars_rhs)

    allx̄ = Set(
        filter(
            !isempty, map(u -> filter(x -> t === nothing || !isequal(x, t), arguments(u)), depvars)
        )
    )
    return depvar.(depvars, [v]), first(allx̄)
end

function legacy_parse_bcs(bcs, v::VariableMap, orders)
    t = v.time
    depvar_ops = v.depvar_ops

    # Create some rules to match which bundary/variable a bc concerns
    u0 = []
    bceqs = []
    ## BC matching rules, returns the variable and parameter the bc concerns

    lower_boundary_rules, upper_boundary_rules = generate_boundary_matching_rules(v, orders)

    boundarymap = Dict([operation(u) => Dict([x => [] for x in all_ivs(v)]) for u in v.ū])

    # Generate initial conditions and bc equations
    for bc in bcs
        boundary = nothing
        # Split out additive terms
        terms = split_terms(bc)
        # * Assume that the BC is defined on the edge of the domain
        # Check whether the bc is on the lower boundary, or periodic, we don't care which depvar/var
        u_ = u__ = x_ = x__ = nothing
        isinterface = false
        interface_orders = []
        for term in terms, r in flatten_vardict(lower_boundary_rules)
            # Need to check the order, higher order should be discretized as a normal upper index
            # As zeroth order goes on the lower index. Needs a zeroth order interface to make sense,
            # but I assume that this will always be the case where there are higher order interfaces.
            #Check if the rule changes the expression
            if subsmatch(term, r)
                # Get the matched variables from the rule
                u_, x_, order = r.second
                # Mark the boundary
                if boundary === nothing
                    boundary = LowerBoundary(u_, t, x_, order, bc, v)
                end
                # do it again for the upper end to check for periodic, but only check the current depvar and indvar
                push!(interface_orders, order)
                for term_ in setdiff(terms, [term]),
                        r_ in flatten_vardict(upper_boundary_rules)

                    if subsmatch(term_, r_)
                        isinterface = true
                        u__, x__, order__ = r_.second
                        @assert ndims(u_, v) == ndims(u__, v) "Invalid Interface Boundary $bc: Dependent variables $(u_) and $(u__) have different numbers of dimensions."

                        push!(interface_orders, order__)
                    end
                end
                # Handle flux condition at interface
            end
        end
        # If this is a condition on flux at an interface, leave it to become an upper boundary
        # Else if it is a simple interface, make interface boundaries
        if isinterface
            if all(==(0), interface_orders)
                boundary = (
                    InterfaceBoundary{Val(false), Val(true)}(u_, u__, x_, x__, bc),
                    InterfaceBoundary{Val(true), Val(false)}(u__, u_, x__, x_, bc),
                )
            else
                boundary = HigherOrderInterfaceBoundary(
                    u__, u_, x__, x_, t, maximum(interface_orders), bc, v
                )
            end
        end
        # repeat for upper boundary
        if boundary === nothing
            for term in terms, r in flatten_vardict(upper_boundary_rules)

                if subsmatch(term, r)
                    u_, x_, order = r.second
                    boundary = UpperBoundary(u_, t, x_, order, bc, v)
                    break
                end
            end
        end
        @assert boundary !== nothing "Boundary condition $bc is not on a boundary of the domain, or is not a valid boundary condition"
        if boundary isa Tuple
            push!(boundarymap[operation(boundary[1].u)][boundary[1].x], boundary[1])
            push!(boundarymap[operation(boundary[2].u)][boundary[2].x], boundary[2])
        else
            push!(boundarymap[operation(boundary.u)][boundary.x], boundary)
        end
    end
    return boundarymap
end
