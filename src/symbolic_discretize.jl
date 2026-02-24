function cardinalize_eqs!(pdesys)
    pdeeqs = get_eqs(pdesys)
    for (i, eq) in enumerate(pdeeqs)
        pdeeqs[i] = eq.lhs - eq.rhs ~ 0
    end
    return
end

function SciMLBase.symbolic_discretize(pdesys::PDESystem, discretization::AbstractDiscretization; checks=true)
    t = get_time(discretization)
    pdesys, complexmap = handle_complex(pdesys)
    cardinalize_eqs!(pdesys)
    pdesys, replaced_vars = make_pdesys_compatible(pdesys)

    ############################
    # System Parsing and Transformation
    ############################
    # Parse the variables in to the right form and store useful information about the system
    v = VariableMap(pdesys, discretization, replaced_vars = replaced_vars)
    # Check for basic interface errors
    interface_errors(pdesys, v, discretization)
    # Extract tspan
    tspan = t !== nothing ? v.intervals[t] : nothing
    # Find the derivative orders in the bcs
    bcorders = Dict(map(x -> x => d_orders(x, get_bcs(pdesys)), all_ivs(v)))
    # Create a map of each variable to their boundary conditions including initial conditions
    boundarymap = parse_bcs(get_bcs(pdesys), v, bcorders)
    # Check that the boundary map is valid
    check_boundarymap(boundarymap, v, discretization)

    # Transform system so that it is compatible with the discretization
    if should_transform(pdesys, discretization, boundarymap)
        pdesys = transform_pde_system!(v, boundarymap, pdesys, discretization)
    end

    pdeeqs = get_eqs(pdesys)
    bcs = get_bcs(pdesys)

    ############################
    # Discretization of system
    ############################
    disc_state = construct_disc_state(discretization)

    # Create discretized space and variables, this is called `s` throughout
    s = construct_discrete_space(v, discretization)
    # Get the interior and variable to solve for each equation
    #TODO: do the interiormap before and independent of the discretization i.e. `s`
    vareqmap = construct_var_equation_mapping(pdeeqs, boundarymap, s, discretization)
    # Get the derivative orders appearing in each equation
    pdeorders = Dict(map(x -> x => d_orders(x, pdeeqs), indvars(v)))
    bcorders = Dict(map(x -> x => d_orders(x, bcs), indvars(v)))
    orders = Dict(map(x -> x => collect(union(pdeorders[x], bcorders[x])), indvars(v)))

    # Generate finite difference weights
    derivweights = construct_differential_discretizer(pdesys, s, discretization, orders)

    # Separate bcs and ics
    ics = t === nothing ? [] :
        mapreduce(u -> boundarymap[u][t], vcat, operation.(depvars(v)))

    bcmap = Dict(
        map(collect(keys(boundarymap))) do u
            u => Dict(
                map(indvars(v)) do x
                    x => boundarymap[u][x]
                end
            )
        end
    )

    ####
    # Loop over equations, Discretizing them and their dependent variables' boundary conditions
    ####
    for pde in pdeeqs
        # Read the dependent variables on both sides of the equation
        depvars_lhs = get_depvars(pde.lhs, v.depvar_ops)
        depvars_rhs = get_depvars(pde.rhs, v.depvar_ops)
        depvars = collect(depvars_lhs ∪ depvars_rhs)
        depvars = filter(u -> !any(map(x -> unwrap_const(safe_unwrap(x)) isa Number, arguments(u))), depvars)

        eqvar = get_eqvar(vareqmap, pde)

        # * Assumes that all variables in the equation have same dimensionality except edgevals
        args = ivs(eqvar, v)
        indexmap = Dict([args[i] => i for i in 1:length(args)])
        # Generate the equations for the interior points
        discretize_equation!(
            disc_state, pde, vareqmap, eqvar, bcmap,
            depvars, s, derivweights, indexmap, discretization
        )
    end

    u0 = generate_ic_defaults(ics, s, discretization)

    # Combine PDE equations and BC equations
    # Pass u0 to generate_metadata for storage (needed for MTK v11 compatibility)
    metadata = generate_metadata(s, discretization, pdesys, boundarymap, complexmap, u0)

    return generate_system(disc_state, s, u0, tspan, metadata, discretization; checks=checks)
end
