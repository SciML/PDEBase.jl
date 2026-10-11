function cardinalize_eqs!(pdesys)
    pdeeqs = get_eqs(pdesys)
    for (i, eq) in enumerate(pdeeqs)
        pdeeqs[i] = eq.lhs - eq.rhs ~ 0
    end
    return
end

function SciMLBase.symbolic_discretize(
        pdesys::PDESystem, discretization::AbstractEquationSystemDiscretization; checks = true
    )
    interface_errors(pdesys, discretization)
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

    return generate_system(disc_state, s, u0, tspan, metadata, discretization; checks = checks)
end

"""
    symbolic_discretize(pdesys::PDESystem, discretization::AbstractOptimizationSystemDiscretization; checks = true)

Driver for discretizations that lower a `PDESystem` into an optimization system, such as
mesh-free methods that fit trial functions to the equations. Every equation and boundary
condition is lowered independently on its own residual domain: there is no notion of an
equation variable, boundary conditions are not required to lie on a boundary of the
domain, and no boundary map is parsed.

The hook sequence is

1. `handle_complex`, `cardinalize_eqs!` and `make_pdesys_compatible` normalize the system;
2. `VariableMap` and `interface_errors` analyze the variables;
3. `construct_disc_state(discretization)` creates the mutable discretization state;
4. `construct_discrete_space(v, pdesys, discretization)` creates the space the residuals
   are evaluated on;
5. `construct_differential_discretizer(pdesys, s, discretization, orders)` chooses how
   derivatives are lowered;
6. `discretize_equation!(disc_state, eq, kind, s, derivweights, discretization)` is
   called for every PDE (`kind = :pde`) and then every boundary condition (`kind = :bc`);
7. `generate_metadata(s, discretization, pdesys, nothing, complexmap, [])` and
   `generate_system(disc_state, s, nothing, nothing, metadata, discretization; checks)`
   assemble the result, which `generate_system` returns as is.
"""
function SciMLBase.symbolic_discretize(
        pdesys::PDESystem, discretization::AbstractOptimizationSystemDiscretization;
        checks = true
    )
    pdesys, complexmap = handle_complex(pdesys)
    cardinalize_eqs!(pdesys)
    pdesys, replaced_vars = make_pdesys_compatible(pdesys)

    v = VariableMap(pdesys, discretization, replaced_vars = replaced_vars)
    interface_errors(pdesys, v, discretization)

    pdeeqs = get_eqs(pdesys)
    bcs = get_bcs(pdesys)

    disc_state = construct_disc_state(discretization)
    s = construct_discrete_space(v, pdesys, discretization)
    orders = Dict(
        map(all_ivs(v)) do x
            x => collect(union(d_orders(x, pdeeqs), d_orders(x, bcs)))
        end
    )
    derivweights = construct_differential_discretizer(pdesys, s, discretization, orders)

    for pde in pdeeqs
        discretize_equation!(disc_state, pde, :pde, s, derivweights, discretization)
    end
    for bc in bcs
        discretize_equation!(disc_state, bc, :bc, s, derivweights, discretization)
    end

    metadata = generate_metadata(s, discretization, pdesys, nothing, complexmap, [])
    return generate_system(disc_state, s, nothing, nothing, metadata, discretization; checks)
end
