"""
    PeriodicMap{hasperiodic}

A structure that tracks periodic boundary conditions for variables in a PDE system.

# Type Parameters
- `hasperiodic`: A `Val` type parameter indicating whether any periodic conditions exist

# Fields
- `map`: A nested dictionary structure mapping dependent variable operations to their
  independent variables and whether each boundary is periodic

# Constructor
    PeriodicMap(bmap, v::VariableMap)

Constructs a `PeriodicMap` from a boundary map and variable map, automatically detecting
which boundaries have periodic conditions.
"""
struct PeriodicMap{hasperiodic}
    map::Any
end

function PeriodicMap(bmap, v)
    map = Dict(
        [
            operation(u) => Dict([x => isperiodic(bmap, u, x) for x in all_ivs(v)])
                for u in v.ū
        ]
    )
    vals = reduce(vcat, collect.(values.(collect(values(map)))))
    hasperiodic = Val(any(p -> p isa Val{true}, vals))
    return PeriodicMap{hasperiodic}(map)
end
