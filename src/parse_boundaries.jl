### INITIAL AND BOUNDARY CONDITIONS ###

#TODO: Retire DiscreteSpace and move to a DiscreteEquation class with corresponding DiscreteVariables which carry all the required information, define substitute methods for these.

"""
    AbstractBoundary

Abstract base type for all boundary condition representations.
"""
abstract type AbstractBoundary end

"""
    AbstractTruncatingBoundary <: AbstractBoundary

Abstract type for boundary conditions at domain boundaries that truncate the domain.
"""
abstract type AbstractTruncatingBoundary <: AbstractBoundary end

"""
    AbstractInterfaceBoundary <: AbstractTruncatingBoundary

Abstract type for boundary conditions at interfaces between different regions or subdomains.
"""
abstract type AbstractInterfaceBoundary <: AbstractTruncatingBoundary end

# The dependent variables of boundary equation `eq`, in their declared form, and the
# non-time arguments of one of its dependent-variable calls, as stored in a boundary.
function _boundary_vars(eq, t, v)
    depvars_lhs = get_depvars(eq.lhs, v.depvar_ops)
    depvars_rhs = get_depvars(eq.rhs, v.depvar_ops)
    depvars = collect(depvars_lhs ∪ depvars_rhs)
    #depvars =  filter(u -> !any(map(x-> x isa Number, arguments(u))), depvars)

    allx̄ = Set(
        filter(
            !isempty, map(u -> filter(x -> t === nothing || !isequal(x, t), arguments(u)), depvars)
        )
    )
    return depvar.(depvars, [v]), first(allx̄)
end

"""
    LowerBoundary <: AbstractTruncatingBoundary

Represents a boundary condition at the lower boundary of a domain.

The public fields describe the boundary after PDEBase has matched the input
condition to a dependent variable and spatial coordinate. Construct instances
through `LowerBoundary(u, t, x, order, eq, v)` while parsing a `PDESystem`.

# Fields
- `u`: The dependent variable
- `x`: The independent variable at this boundary
- `depvars`: All dependent variables involved in the boundary condition
- `indvars`: All independent variables involved in the boundary condition
- `eq`: The boundary condition equation
- `order`: The order of the derivative in the boundary condition

# Arguments
- `u`: The dependent-variable expression at the boundary.
- `t`: The time variable, or `nothing` for a steady-state system.
- `x`: The independent variable whose lower endpoint is represented.
- `order`: The derivative order associated with the boundary condition.
- `eq`: The symbolic boundary equation.
- `v::VariableMap`: The variable map used to reconstruct dependent variables.

# Returns
- `LowerBoundary`: A normalized boundary representation.

# Examples
```julia
boundary = LowerBoundary(u(t, 0), t, x, 0, u(t, 0) ~ 0, v)
boundary.x == x
```
"""
struct LowerBoundary <: AbstractTruncatingBoundary
    u::Any
    x::Any
    depvars::Any
    indvars::Any
    eq::Any
    order::Any
    function LowerBoundary(u, t, x, order, eq, v)
        depvars, indvars = _boundary_vars(eq, t, v)
        return new(u, x, depvars, indvars, eq, order)
    end
end

"""
    UpperBoundary <: AbstractTruncatingBoundary

Represents a boundary condition at the upper boundary of a domain.

Construct instances through `UpperBoundary(u, t, x, order, eq, v)` while
parsing a `PDESystem`.

# Fields
- `u`: The dependent variable
- `x`: The independent variable at this boundary
- `depvars`: All dependent variables involved in the boundary condition
- `indvars`: All independent variables involved in the boundary condition
- `eq`: The boundary condition equation
- `order`: The order of the derivative in the boundary condition

# Arguments
- `u`: The dependent-variable expression at the boundary.
- `t`: The time variable, or `nothing` for a steady-state system.
- `x`: The independent variable whose upper endpoint is represented.
- `order`: The derivative order associated with the boundary condition.
- `eq`: The symbolic boundary equation.
- `v::VariableMap`: The variable map used to reconstruct dependent variables.

# Returns
- `UpperBoundary`: A normalized boundary representation.

# Examples
```julia
boundary = UpperBoundary(u(t, 1), t, x, 0, u(t, 1) ~ 0, v)
boundary.x == x
```
"""
struct UpperBoundary <: AbstractTruncatingBoundary
    u::Any
    x::Any
    depvars::Any
    indvars::Any
    eq::Any
    order::Any
    function UpperBoundary(u, t, x, order, eq, v)
        depvars, indvars = _boundary_vars(eq, t, v)
        return new(u, x, depvars, indvars, eq, order)
    end
end

# Note that it is assumed throughout MOL that the variables in an interface BC have the same argument signature,
# differing in only one variable which is that of the interface. This is not checked here, but will cause errors if it is not true.
# Interfaces are assumed to be on the lower boundary of the domain.

"""
    InterfaceBoundary{IsUpper_u, IsUpper_u2} <: AbstractInterfaceBoundary

Represents a boundary condition at an interface between two regions or subdomains.

Use the ordinary five-field constructor after the two sides have been
classified. The type parameters record whether each side is at its upper
endpoint; they are `Val` types so dispatch can use the orientation.

# Type Parameters
- `IsUpper_u`: `Val{true}` if the first variable is at the upper boundary, `Val{false}` otherwise
- `IsUpper_u2`: `Val{true}` if the second variable is at the upper boundary, `Val{false}` otherwise

# Fields
- `u`: The first dependent variable at the interface
- `u2`: The second dependent variable at the interface
- `x`: The independent variable for the first variable
- `x2`: The independent variable for the second variable
- `eq`: The interface boundary condition equation

# Note
It is assumed that the variables in an interface boundary condition have the same argument
signature, differing in only one variable which defines the interface. This assumption is
not validated but will cause errors if violated.

# Arguments
- `u`, `u2`: Dependent-variable expressions on the two sides.
- `x`, `x2`: The interface coordinates on the two sides.
- `eq`: The symbolic interface equation.

# Returns
- `InterfaceBoundary`: An orientation-aware interface condition.

# Examples
```julia
boundary = InterfaceBoundary{Val(true), Val(false)}(u, u2, x, x2, u ~ u2)
isinterface(boundary) == Val(true)
```
"""
struct InterfaceBoundary{IsUpper_u, IsUpper_u2} <: AbstractInterfaceBoundary
    u::Any
    u2::Any
    x::Any
    x2::Any
    eq::Any
end

"""
    HigherOrderInterfaceBoundary <: AbstractInterfaceBoundary

Represents a higher-order boundary condition at an interface between regions.

Higher-order interface boundary conditions involve derivatives and are handled differently
from simple interface conditions. These are typically flux conditions or derivative matching
conditions at interfaces.

Construct instances through `HigherOrderInterfaceBoundary(u, u2, x, x2, t,
order, eq, v)` so the dependent and independent variable fields are normalized
from `eq`.

# Fields
- `u`: The first dependent variable at the interface
- `u2`: The second dependent variable at the interface
- `x`: The independent variable for the first variable
- `x2`: The independent variable for the second variable
- `depvars`: All dependent variables involved in the boundary condition
- `indvars`: All independent variables involved in the boundary condition
- `eq`: The interface boundary condition equation
- `order`: The maximum order of derivatives in the boundary condition

# Arguments
- `u`, `u2`: Dependent-variable expressions on the two sides.
- `x`, `x2`: The interface coordinates on the two sides.
- `t`: The time variable, or `nothing` for a steady-state system.
- `order`: The maximum derivative order represented by `eq`.
- `eq`: The symbolic interface equation.
- `v::VariableMap`: The variable map used to reconstruct dependent variables.

# Returns
- `HigherOrderInterfaceBoundary`: A normalized higher-order interface.
"""
struct HigherOrderInterfaceBoundary <: AbstractInterfaceBoundary
    u::Any
    u2::Any
    x::Any
    x2::Any
    depvars::Any
    indvars::Any
    eq::Any
    order::Any
    function HigherOrderInterfaceBoundary(u, u2, x, x2, t, order, eq, v)
        depvars, indvars = _boundary_vars(eq, t, v)
        return new(u, u2, x, x2, depvars, indvars, eq, order)
    end
end

"""
    PieceBoundary <: AbstractBoundary

Represents a boundary condition on a named [`BoundaryPiece`](@ref) of the spatial
domain, written with traces such as `u(t, Γ) ~ g` or `Dn_Γ(u(t, Γ)) ~ h`.

`parse_bcs` constructs these for an `AbstractMeshDiscretization`, keyed by the piece in
the boundary map. Construct instances through `PieceBoundary(u, piece, eq, v)`.

# Fields
- `u`: The first trace in the condition, such as `u(t, Γ)`
- `piece`: The boundary piece the condition is on
- `depvars`: All dependent variables involved in the boundary condition
- `eq`: The boundary condition equation
- `kind`: The syntactic kind of the condition from `PDEBase.bckind`:
  `PDEBase.Essential()`, `PDEBase.Natural()` or `PDEBase.Robin()`
- `order`: The highest derivative order applied to a trace, where a
  [`NormalDerivative`](@ref) counts one

# Arguments
- `u`: The trace the condition is filed under.
- `piece::BoundaryPiece`: The piece of the trace.
- `eq::Equation`: The symbolic boundary equation.
- `v::VariableMap`: The variable map used to reconstruct dependent variables.

# Returns
- `PieceBoundary`: A normalized boundary representation.

# Examples
```julia
Γ = boundarypiece(Ω, :wall)
boundary = PieceBoundary(u(t, Γ), Γ, NormalDerivative(Γ)(u(t, Γ)) ~ 0, v)
boundary.kind == PDEBase.Natural()
boundary.order == 1
```
"""
struct PieceBoundary <: AbstractBoundary
    u::Any
    piece::BoundaryPiece
    depvars::Any
    eq::Equation
    kind::BCKind
    order::Any
    function PieceBoundary(u, piece, eq, v)
        depvars, _ = _boundary_vars(eq, v.time, v)
        return new(u, piece, depvars, eq, bckind(eq), max(_trace_order(eq.lhs), _trace_order(eq.rhs)))
    end
end

# The highest spatial derivative order applied to a trace in `ex`. As in `bckind`, a
# `Differential` in a variable the trace keeps, such as `t` in `Dt(u(t, Γ))`, does not count.
function _trace_order(ex, derivatives = ())
    ex = safe_unwrap(ex)
    iscall(ex) || return 0
    args = arguments(ex)
    if any(_ispiece, args)
        spatial = filter(op -> _is_spatial_derivative(op, args), collect(derivatives))
        return sum(op -> op isa NormalDerivative ? 1 : op.order, spatial; init = 0)
    end
    op = operation(ex)
    inner = op isa Union{Differential, NormalDerivative} ? (derivatives..., op) : derivatives
    return maximum(a -> _trace_order(a, inner), args; init = 0)
end

function Base.isequal(i1::InterfaceBoundary, i2::InterfaceBoundary)
    front = (
        isequal(i1.u, i2.u) & isequal(i1.u2, i2.u2) & isequal(i1.x, i2.x) &
            isequal(i1.x2, i2.x2)
    )
    back = (
        isequal(i1.u, i2.u2) & isequal(i1.u2, i2.u) & isequal(i1.x, i2.x2) &
            isequal(i1.x2, i2.x)
    )
    return front | back
end

"""
    getvars(b::AbstractBoundary)

Return the dependent variable and independent boundary variable represented by
boundary object `b`.

# Arguments
- `b::AbstractBoundary`: A lower, upper, or interface boundary.

# Returns
- `Tuple`: `(b.u, b.x)`.

# Examples
```julia
u_at_boundary, x = getvars(boundary)
```
"""
getvars(b::AbstractBoundary) = (b.u, b.x)

"""
    isperiodic(b1::InterfaceBoundary, b2::InterfaceBoundary)

Return `true` when interface boundaries `b1` and `b2` form matching periodic
boundary conditions.

# Arguments
- `b1`, `b2`: Interface boundaries to compare.

# Returns
- `Bool`: Whether the variable operations and interface coordinates are paired.
"""
function isperiodic(
        b1::InterfaceBoundary{b1u, b1u2},
        b2::InterfaceBoundary{b2u, b2u2}
    ) where {b1u, b1u2, b2u, b2u2}
    us_equal = isequal(operation(b1.u), operation(b2.u2)) &&
        isequal(operation(b2.u), operation(b1.u2))
    xs_equal = isequal(b1.x, b2.x2) && isequal(b1.x2, b2.x)
    return us_equal && xs_equal
end

"""
    isinterface(b)

Return `Val(true)` when `b` is an `InterfaceBoundary`, otherwise `Val(false)`.

# Arguments
- `b`: A boundary representation.

# Returns
- `Val{true}` for `InterfaceBoundary` and `Val{false}` otherwise.
"""
@inline function isinterface(b)
    if b isa InterfaceBoundary
        return Val(true)
    else
        return Val(false)
    end
end

"""
    filter_interfaces(bs)

Return the interface boundaries contained in the collection `bs`.

# Arguments
- `bs`: An iterable of boundary representations.

# Returns
- A vector containing only `InterfaceBoundary` values.
"""
filter_interfaces(bs) = filter(b -> b isa InterfaceBoundary, bs)

"""
    haslowerupper(bs, x)

Return booleans indicating whether interface boundaries in `bs` include lower
and upper boundaries for independent variable `x`.

# Arguments
- `bs`: An iterable of boundary representations.
- `x`: The independent variable to inspect.

# Returns
- `Tuple{Bool, Bool}`: `(has_lower, has_upper)`.
"""
function haslowerupper(bs, x)
    haslower = false
    hasupper = false
    for b in filter_interfaces(bs)
        if isequal(b.x, x)
            if isupper(b)
                hasupper = true
            else
                haslower = true
            end
        end
    end
    return haslower, hasupper
end

"""
    has_interfaces(bmps)

Return `true` when a nested boundary map contains at least one
`InterfaceBoundary`.

# Arguments
- `bmps`: A nested boundary map indexed by dependent variable and independent
  variable.

# Returns
- `Bool`: Whether at least one interface boundary is present.
"""
function has_interfaces(bmps)
    return any(b -> b isa InterfaceBoundary, reduce(vcat, reduce(vcat, collect.(values.(collect(values(bmps)))))))
end

# indexes for Iedge depending on boundary type
"""
    isupper(b)

Return `true` when boundary object `b` corresponds to the upper side of a
domain or interface.

# Arguments
- `b`: A boundary representation.

# Returns
- `Bool`: Whether `b` is oriented to the upper side.
"""
isupper(::LowerBoundary) = false
isupper(::UpperBoundary) = true
function isupper(::InterfaceBoundary{IsUpper_u}) where {IsUpper_u}
    return IsUpper_u isa Val{true} ? true : false
end
isupper(::HigherOrderInterfaceBoundary) = true

"""
    flatten_vardict(bmps)

Flatten a nested boundary map into a vector of boundary objects.

# Arguments
- `bmps`: A nested boundary map indexed by dependent variable and independent
  variable.

# Returns
- `Vector`: All boundary objects in the map, preserving map traversal order.

# Examples
```julia
boundaries = flatten_vardict(boundarymap)
```
"""
flatten_vardict(bmps) = reduce(vcat, reduce(vcat, collect.(values.(collect(values(bmps))))))

function _boundary_rules(v, orders, u, x, val)
    args = v.args[operation(u)]
    args = substitute.(args, (x => val,))
    varrule = operation(u)(args...) => [operation(u)(args...), x, 0]

    spacerules = [
        (Differential(x)^d)(operation(u)(args...)) => [
            operation(u)(args...), x, d,
        ] for d in reverse(orders[x])
    ]

    if v.time !== nothing && x !== v.time
        timerules = [
            Differential(v.time)(operation(u)(args...)) => [
                operation(u)(args...), x, d,
            ] for d in reverse(orders[v.time])
        ]
        return vcat(spacerules, timerules, varrule)
    else
        return vcat(spacerules, varrule)
    end
end

function generate_boundary_matching_rules(v, orders)
    # TODO: Check for bc equations of multiple variables
    lowerboundary(x) = v.intervals[x][1]
    upperboundary(x) = v.intervals[x][2]

    # Rules to match boundary conditions on the lower boundaries
    lower = Dict(
        [
            operation(u) => Dict(
                [
                    x => _boundary_rules(
                        v, orders, u, x, lowerboundary(x)
                    )
                        for x in all_ivs(u, v)
                ]
            ) for u in v.ū
        ]
    )

    upper = Dict(
        [
            operation(u) => Dict(
                [
                    x => _boundary_rules(
                        v, orders, u, x, upperboundary(x)
                    )
                        for x in all_ivs(u, v)
                ]
            ) for u in v.ū
        ]
    )

    return (lower, upper)
end
"""
    parse_bcs(bcs, v::VariableMap, orders, disc::AbstractDiscretization)
    parse_bcs(bcs, v::VariableMap, orders)

Create the boundary map of a `PDESystem`'s boundary and initial conditions, used later when
discretizing the boundary condition equations. Public but not exported, so call it as
`PDEBase.parse_bcs`. `symbolic_discretize` calls the four-argument form; the three-argument
form is kept for compatibility and equals the method for a plain
`AbstractEquationSystemDiscretization`.

The map is `boundarymap[operation(u)][key]`, a vector of boundary objects for each dependent
variable `u`. `orders[x]` lists the derivative orders in `x` that appear in `bcs`. The
four-argument form dispatches on the discretization family:

  - `AbstractEquationSystemDiscretization` (the default): matches each condition by
    substituting interval endpoints. The key is the independent variable whose endpoint the
    condition is on, and initial conditions are under the time variable.
  - `AbstractBoxDiscretization`: the same keys and boundary objects as the default, found by
    classifying each call of a dependent variable with [`classify_call`](@ref) and
    [`location`](@ref PDEBase.location). An endpoint within the tolerance of
    [`classify_argument`](@ref) counts as the endpoint, and a trace `u(t, Γ)` on a face of
    a box is the same face as its endpoint spelling.
  - `AbstractMeshDiscretization`: a [`PieceBoundary`](@ref) for each condition written
    with traces, under the key of its [`BoundaryPiece`](@ref); every piece of the spatial
    domain has a key. Initial conditions are `LowerBoundary`s under the time variable.
    Other conditions are an error.
"""
function parse_bcs(bcs, v::VariableMap, orders)
    lower_boundary_rules, upper_boundary_rules = generate_boundary_matching_rules(v, orders)
    rules = (lo = flatten_vardict(lower_boundary_rules), hi = flatten_vardict(upper_boundary_rules))
    ## BC matching rules, returns the variable and parameter the bc concerns
    matches(term, side) = [r.second for r in rules[side] if subsmatch(term, r)]
    return _parse_bcs(bcs, v, matches)
end

parse_bcs(bcs, v::VariableMap, orders, ::AbstractEquationSystemDiscretization) = parse_bcs(bcs, v, orders)

function parse_bcs(bcs, v::VariableMap, orders, ::AbstractBoxDiscretization)
    _check_traced_conditions(bcs, v)
    return _parse_bcs(map(bc -> _snap_endpoints(bc, v), bcs), v, _face_matcher(v, orders))
end

function parse_bcs(bcs, v::VariableMap, orders, disc::AbstractMeshDiscretization)
    t = v.time
    Ω = spatial_domain(v)
    keys_ = t === nothing ? Any[] : Any[t]
    Ω !== nothing && applicable(boundarypieces, Ω) && append!(keys_, boundarypieces(Ω))
    boundarymap = Dict([operation(u) => Dict{Any, Any}([k => [] for k in keys_]) for u in v.ū])
    matches = _face_matcher(v, orders)
    for bc in bcs
        class = classify_condition(v, bc)
        trace = _first_trace(bc)
        if class.kind === :initial
            _parse_bc!(boundarymap, _snap_endpoints(bc, v), v, matches)
        elseif class.kind === :boundary && trace !== nothing
            piece = location(classify_call(v, trace)[2]).piece
            push!(get!(() -> [], boundarymap[operation(trace)], piece), PieceBoundary(trace, piece, bc, v))
        else
            throw(
                ArgumentError(
                    "boundary condition $(bc) is not written with traces on boundary pieces, which " *
                        "$(nameof(typeof(disc))) requires; write u(t, Γ) for a piece Γ of the spatial domain"
                )
            )
        end
    end
    return boundarymap
end

# Shared by every `parse_bcs` method that files conditions under interval endpoints:
# `matches(term, side)` returns the `(u, x, order)` of every boundary rule that `term`
# matches on the `:lo` or `:hi` endpoints, in rule order, where `u` is the dependent
# variable with `x` replaced by the endpoint and `order` the derivative order of the rule.
function _parse_bcs(bcs, v::VariableMap, matches)
    boundarymap = Dict([operation(u) => Dict([x => [] for x in all_ivs(v)]) for u in v.ū])

    # Generate initial conditions and bc equations
    for bc in bcs
        _parse_bc!(boundarymap, bc, v, matches)
    end
    return boundarymap
end

function _parse_bc!(boundarymap, bc, v::VariableMap, matches)
    t = v.time
    boundary = nothing
    # Split out additive terms
    terms = split_terms(bc)
    # * Assume that the BC is defined on the edge of the domain
    # Check whether the bc is on the lower boundary, or periodic, we don't care which depvar/var
    u_ = u__ = x_ = x__ = nothing
    isinterface = false
    interface_orders = []
    for term in terms, m in matches(term, :lo)
        # Need to check the order, higher order should be discretized as a normal upper index
        # As zeroth order goes on the lower index. Needs a zeroth order interface to make sense,
        # but I assume that this will always be the case where there are higher order interfaces.
        # Get the matched variables from the rule
        u_, x_, order = m
        # Mark the boundary
        if boundary === nothing
            boundary = LowerBoundary(u_, t, x_, order, bc, v)
        end
        # do it again for the upper end to check for periodic, but only check the current depvar and indvar
        push!(interface_orders, order)
        for term_ in setdiff(terms, [term]), m_ in matches(term_, :hi)
            isinterface = true
            u__, x__, order__ = m_
            @assert ndims(u_, v) == ndims(u__, v) "Invalid Interface Boundary $bc: Dependent variables $(u_) and $(u__) have different numbers of dimensions."

            push!(interface_orders, order__)
        end
        # Handle flux condition at interface
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
        for term in terms, m in matches(term, :hi)
            u_, x_, order = m
            boundary = UpperBoundary(u_, t, x_, order, bc, v)
            break
        end
    end
    @assert boundary !== nothing "Boundary condition $bc is not on a boundary of the domain, or is not a valid boundary condition"
    if boundary isa Tuple
        push!(boundarymap[operation(boundary[1].u)][boundary[1].x], boundary[1])
        push!(boundarymap[operation(boundary[2].u)][boundary[2].x], boundary[2])
    else
        push!(boundarymap[operation(boundary.u)][boundary.x], boundary)
    end
    return boundarymap
end

# The classifier's counterpart of the rules of `generate_boundary_matching_rules`, with the
# same rule order: for each dependent variable `u` and each of its arguments `x`, a call
# of `u` on the face of `x` matches with order `d` when its derivative order in `x` (see
# `_face_order`) is `d ∈ reverse(orders[x])`, then with order `d ∈ reverse(orders[t])` when
# it is the argument of `Differential(t)` for a spatial `x`, and then with order 0 wherever
# it appears.
function _face_matcher(v::VariableMap, orders)
    # Constructed as the rule dictionaries are, so that it iterates in their order.
    keyorder = Dict([operation(u) => Dict([x => nothing for x in all_ivs(u, v)]) for u in v.ū])
    candidates = [(op, x) for (op, xs) in keyorder for x in keys(xs)]
    return function (term, side)
        calls = _calls_with_parents!(Tuple{Any, Any}[], term, (), v)
        found = []
        for (op, x) in candidates
            face = BoxFace(x, side)
            chains = [c for (c, call) in calls if isequal(operation(call), op) && isequal(_face_location(v, call), face)]
            isempty(chains) && continue
            u = _face_call(v, op, x, v.intervals[x][side === :lo ? 1 : 2])
            for d in reverse(orders[x])
                any(c -> _face_order(v, c, x, face) == d, chains) && push!(found, (u, x, d))
            end
            if v.time !== nothing && x !== v.time
                for d in reverse(orders[v.time])
                    any(c -> !isempty(c) && isequal(first(c), Differential(v.time)), chains) && push!(found, (u, x, d))
                end
            end
            push!(found, (u, x, 0))
        end
        return found
    end
end

_face_call(v, op, x, val) = op(substitute.(v.args[op], (x => val,))...)

# `bc` with every call on a box face that spells its endpoint only to within the tolerance of
# `classify_argument`, as `u(t, 0.1 + 0.2)` on `Interval(0.0, 0.3)`, rewritten to the exact
# endpoint, so that the boundary's `eq` and `indvars` agree with its `u`. A condition with
# no such call is returned as is.
function _snap_endpoints(bc, v::VariableMap)
    lhs, rhs = _get_lhs_rhs(bc)
    rules = Dict()
    for (_, call) in _calls_with_parents!(_calls_with_parents!(Tuple{Any, Any}[], lhs, (), v), rhs, (), v)
        any(_ispiece, arguments(call)) && continue
        face = _face_location(v, call)
        face isa BoxFace || continue
        exact = _face_call(v, operation(call), face.iv, v.intervals[face.iv][face.side === :lo ? 1 : 2])
        isequal(call, exact) || (rules[call] = exact)
    end
    isempty(rules) && return bc
    lhs, rhs = substitute_traces(lhs, rules), substitute_traces(rhs, rules)
    return bc isa Pair ? lhs => rhs : lhs ~ rhs
end

_is_face_normal(v, op, face) = op isa NormalDerivative && isequal(location(v, [Piece(op.piece)]), face)

# The derivative order in `x` of a call on `face`, given the operations enclosing it,
# innermost first. Without a normal derivative it is the order of the call's immediate
# `Differential(x, d)`, as the substitution rules match `Differential(x)^d` applied to the
# call, and 0 otherwise. A run of derivatives in `x` and normal derivatives on `face`
# directly enclosing the call adds up, each normal derivative counting 1, as nested
# `Differential(x)`s fold into one `Differential(x, d)` in the endpoint spelling.
function _face_order(v, chain, x, face)
    isempty(chain) && return 0
    inx(op) = op isa Differential && isequal(op.x, x)
    total, normal = 0, false
    for op in chain
        if inx(op)
            total += op.order
        elseif _is_face_normal(v, op, face)
            total += 1
            normal = true
        else
            break
        end
    end
    normal && return total
    return inx(first(chain)) ? first(chain).order : 0
end

# Every call of a dependent variable in `ex`, with the operations enclosing it, innermost
# first.
function _calls_with_parents!(found, ex, enclosing, v)
    ex = safe_unwrap(ex)
    iscall(ex) || return found
    op = operation(ex)
    haskey(v.args, op) && push!(found, (enclosing, ex))
    for arg in arguments(ex)
        _calls_with_parents!(found, arg, (op, enclosing...), v)
    end
    return found
end

# The location of a call, or `nothing` when it does not classify: such a call is not on a
# face, as it matches no rule of the substitution matcher.
function _face_location(v, call)
    return try
        location(v, classify_call(v, call)[2])
    catch e
        e isa Union{ArgumentError, InexactError} || rethrow()
        nothing
    end
end

# The first trace `u(..., Γ, ...)` in `ex`, or `nothing`.
function _first_trace(ex)
    ex = safe_unwrap(ex)
    iscall(ex) || return nothing
    any(_ispiece, arguments(ex)) && return ex
    for arg in arguments(ex)
        found = _first_trace(arg)
        found === nothing || return found
    end
    return nothing
end
function _first_trace(eq::Union{Equation, Pair})
    lhs, rhs = _get_lhs_rhs(eq)
    found = _first_trace(lhs)
    return found === nothing ? _first_trace(rhs) : found
end

# Conditions written with traces must classify, so that a normal derivative on another
# piece than its trace is an error rather than an unmatched condition.
function _check_traced_conditions(bcs, v)
    for bc in bcs
        any(ex -> _find_trace(ex) !== nothing, _get_lhs_rhs(bc)) && classify_condition(v, bc)
    end
    return nothing
end
