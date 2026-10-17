"""
    unitindices(N::Int)

Create one unit `CartesianIndex` for each of the `N` dimensions.

# Arguments
- `N::Int`: The number of dimensions.

# Returns
- A vector of `CartesianIndex` values, or `CartesianIndex()` when `N == 0`.

# Examples
```julia
unitindices(2) == [CartesianIndex(1, 0), CartesianIndex(0, 1)]
```
"""
function unitindices(N::Int) #create unit CartesianIndex for each dimension
    null = zeros(Int, N)
    if N == 0
        return CartesianIndex()
    else
        return map(1:N) do i
            unit_i = copy(null)
            unit_i[i] = 1
            CartesianIndex(Tuple(unit_i))
        end
    end
end

"""
    unitindex(N, j)

Create a unit `CartesianIndex` in dimension `j` of an `N`-dimensional index.

# Arguments
- `N`: The number of dimensions.
- `j`: The dimension containing the unit entry.

# Returns
- `CartesianIndex`: An index with one in dimension `j` and zero elsewhere.
"""
unitindex(N, j) = CartesianIndex(ntuple(i -> i == j, N))

function remove(args, t)
    return filter(x -> t === nothing || !isequal(safe_unwrap(x), safe_unwrap(t)), args)
end
remove(v::AbstractVector, a::Number) = filter(x -> !isequal(x, a), v)

"""
    d_orders(x, pdeeqs)

Return derivative orders with respect to independent variable `x` that appear
in the PDE equations or boundary conditions `pdeeqs`.

# Arguments
- `x`: The independent variable to inspect.
- `pdeeqs`: An iterable of `Equation` or `Pair` values.

# Returns
- A descending vector of unique derivative orders.

# Examples
```julia
d_orders(x, [Dx(Dx(u)) ~ 0]) == [2]
```
"""
function d_orders(x, pdeeqs)
    # Handle both Equation (has lhs, rhs) and Pair (has first, second) types
    _get_rhs(pde) = safe_unwrap(pde isa Pair ? pde.second : pde.rhs)
    _get_lhs(pde) = safe_unwrap(pde isa Pair ? pde.first : pde.lhs)
    orders = Set{Int}()
    for pde in pdeeqs
        union!(orders, differential_order(_get_rhs(pde), safe_unwrap(x)))
        union!(orders, differential_order(_get_lhs(pde), safe_unwrap(x)))
    end
    return sort!(collect(orders); rev = true)
end

insert(args...) = insert!(args[1], args[2:end]...)

@inline function sym_dot(a, b)
    return mapreduce((+), zip(a, b)) do (a_, b_)
        a_ * b_
    end
end

show_verbosemode(x) = false

vcat!(a::AbstractArray, b::AbstractArray) = append!(a, b)
vcat!(a::AbstractArray, b...) = append!(a, vcat(b...))
vcat!(a::AbstractArray, b) = push!(a, b)
