
"""
A function that creates a tuple of CartesianIndices of unit length and `N` dimensions, one pointing along each dimension.
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
Get a unit `CartesianIndex` in dimension `j` of length `N`.
"""
unitindex(N, j) = CartesianIndex(ntuple(i -> i == j, N))

remove(args, t) = filter(x -> t === nothing || !isequal(safe_unwrap(x), safe_unwrap(t)), args)
remove(v::AbstractVector, a::Number) = filter(x -> !isequal(x, a), v)



d_orders(x, pdeeqs) = reverse(sort(collect(union((differential_order(pde.rhs, safe_unwrap(x)) for pde in pdeeqs)..., (differential_order(pde.lhs, safe_unwrap(x)) for pde in pdeeqs)...))))

insert(args...) = insert!(args[1], args[2:end]...)

@inline function sym_dot(a, b)
    mapreduce((+), zip(a, b)) do (a_, b_)
        a_ * b_
    end
end

vcat!(a::AbstractArray, b::AbstractArray) = append!(a, b)
vcat!(a::AbstractArray, b::AbstractArray...) = append!(a, vcat(b...))
vcat!(a::AbstractArray, b...) = append!(a, vcat(b...))
vcat!(a::AbstractArray, b) = push!(a, b)
