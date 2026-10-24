"""
    PDEBase.StrongForm

Symbolic analysis of a PDE system's strong form, shared by discretizers that need
the flux structure of each equation (finite elements, collocation, discontinuous
Galerkin). Each equation `e = lhs - rhs` is checked against an operator whitelist
and split into a time-derivative part and a spatial part, ready to be written as
`e = m⋅u̇ + Σᵢ ∂ᵢGᵢ + R`.
"""
module StrongForm

import SciMLPublic: @public
using Symbolics: Symbolics, Differential, Equation, Num, unwrap
using SymbolicUtils: SymbolicUtils, arguments, iscall, operation, symtype

include("types.jl")
include("whitelist.jl")
include("split_time.jl")

@public StateSlots, IBPPolicy, ProductRule, Weighted, Reject, EssentialBC, NaturalBC,
    PeriodicBC, StrongFormAnalysis, check_whitelist, split_time

end # module StrongForm
