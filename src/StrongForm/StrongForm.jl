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
using Symbolics: Symbolics, Differential, Equation, Num, substitute, unwrap, wrap
using SymbolicUtils: SymbolicUtils, arguments, iscall, operation, symtype, unwrap_const
using TermInterface: maketerm, metadata
using ModelingToolkitBase: get_eqs
using ..PDEBase: VariableMap, LowerBoundary, UpperBoundary, InterfaceBoundary, _flatten_eqs,
    chain_flatten_array_variables, depvars, indvars, isupper, replaced_vars

include("types.jl")
include("whitelist.jl")
include("split_time.jl")
include("split_divergence.jl")
include("slots.jl")
include("boundary.jl")
include("pairing.jl")
include("analyse.jl")

export analyse

@public StateSlots, IBPPolicy, ProductRule, Weighted, Reject, EssentialBC, NaturalBC,
    PeriodicBC, StrongFormAnalysis, check_whitelist, split_time, split_divergence,
    distribute_derivatives, NONSMOOTH_FUNCTIONS, DISCONTINUOUS_FUNCTIONS, SlotContext, slot_volume, check_evaluable

end # module StrongForm
