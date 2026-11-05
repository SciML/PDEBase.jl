# Strong-Form Analysis

`PDEBase.StrongForm` analyses the strong form of a `PDESystem` for discretizers that need the
flux structure of each equation, such as finite-element, collocation or discontinuous Galerkin
methods. For each equation `e = lhs - rhs` it finds the mass part, the flux `G` and the
remainder `R` with `e = m⋅u̇ + Σᵢ ∂ᵢGᵢ + R`, all written in symbolic quadrature-point slots.

```@docs
PDEBase.StrongForm
PDEBase.StrongForm.analyse
PDEBase.StrongForm.StrongFormAnalysis
PDEBase.StrongForm.StateSlots
PDEBase.StrongForm.SlotContext
PDEBase.StrongForm.declared_fields
PDEBase.StrongForm.IBPPolicy
PDEBase.StrongForm.ProductRule
PDEBase.StrongForm.Weighted
PDEBase.StrongForm.Reject
PDEBase.StrongForm.EssentialBC
PDEBase.StrongForm.NaturalBC
PDEBase.StrongForm.PeriodicBC
PDEBase.StrongForm.check_whitelist
PDEBase.StrongForm.split_time
PDEBase.StrongForm.distribute_derivatives
PDEBase.StrongForm.split_divergence
PDEBase.StrongForm.NONSMOOTH_FUNCTIONS
PDEBase.StrongForm.DISCONTINUOUS_FUNCTIONS
PDEBase.StrongForm.slot_volume
PDEBase.StrongForm.check_evaluable
```

## Boundary conditions

`analyse` classifies each spatial boundary condition of the `parse_bcs` boundary map as
essential (`u_f = g`), natural (the facet flux `Φ` of the one equation whose conormal flux the
condition is proven proportional to) or periodic. A box discretization's conditions are keyed by
their `PDEBase.BoxFace`, and the `PieceBoundary` entries of a mesh discretization's map by their
`BoundaryPiece`. Traces are located with the classifier and put into slots with
`PDEBase.substitute_traces`, never derivative-expanded first: `NormalDerivative(Γ)(u_f(t, Γ))`
becomes `N·G_f`, with the numeric normal `±e_iv` on a face of the box and the symbolic normal
slots `N` on any other piece, whose normal the backend supplies per quadrature point.

A condition may also write the normal itself with `n = PDEBase.normal(Γ)`: component `n[d]`
(`d` is the position in `indvars(v)`) becomes `N[d]` on a piece and the numeric `±1` or `0` on a
face of the box. This expresses terms that `NormalDerivative` cannot, such as `Σ_d n_d ∂_c E_d`
of a vector-valued field or the tangential impedance term `n × (n × E)`. A Cartesian derivative
such as `Dx(u(t, Γ))` is accepted on a face of the box, and on any other piece only in a
condition that also contains `normal(Γ)` of that piece, where it becomes the gradient slot
`G`; it is then the interior gradient evaluated on the facet, not a surface gradient. Without
the normal it is an error that suggests `NormalDerivative`. A facet normal is also an error in
a volume equation, under a derivative, on another piece than the condition's, and in a
condition with no derivative of a trace (such a condition would be essential, and an essential
condition cannot depend on the normal). The conormal matching needs no change: in a condition
such as the absorbing-boundary term `n × curl E + γ n × (n × E) ~ g`, the derivative terms
are matched against the conormal flux of one equation, and the normal-dependent impedance term
lands in the facet flux `Φ`, a product of `N` and `U`. Matching the curl-type derivative terms
of a two-dimensional in-plane vector problem is ambiguous, since both equations' conormal
coefficient vectors are parallel.
`pair_equations` is the helper a discretizer uses to give each equation its test field.

```@docs
PDEBase.StrongForm.boundary_conditions
PDEBase.StrongForm.slot_condition
PDEBase.StrongForm.restrict_to_face
PDEBase.StrongForm.match_flux
PDEBase.StrongForm.route_natural
PDEBase.StrongForm.essential_condition
PDEBase.StrongForm.periodic_condition
PDEBase.StrongForm.pair_equations
```
