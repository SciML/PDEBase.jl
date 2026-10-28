# Strong-Form Analysis

`PDEBase.StrongForm` analyses the strong form of a `PDESystem` for discretizers that need the
flux structure of each equation, such as finite-element, collocation or discontinuous Galerkin
methods. For each equation `e = lhs - rhs` it finds the mass part, the flux `G` and the
remainder `R` with `e = m⋅u̇ + Σᵢ ∂ᵢGᵢ + R`, all written in symbolic quadrature-point slots.

```@docs
PDEBase.StrongForm
PDEBase.StrongForm.StrongFormAnalysis
PDEBase.StrongForm.StateSlots
PDEBase.StrongForm.IBPPolicy
PDEBase.StrongForm.ProductRule
PDEBase.StrongForm.Weighted
PDEBase.StrongForm.Reject
PDEBase.StrongForm.EssentialBC
PDEBase.StrongForm.NaturalBC
PDEBase.StrongForm.PeriodicBC
PDEBase.StrongForm.check_whitelist
PDEBase.StrongForm.split_time
```
