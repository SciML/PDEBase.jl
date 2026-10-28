# [Argument Classification](@id classify)

PDEBase recognises where a call of a dependent variable lies by classifying each of its
arguments. `u(t, x, y)`, the endpoint form `u(t, 0, y)` and the trace `u(t, Γ)` go
through the same code: an argument is the declared variable (`Free`), an interval
endpoint (`BoxFace`), or a boundary piece standing for all the spatial arguments
(`Piece`). A piece that is a face of the box normalises to the `BoxFace` of its endpoint
form, so on `(x, y) ∈ Rectangle(0 .. 1, 0 .. 1)` the calls `u(t, 0, y)` and
`u(t, boundarypiece(Ω, :left))` have the same location.

```julia
v = VariableMap(pdesys, disc)
u_op, classes = classify_call(v, u(t, 0, y))  # [Free(t), BoxFace(x, :lo), Free(y)]
PDEBase.location(v, classes)                  # BoxFace(x, :lo)
classify_condition(v, u(t, 1, y) ~ u(t, 0, y)).kind  # :periodic
```

A trace keeps the time argument: on a time-dependent `u(t, x, y)` it is `u(t, Γ)`, and
`u(Γ)` is an error.

## Substitution through operators

Plain `substitute` does not descend into `Differential` or `NormalDerivative`
applications, so it leaves the trace in `Dx(u(t, Γ))` and `Dn(u(t, Γ))` unchanged.
Rewrite traces with `substitute_traces`, then check with `assert_no_traces` before any
`expand_derivatives`, which turns a remaining `Dx(u(t, Γ))` into `0`.

## Box lowering

`to_box_bcs` turns every trace condition on the faces of an interval or box into the
endpoint form box discretizations parse: `u(t, Γ_left) ~ g` becomes `u(t, 0.0, y) ~ g`,
`Dn_Γ(u(t, Γ))` becomes `-Dx(...)` on a lower face and `+Dx(...)` on an upper face, and
conditions already in endpoint form pass through. `symbolic_discretize` applies it
automatically for every `AbstractBoxDiscretization` through `lower_boundaries`. An interval
is the one-dimensional box, with pieces `:x1_lo` and `:x1_hi`.

```julia
Ω = Rectangle(0.0 .. 1.0, 0.0 .. 1.0)
Γ = boundarypiece(Ω, :right)
bcs = [u(0, x, y) ~ 0, u(t, 0, y) ~ sin(t), NormalDerivative(Γ)(u(t, Γ)) ~ 0]
@named pdesys = PDESystem(eqs, bcs, [(x, y) ∈ Ω, t ∈ Interval(0.0, 1.0)], [t, x, y], [u(t, x, y)])
get_bcs(to_box_bcs(pdesys))  # [u(0, x, y) ~ 0, u(t, 0, y) ~ sin(t), Dx(u(t, 1.0, y)) ~ 0]
```

```@docs
to_box_bcs
```

## Classifier

```@docs
classify_argument
classify_call
classify_condition
PDEBase.location
PDEBase.boxface
PDEBase.normal
substitute_traces
assert_no_traces
```

## Classes

```@docs
PDEBase.ArgumentClass
PDEBase.Free
PDEBase.BoxFace
PDEBase.Piece
PDEBase.Interior
PDEBase.ConditionClass
```
