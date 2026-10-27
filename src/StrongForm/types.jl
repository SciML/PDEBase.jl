"""
    StateSlots(nfields, dim, nparams)

Symbolic quadrature-point slots that every strong-form output is written in.

# Fields
- `U`: state, one slot per field.
- `G`: gradient, field-major: `G[(f - 1) * dim + d]` is `∂u_f/∂x_d`.
- `Ud`: time derivative `∂ₜu_f`, one slot per field.
- `H`: Hessian, `H[((f - 1) * dim + d - 1) * dim + e]` is `∂²u_f/∂x_d∂x_e`.
- `X`: spatial coordinates, `T`: time, `N`: outward unit normal, `P`: parameters.
"""
struct StateSlots
    U::Vector{Num}
    G::Vector{Num}
    Ud::Vector{Num}
    H::Vector{Num}
    X::Vector{Num}
    T::Num
    N::Vector{Num}
    P::Vector{Num}
end

function StateSlots(nfields, dim, nparams)
    slots(name, n) = [Symbolics.variable(name, i) for i in 1:n]
    return StateSlots(
        slots(:U, nfields), slots(:G, nfields * dim), slots(:Ud, nfields),
        slots(:H, nfields * dim * dim), slots(:X, dim), Symbolics.variable(:T),
        slots(:N, dim), slots(:P, nparams)
    )
end

"""
    IBPPolicy

How a non-divergence second-order term `c⋅∂ᵢw` is integrated by parts: [`ProductRule`](@ref),
[`Weighted`](@ref) or [`Reject`](@ref).
"""
abstract type IBPPolicy end

"""
    ProductRule()

Default [`IBPPolicy`](@ref): rewrite `c⋅∂ᵢw = ∂ᵢ(c⋅w) - (∂ᵢc)⋅w`, adding `c⋅w` to the flux
and `-(∂ᵢc)⋅w` to the remainder. The correction term makes the stiffness matrix non-symmetric.
"""
struct ProductRule <: IBPPolicy end

"""
    Weighted(w)

[`IBPPolicy`](@ref) that multiplies the equation by the weight `w(x)` before integrating by
parts, e.g. `w = r²` for `1/r² ∂ᵣ(r² ∂ᵣu)`.
"""
struct Weighted{W} <: IBPPolicy
    w::W
end

"""
    Reject()

[`IBPPolicy`](@ref) that errors on any term that would need the product rule.
"""
struct Reject <: IBPPolicy end

"""
    EssentialBC(piece, field, g)

Essential condition `u_field = g(X, T, P)` on the boundary key `piece`; `g` is free of state.
`field` is 1-based.
"""
struct EssentialBC{K}
    piece::K
    field::Int
    g::Num
end

"""
    NaturalBC(piece, eq, Φ)

Matched boundary flux `Φ` of equation `eq` (1-based) on the boundary key `piece`, in slot form.
It contributes `∮ Φ v` to the weak form of `eq`.
"""
struct NaturalBC{K}
    piece::K
    eq::Int
    Φ::Num
end

"""
    PeriodicBC(master, mirror, field, shift)

Periodic identification of `field` (1-based) on the boundary key `mirror` with `master`,
related by the translation `shift`.
"""
struct PeriodicBC{K, S}
    master::K
    mirror::K
    field::Int
    shift::S
end

"""
    StrongFormAnalysis

The strong form of a PDE system in slot form ([`StateSlots`](@ref)). For equation `e`,

    strong[e] = mass[e] + Σᵢ ∂ᵢ flux[e][i] + remainder[e]

where `mass[e]` is the whole time-derivative part `Σ_f m_ef⋅Ud_f`, not a bare coefficient,
so cross-field mass is representable. `conormal[e]` is `∂(flux[e]⋅N)/∂G`, field-major like
`G`. Field `f` is the `f`-th declared dependent variable of the system (see
[`declared_fields`](@ref)), whatever order the equations mention them in. Boundary data is
keyed by the boundary key type `K`. `essential`, `natural` and `periodic` list their
conditions in the order of the system's boundary conditions; the analysis ranks nothing, so
which essential condition wins where faces meet is the discretizer's choice.

`mass`, `flux`, `remainder` and `conormal` never hold an unevaluated derivative. `strong` may:
the strong residual of `∂ₓ(a(x)⋅∂ₓu)` with `a` registered without a derivative rule holds
`Differential(X₁)(a(X₁))`. Consumers must call [`check_evaluable`](@ref) on it before
`build_function`, which would compile such a derivative to zero.
"""
struct StrongFormAnalysis{K}
    slots::StateSlots
    nfields::Int
    dim::Int
    strong::Vector{Num}
    mass::Vector{Num}
    flux::Vector{Vector{Num}}
    remainder::Vector{Num}
    conormal::Vector{Vector{Num}}
    essential::Vector{EssentialBC{K}}
    natural::Vector{NaturalBC{K}}
    periodic::Vector{PeriodicBC{K}}
end

# Change the boundary-key type, e.g. so a `Vector{PeriodicBC{Any}}` accepts a
# `PeriodicBC{BoxFace}`.
Base.convert(::Type{EssentialBC{K}}, b::EssentialBC) where {K} = EssentialBC{K}(b.piece, b.field, b.g)
Base.convert(::Type{NaturalBC{K}}, b::NaturalBC) where {K} = NaturalBC{K}(b.piece, b.eq, b.Φ)
Base.convert(::Type{PeriodicBC{K}}, b::PeriodicBC) where {K} =
    PeriodicBC{K, typeof(b.shift)}(b.master, b.mirror, b.field, b.shift)
