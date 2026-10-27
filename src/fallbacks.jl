### UPSTREAM SYNTAX FALLBACKS ###

# Each definition below is used only while the installed upstream package lacks it, so
# that it disappears once Symbolics or DomainSets releases its own. Keep the semantics
# identical to the upstream proposals: JuliaSymbolics/Symbolics.jl#2115 for
# `NormalDerivative`, and the DomainSets `BoundaryPiece` proposal.

@static if isdefined(Symbolics, :NormalDerivative)
    using Symbolics: NormalDerivative
else
    @doc """
        NormalDerivative(Γ)

    Represents the normal derivative operator on a boundary piece `Γ`. `Dn =
    NormalDerivative(Γ)` applied to a trace `u(t, Γ)` (a call of a dependent variable whose
    spatial arguments are replaced by the domain `Γ`) gives `Dn(u(t, Γ))`, printed as
    `Dn_Γ(u(t, Γ))`: the trace of the normal derivative. That is, differentiate `u` along
    the outward normal of `Γ`, then evaluate on `Γ`.

    Like `Differential`, applying a `NormalDerivative` only builds a symbolic term, with
    the operator as its operation and the trace as its single argument. `expand_derivatives`
    leaves the term alone. This type does not check that the piece of the operator matches
    the piece in the trace; that is up to the consumer. Both are recoverable from the term:
    the operator's piece is `operation(term).piece`, and the trace's is among
    `arguments(only(arguments(term)))`.

    This is PDEBase's fallback for the definition proposed in
    JuliaSymbolics/Symbolics.jl#2115, used only while the installed Symbolics lacks
    `NormalDerivative`; once Symbolics defines it, `PDEBase.NormalDerivative` is
    `Symbolics.NormalDerivative`. Two behaviours of the upstream version need changes inside
    Symbolics and are absent from the fallback:

      - `expand_derivatives(Dt(Dn(u(t, Γ))))` returns `0` rather than staying unevaluated,
        because Symbolics treats every unknown operator application as a fresh variable.
        Never derivative-expand an equation holding traces.
      - `Symbolics.substitute_in_deriv` does not descend into `Dn(...)`. Replace traces with
        PDEBase's own trace substitution instead.

    # Fields
    - `piece`: The boundary piece whose outward normal the derivative is taken along.
    """
    struct NormalDerivative{D <: DomainSets.Domain} <: SymbolicUtils.Operator
        piece::D
    end
    # Upstream defines this for `BasicSymbolic{VartypeT}` only; `VartypeT` is not public,
    # so the vartype is taken from the argument, which is identical for Symbolics' terms.
    function (Dn::NormalDerivative)(x::SymbolicUtils.BasicSymbolic{T}) where {T}
        return SymbolicUtils.BSImpl.Term{T}(
            Dn, SymbolicUtils.ArgsT{T}((x,));
            type = SymbolicUtils.symtype(x), shape = SymbolicUtils.shape(x)
        )
    end
    (Dn::NormalDerivative)(x::Num) = Num(Dn(unwrap(x)))
    (Dn::NormalDerivative)(x::Symbolics.Arr{T, N}) where {T, N} = Symbolics.Arr{T, N}(Dn(unwrap(x)))
    SymbolicUtils.isbinop(::NormalDerivative) = false

    Base.:(==)(a::NormalDerivative, b::NormalDerivative) = isequal(a, b)
    Base.isequal(a::NormalDerivative, b::NormalDerivative) = isequal(a.piece, b.piece)
    Base.hash(a::NormalDerivative, h::UInt) = hash(a.piece, hash(NormalDerivative, h))

    function Base.show(io::IO, Dn::NormalDerivative)
        label = hasproperty(Dn.piece, :name) ? Dn.piece.name : Dn.piece
        return print(io, "Dn_", label)
    end
end

@static if isdefined(DomainSets, :BoundaryPiece)
    using DomainSets: BoundaryPiece
else
    @doc """
        BoundaryPiece(parent, name::Symbol, geometry = nothing)

    A named piece of the boundary of the domain `parent`. `geometry` is the piece as an
    analytic sub-domain, or `nothing` when the piece is known only by name, as for a facet
    set of a mesh. A piece is used in place of the spatial arguments of a dependent
    variable to write a trace, `u(t, Γ)`, and as the argument of a [`NormalDerivative`](@ref).
    It prints as `Γ_<name>`, and `parent(Γ)` returns the parent domain.

    Pieces compare and hash structurally, by parent, name and geometry. `isless` orders
    pieces by name, then by the hashes of parent and geometry, so that Symbolics can order
    the terms of a sum of traces such as `u(t, Γa) - u(t, Γb)`.

    This is PDEBase's fallback for the type proposed for DomainSets, used only while the
    installed DomainSets lacks `BoundaryPiece`.

    # Fields
    - `parent`: The domain whose boundary this is a piece of.
    - `name`: The name of the piece.
    - `geometry`: The piece as an analytic domain, or `nothing`.
    """
    struct BoundaryPiece{T, P, G} <: DomainSets.Domain{T}
        parent::P
        name::Symbol
        geometry::G
    end
    function BoundaryPiece(parent::DomainSets.Domain{T}, name::Symbol, geometry = nothing) where {T}
        return BoundaryPiece{T, typeof(parent), typeof(geometry)}(parent, name, geometry)
    end

    Base.parent(p::BoundaryPiece) = p.parent
    Base.show(io::IO, p::BoundaryPiece) = print(io, "Γ_", p.name)

    Base.:(==)(a::BoundaryPiece, b::BoundaryPiece) = isequal(a, b)
    function Base.isequal(a::BoundaryPiece, b::BoundaryPiece)
        return a.name === b.name && isequal(a.parent, b.parent) &&
            isequal(a.geometry, b.geometry)
    end
    function Base.hash(p::BoundaryPiece, h::UInt)
        return hash(p.geometry, hash(p.parent, hash(p.name, hash(BoundaryPiece, h))))
    end
    # SymbolicUtils orders sum terms by comparing same-typed constant arguments with `<`.
    # Lexicographic on a key, so a strict weak order; consistent with `isequal` because
    # `hash` is (IntervalSets' empty intervals are the known exception: `isequal` but
    # unequal hashes). Stable across processes except for parents hashed by identity.
    _order_key(p::BoundaryPiece) = (p.name, hash(p.parent), hash(p.geometry))
    Base.isless(a::BoundaryPiece, b::BoundaryPiece) = isless(_order_key(a), _order_key(b))
end

# Face names of a box, in the order of `components(boundary(d))`: a rectangle's
# bottom, right, top and left faces, and in other dimensions the lower then upper face of
# each axis, `:x1_lo, :x1_hi, :x2_lo, ...`.
const _RECTANGLE_FACES = (:bottom, :right, :top, :left)
const _RECTANGLE_FACE_AXES = (bottom = (2, :lo), right = (1, :hi), top = (2, :hi), left = (1, :lo))

function _box_face_names(n::Int)
    n == 2 && return collect(_RECTANGLE_FACES)
    return [Symbol("x", axis, "_", side) for axis in 1:n for side in (:lo, :hi)]
end

# `(axis, side)` of a face of a box, with `side` either `:lo` or `:hi`.
function _box_face(p::BoundaryPiece)
    haskey(_RECTANGLE_FACE_AXES, p.name) && return _RECTANGLE_FACE_AXES[p.name]
    m = match(r"^x(\d+)_(lo|hi)$", String(p.name))
    m === nothing && throw(ArgumentError("$(p) is not a face of a box"))
    return (parse(Int, m[1]), Symbol(m[2]))
end

@static if isdefined(DomainSets, :boundarypieces)
    using DomainSets: boundarypieces
else
    @doc """
        boundarypieces(d)

    Return the named [`BoundaryPiece`](@ref)s that make up the boundary of `d`.

    For a box (`DomainSets.HyperRectangle`) the pieces are its faces, in the order of
    `components(boundary(d))`, with each face as the piece's geometry. A rectangle's faces
    are named `:bottom, :right, :top, :left`. In any other dimension `n` they are named by
    axis index, lower then upper per axis: `:x1_lo, :x1_hi, :x2_lo, …, :xn_hi`. An interval
    is the one-dimensional box: its pieces are its endpoints `:x1_lo` and `:x1_hi`, with
    the endpoint as a `DomainSets.Point`.
    """
    function boundarypieces end

    function boundarypieces(d::Union{DomainSets.HyperRectangle, DomainSets.AbstractInterval})
        faces = DomainSets.components(DomainSets.boundary(d))
        names = _box_face_names(DomainSets.dimension(d))
        return [BoundaryPiece(d, name, face) for (name, face) in zip(names, faces)]
    end
end

@static if isdefined(DomainSets, :boundarypiece)
    using DomainSets: boundarypiece
else
    @doc """
        boundarypiece(d, name::Symbol)

    Return the [`BoundaryPiece`](@ref) of the boundary of `d` named `name`, one of the
    names of [`boundarypieces(d)`](@ref boundarypieces). An unknown name is an error that
    lists the valid names.
    """
    function boundarypiece(d, name::Symbol)
        pieces = boundarypieces(d)
        for p in pieces
            p.name === name && return p
        end
        valid = join((repr(p.name) for p in pieces), ", ")
        throw(ArgumentError("$(d) has no boundary piece named $(repr(name)); valid names are $(valid)"))
    end
end
