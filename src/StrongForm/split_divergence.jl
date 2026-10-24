"""
    NONSMOOTH_FUNCTIONS

Functions that are not differentiable in the symbolic sense. Their Symbolics derivatives are 0 or
`ifelse(...)`, so a product-rule correction built from them would vanish silently (spec §5.1, M1).
Comparisons and Bool operators are included; they only occur inside `ifelse` conditions. A
derivative that comes out piecewise or `NaN` is rejected as well, catching functions not listed.
"""
const NONSMOOTH_FUNCTIONS = (
    ifelse, abs, max, min, sign, floor, ceil, round, trunc, signbit, clamp,
    mod, rem, mod1, rem2pi, div, fld, cld,
    <, <=, >, >=, ==, !=, &, |, !, xor,
)

const ORDINALS = ("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth")
ordinal(n) = n <= length(ORDINALS) ? ORDINALS[n] : "$(n)th"

depends_on(ex, x) = isequal(ex, x) || iscall(ex) && any(a -> depends_on(a, x), arguments(ex))

is_spatial_derivative(ex, sivs) =
    iscall(ex) && operation(ex) isa Differential && any(isequal(operation(ex).x), sivs)

"Total spatial order of the deepest derivative chain in `ex`, read from `Differential.order`."
function spatial_order(ex, sivs)
    iscall(ex) || return 0
    own = is_spatial_derivative(ex, sivs) ? operation(ex).order : 0
    return own + maximum((spatial_order(a, sivs) for a in arguments(ex)); init = 0)
end

"""
    distribute_derivatives(ex)

Rewrite every `D(a + b)` in `ex` as `D(a) + D(b)`, expanding derivative arguments first, so that
each additive term of the result is a product of factors. Never derivative-expands, so it is
safe on expressions holding traces.
"""
function distribute_derivatives(ex)
    ex = unwrap(ex)
    iscall(ex) || return ex
    op = operation(ex)
    args = map(distribute_derivatives, arguments(ex))
    op isa Differential || return maketerm(typeof(ex), op, args, metadata(ex))
    terms = filter(t -> !(unwrap_const(t) isa Number), addterms(unwrap(Symbolics.expand(only(args)))))
    return unwrap(sum(op, terms; init = Num(0)))
end

"Is `c`, as a function of `x`, built only from functions with symbolic derivatives?"
function is_smooth(c, x)
    depends_on(c, x) || return true
    iscall(c) || return true
    op, args = operation(c), arguments(c)
    any(f -> op === f, NONSMOOTH_FUNCTIONS) && return false
    op isa Function && op !== getindex &&
        !has_derivative_rules(op, args, a -> depends_on(a, x)) && return false
    return all(a -> is_smooth(a, x), args)
end

"Does `ex` hold an `ifelse` depending on `x`, or a `NaN`, as derivatives of jumps do?"
function is_piecewise(ex, x)
    iscall(ex) || return (v = unwrap_const(ex); v isa AbstractFloat && isnan(v))
    operation(ex) === ifelse && depends_on(ex, x) && return true
    return any(a -> is_piecewise(a, x), arguments(ex))
end

not_in_divergence_form(term) =
    error("second derivative not in a divergence-compatible position: $term")

"""
    split_divergence(spatial, sivs, policy::IBPPolicy = ProductRule()) -> (G, R, corrected)

Stage b of the strong-form analysis. Writes the sum of the spatial terms from [`split_time`](@ref)
as `Σᵢ ∂ᵢG[i] + sum(R)` over the spatial variables `sivs`. A term without a second-order factor
goes to `R` unchanged. A term `c⋅Dᵢ(w′)` with one second-order factor adds `c⋅w′` to `G[i]` and,
when `c` depends on `xᵢ`, the product-rule correction `-(Dᵢc)⋅w′` to `R`, unless the corrections
of all terms simplify to zero together; `corrected` reports whether corrections were added.
Under `Reject()` a nonzero correction is an error.

Errors on operators of order three or more, on second derivatives outside a divergence-compatible
position (inside a nonlinearity, two in one product, or a coefficient depending on `∇u`), and on a
non-differentiable coefficient that would need the product rule.
"""
function split_divergence(spatial, sivs, policy::IBPPolicy = ProductRule())
    sivs = unwrap.(sivs)
    G = Num[0 for _ in sivs]
    R = Num[]
    needy, corrections = Num[], Num[]
    for term in unwrap.(spatial)
        n = spatial_order(term, sivs)
        n > 2 && error("$(ordinal(n))-order operator needs a mixed or C¹ formulation: $term")
        fs = factors(term)
        idx = findall(f -> is_spatial_derivative(f, sivs) && spatial_order(f, sivs) >= 2, fs)
        if isempty(idx)
            spatial_order(term, sivs) >= 2 && not_in_divergence_form(term)
            push!(R, term)
            continue
        end
        length(idx) == 1 || not_in_divergence_form(term)
        Dw = fs[only(idx)]
        c = prod(fs[j] for j in eachindex(fs) if j != only(idx); init = 1)
        spatial_order(c, sivs) == 0 || not_in_divergence_form(term)
        op, w = operation(Dw), only(arguments(Dw))
        inner = op.order == 1 ? w : unwrap(Differential(op.x, op.order - 1)(w))
        G[findfirst(isequal(op.x), sivs)] += c * inner
        depends_on(c, op.x) || continue
        dc = Symbolics.expand_derivatives(Differential(op.x)(c))
        (is_smooth(c, op.x) && !is_piecewise(dc, op.x)) ||
            error("coefficient $c is not differentiable; write the equation in divergence form or register the derivative")
        push!(needy, term)
        push!(corrections, -dc * inner)
    end
    # Corrections of terms split apart by `expand` may cancel, as for (sin² + cos²)⋅Δu.
    isequal(Symbolics.simplify(sum(corrections; init = Num(0))), Num(0)) && return G, R, false
    policy isa Reject &&
        error("$(sum(needy)) needs the product rule, which the Reject() integration-by-parts policy forbids")
    return G, vcat(R, corrections), true
end
