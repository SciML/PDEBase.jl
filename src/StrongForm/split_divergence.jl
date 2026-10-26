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

"""
    DISCONTINUOUS_FUNCTIONS

The [`NONSMOOTH_FUNCTIONS`](@ref) that can jump: all of them except the continuous `abs`, `max`,
`min` and `clamp`. A derivative that gets expanded may act on a continuous non-smooth function
(the result is correct almost everywhere) but not on a jump, whose delta it would drop.
"""
const DISCONTINUOUS_FUNCTIONS = Tuple(f for f in NONSMOOTH_FUNCTIONS if f ∉ (abs, max, min, clamp))

"""
Is `c`, as a function of `x`, built only from functions with symbolic derivatives and none of
`excluded`?
"""
function is_smooth(c, x, excluded = NONSMOOTH_FUNCTIONS)
    depends_on(c, x) || return true
    iscall(c) || return true
    op, args = operation(c), arguments(c)
    any(f -> op === f, excluded) && return false
    op isa Function && op !== getindex &&
        !has_derivative_rules(op, args, a -> depends_on(a, x)) && return false
    return all(a -> is_smooth(a, x, excluded), args)
end

has_nan(ex) = iscall(ex) ? any(has_nan, arguments(ex)) :
    (v = unwrap_const(ex); v isa AbstractFloat && isnan(v))

"Does `ex` hold an `ifelse` depending on `x`, or a `NaN`, as derivatives of jumps do?"
function is_piecewise(ex, x)
    iscall(ex) || return has_nan(ex)
    operation(ex) === ifelse && depends_on(ex, x) && return true
    return any(a -> is_piecewise(a, x), arguments(ex))
end

not_in_divergence_form(term) =
    error("second derivative not in a divergence-compatible position: $term")

function not_differentiable(c, depvars)
    what = is_state(c, depvars) ? "expression" : "coefficient"
    return error("$what $c is not differentiable; write the equation in divergence form or register the derivative")
end

derivative_in(c, x) = Symbolics.expand_derivatives(Differential(x)(c))

"Smooth enough for a product-rule correction: free of every [`NONSMOOTH_FUNCTIONS`](@ref)."
differentiable_in(c, x) = is_smooth(c, x) && !is_piecewise(derivative_in(c, x), x)

"Continuous enough to be differentiated a.e.: free of [`DISCONTINUOUS_FUNCTIONS`](@ref)."
continuous_in(c, x) = is_smooth(c, x, DISCONTINUOUS_FUNCTIONS) && !has_nan(derivative_in(c, x))

"""
    check_expandable(ex, sivs, depvars = ())

Error unless every spatial derivative in `ex`, all of which get expanded when `ex` is put into
slots, acts on an argument without jumps: a jump's delta would be dropped silently.
"""
function check_expandable(ex, sivs, depvars = ())
    ex = unwrap(ex)
    iscall(ex) || return nothing
    if is_spatial_derivative(ex, sivs)
        arg = only(arguments(ex))
        continuous_in(arg, operation(ex).x) || not_differentiable(arg, depvars)
    end
    foreach(a -> check_expandable(a, sivs, depvars), arguments(ex))
    return nothing
end

"""
    split_divergence(spatial, sivs, policy::IBPPolicy = ProductRule(); depvars = ()) -> (G, R, corrected)

Stage b of the strong-form analysis. Writes the sum of the spatial terms from [`split_time`](@ref)
as `Σᵢ ∂ᵢG[i] + sum(R)` over the spatial variables `sivs`. A term without a second-order factor
goes to `R` unchanged. A term `c⋅Dᵢ(w′)` with one second-order factor adds `c⋅w′` to `G[i]` and,
when `c` depends on `xᵢ`, the product-rule correction `-(Dᵢc)⋅w′` to `R`, unless the corrections
of all terms simplify to zero together; `corrected` reports whether corrections were added.
Under `Reject()` a nonzero correction is an error.

Errors on operators of order three or more, on second derivatives outside a divergence-compatible
position (inside a nonlinearity, two in one product, or a coefficient depending on `∇u`), and on a
non-differentiable coefficient that would need the product rule (any of
[`NONSMOOTH_FUNCTIONS`](@ref)), and on a jump ([`DISCONTINUOUS_FUNCTIONS`](@ref)) under any other
derivative that will be expanded (first-order terms, the inner `w′`): only `c⋅w′` in the flux
stays undifferentiated, so a jump is accepted there alone. `depvars` only words the messages.
"""
function split_divergence(spatial, sivs, policy::IBPPolicy = ProductRule(); depvars = ())
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
            check_expandable(term, sivs, depvars)
            push!(R, term)
            continue
        end
        length(idx) == 1 || not_in_divergence_form(term)
        Dw = fs[only(idx)]
        c = prod(fs[j] for j in eachindex(fs) if j != only(idx); init = 1)
        spatial_order(c, sivs) == 0 || not_in_divergence_form(term)
        op, w = operation(Dw), only(arguments(Dw))
        inner = op.order == 1 ? w : unwrap(Differential(op.x, op.order - 1)(w))
        check_expandable(inner, sivs, depvars)
        G[findfirst(isequal(op.x), sivs)] += c * inner
        depends_on(c, op.x) || continue
        differentiable_in(c, op.x) || not_differentiable(c, depvars)
        push!(needy, term)
        push!(corrections, -Symbolics.expand_derivatives(Differential(op.x)(c)) * inner)
    end
    # Corrections of terms split apart by `expand` may cancel, as for (sin² + cos²)⋅Δu.
    isequal(Symbolics.simplify(sum(corrections; init = Num(0))), Num(0)) && return G, R, false
    policy isa Reject &&
        error("$(sum(needy)) needs the product rule, which the Reject() integration-by-parts policy forbids")
    return G, vcat(R, corrections), true
end
