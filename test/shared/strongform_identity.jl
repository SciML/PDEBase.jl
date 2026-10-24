# Reconstruction identity for PDEBase.StrongForm analyses, shared by the StrongForm tests and
# the PDESystemLibrary gate: strong[e] = mass[e] + Σᵢ ∂ᵢ flux[e][i] + remainder[e].
using PDEBase.StrongForm: StrongFormAnalysis, StateSlots
using Symbolics: Symbolics, Num, substitute

"""
    total_derivative(ex, i, slots, dim)

`∂ᵢ ex` for `ex` in slot form, by the chain rule through the state:
`∂/∂Xᵢ + Σ_f ∂/∂U_f ⋅ G[f, i] + Σ_{f, d} ∂/∂G[f, d] ⋅ H[f, d, i]`.
"""
function total_derivative(ex, i, slots::StateSlots, dim)
    out = Symbolics.derivative(ex, slots.X[i])
    for f in eachindex(slots.U)
        out += Symbolics.derivative(ex, slots.U[f]) * slots.G[(f - 1) * dim + i]
        for d in 1:dim
            out += Symbolics.derivative(ex, slots.G[(f - 1) * dim + d]) *
                slots.H[((f - 1) * dim + d - 1) * dim + i]
        end
    end
    return out
end

"""
    random_slot_values(a::StrongFormAnalysis) -> Dict

Values in [0.5, 1.5] for every slot, with a symmetric Hessian so that the order of mixed
partials does not matter.
"""
function random_slot_values(a::StrongFormAnalysis)
    s, dim = a.slots, a.dim
    vals = Dict{Num, Float64}(x => 0.5 + rand() for x in vcat(s.U, s.G, s.Ud, s.X, s.N, s.P, [s.T]))
    for f in 1:(a.nfields), d in 1:dim, e in d:dim
        h = 0.5 + rand()
        vals[s.H[((f - 1) * dim + d - 1) * dim + e]] = h
        vals[s.H[((f - 1) * dim + e - 1) * dim + d]] = h
    end
    return vals
end

slot_value(ex, vals) = Float64(Symbolics.value(substitute(Num(ex), vals; fold = Val(true))))

"""
    reconstruction_error(a::StrongFormAnalysis, e; nsamples = 5) -> Float64

Largest relative mismatch of `strong[e] = mass[e] + Σᵢ ∂ᵢ flux[e][i] + remainder[e]` over
`nsamples` random slot values.
"""
function reconstruction_error(a::StrongFormAnalysis, e; nsamples = 5)
    rebuilt = a.mass[e] + a.remainder[e] +
        sum(total_derivative(a.flux[e][i], i, a.slots, a.dim) for i in 1:(a.dim))
    return maximum(1:nsamples) do _
        vals = random_slot_values(a)
        strong = slot_value(a.strong[e], vals)
        abs(strong - slot_value(rebuilt, vals)) / max(1.0, abs(strong))
    end
end
