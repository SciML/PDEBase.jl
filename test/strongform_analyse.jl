using PDEBase
using PDEBase.StrongForm
using PDEBase.StrongForm: ProductRule, Reject, Weighted, check_evaluable, split_divergence
using IntervalSets: var".."
using ModelingToolkit: @named, @parameters, PDESystem, get_bcs
using Symbolics
using SymbolicUtils: operation
using Test

include(joinpath(@__DIR__, "shared", "strongform_identity.jl"))

struct TimeDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
end
PDEBase.get_time(d::TimeDiscretization) = d.time

function prepared(eqs, dvs, ivs; t = nothing, ps = Num[], bcs = Equation[])
    @named sys = PDESystem(eqs, bcs, [iv ∈ (0.1 .. 1.0) for iv in ivs], ivs, dvs, ps)
    # The preprocessing `symbolic_discretize` applies before the discretizer hooks.
    PDEBase.cardinalize_eqs!(sys)
    sys, replaced = PDEBase.make_pdesys_compatible(sys)
    return sys, PDEBase.VariableMap(sys, TimeDiscretization(t); replaced_vars = replaced)
end

function analysed(eqs, dvs, ivs; ibp = Dict(), kw...)
    sys, v = prepared(eqs, dvs, ivs; kw...)
    return analyse(sys, v, Dict(); ibp), v
end

same(a, b) = isequal(Symbolics.expand(Num(a) - Num(b)), Num(0))

opaque_coefficient(x) = 1 + x^2
@register_symbolic opaque_coefficient(x)

@parameters k
@variables t x y
@variables u(..) w(..) a(..)
Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
U, W = u(t, x, y), w(t, x, y)

@testset "Cases (i), (ii), (ii-b), (ii-c) in slot form" begin
    # (i) constant coefficient
    an, _ = @test_logs analysed([Dt(U) ~ Dx(k * Dx(U)) + Dy(k * Dy(U)) + x * y + 1], [U], [t, x, y]; t, ps = [k])
    s = an.slots
    X, G, N, P = s.X, s.G, s.N, s.P
    @test (an.nfields, an.dim) == (1, 2)
    @test same(an.mass[1], s.Ud[1])
    @test all(same.(an.flux[1], [-P[1] * G[1], -P[1] * G[2]]))
    @test same(an.remainder[1], -X[1] * X[2] - 1)
    @test all(same.(an.conormal[1], [-P[1] * N[1], -P[1] * N[2]]))
    @test same(an.strong[1], s.Ud[1] - P[1] * s.H[1] - P[1] * s.H[4] - X[1] * X[2] - 1)
    @test reconstruction_error(an, 1) < 1.0e-12

    # (ii) variable coefficient in divergence form
    kxy = 1 + x^2 * y
    an, _ = @test_logs analysed([Dt(U) ~ Dx(kxy * Dx(U)) + Dy(kxy * Dy(U)) + x * y + 1], [U], [t, x, y]; t)
    X, G = an.slots.X, an.slots.G
    @test all(same.(an.flux[1], -(1 + X[1]^2 * X[2]) .* G))
    @test same(an.remainder[1], -X[1] * X[2] - 1)
    @test reconstruction_error(an, 1) < 1.0e-12

    # (ii-b) non-divergence a(x, y)Δu with advection and reaction: the product rule adds ∇a⋅∇u
    axy = 1 + x^2 + y
    eq = Dt(U) ~ axy * (Dx(Dx(U)) + Dy(Dy(U))) - 0.7 * Dx(U) - 0.3 * Dy(U) - 2.5 * U + x * y + 1
    an, _ = @test_logs (:warn, r"equations 1 \(ProductRule\), so .* non-symmetric") analysed([eq], [U], [t, x, y]; t)
    s = an.slots
    X, G = s.X, s.G
    @test all(same.(an.flux[1], -(1 + X[1]^2 + X[2]) .* G))
    @test same(an.remainder[1], 2X[1] * G[1] + G[2] + 0.7G[1] + 0.3G[2] + 2.5s.U[1] - X[1] * X[2] - 1)
    @test reconstruction_error(an, 1) < 1.0e-12
    @test_throws r"needs the product rule, which the Reject\(\) integration-by-parts policy forbids" analysed(
        [eq], [U], [t, x, y]; t, ibp = Dict(1 => Reject())
    )

    # (ii-c) anisotropic tensor with cross terms
    k11, k12, k22 = 2 + x, 0.4 + 0.1y, 1 + y^2
    eq = Dt(U) ~ Dx(k11 * Dx(U) + k12 * Dy(U)) + Dy(k12 * Dx(U) + k22 * Dy(U)) + x * y + 1
    an, _ = @test_logs analysed([eq], [U], [t, x, y]; t, ibp = Dict(eq => Reject()))
    X, G = an.slots.X, an.slots.G
    K11, K12, K22 = 2 + X[1], 0.4 + 0.1X[2], 1 + X[2]^2
    @test all(same.(an.flux[1], [-(K11 * G[1] + K12 * G[2]), -(K12 * G[1] + K22 * G[2])]))
    @test reconstruction_error(an, 1) < 1.0e-12
end

@testset "Case (iii) and the collocation strong kernel at the exact solution" begin
    V = u(x, y)
    ustar = sin(π * x) * cos(π * y) + x^2 * y
    kstar = 1 + ustar^2
    f = expand_derivatives(-(Dx(kstar * Dx(ustar)) + Dy(kstar * Dy(ustar))))
    an, _ = analysed([0 ~ Dx((1 + V^2) * Dx(V)) + Dy((1 + V^2) * Dy(V)) + f], [V], [x, y])
    s = an.slots
    @test all(same.(an.flux[1], -(1 + s.U[1]^2) .* s.G))
    @test isequal(an.mass[1], Num(0))
    @test reconstruction_error(an, 1) < 1.0e-12

    d(ex, ivs...) = foldl((acc, iv) -> expand_derivatives(Differential(iv)(acc)), ivs; init = ex)
    exact = [ustar, d(ustar, x), d(ustar, y), d(ustar, x, x), d(ustar, x, y), d(ustar, y, x), d(ustar, y, y)]
    worst = maximum(1:20) do _
        point = Dict(x => rand(), y => rand())
        u0, gx, gy, hxx, hxy, hyx, hyy = [slot_value(substitute(ex, point), Dict()) for ex in exact]
        vals = Dict(
            s.X[1] => point[x], s.X[2] => point[y], s.U[1] => u0, s.G[1] => gx, s.G[2] => gy,
            s.H[1] => hxx, s.H[2] => hxy, s.H[3] => hyx, s.H[4] => hyy,
        )
        abs(slot_value(an.strong[1], vals))
    end
    @test worst <= 1.0e-13
end

@testset "Appendix A accepted rows" begin
    # Second-order mixed partials
    an, _ = analysed([Dt(U) ~ Dx(Dy(U)) + Dx(Dx(U)) + Dy(Dy(U))], [U], [t, x, y]; t)
    G = an.slots.G
    @test all(same.(an.flux[1], [-G[2] - G[1], -G[2]]))
    @test reconstruction_error(an, 1) < 1.0e-12

    # Multi-field with a coupled flux ∇⋅(D_u∇u + χ∇w)
    @parameters du χ
    eqs = [
        Dt(U) ~ Dx(du * Dx(U) + χ * Dx(W)) + Dy(du * Dy(U) + χ * Dy(W)),
        Dt(W) ~ Dx(Dx(W)) + Dy(Dy(W)) - U * W,
    ]
    an, v = analysed(eqs, [U, W], [t, x, y]; t, ps = [du, χ])
    s = an.slots
    fu, fw = 1, 2  # declared order
    gu, gw = s.G[(2fu - 1):(2fu)], s.G[(2fw - 1):(2fw)]
    @test an.nfields == 2
    @test all(same.(an.flux[1], -s.P[1] .* gu .- s.P[2] .* gw))
    @test same(an.mass[1], s.Ud[fu]) && same(an.mass[2], s.Ud[fw])
    @test length(an.conormal[1]) == 4
    @test all(e -> reconstruction_error(an, e) < 1.0e-12, 1:2)

    # Algebraic row and wave equation as a first-order system, in one dimension
    U1, W1 = u(t, x), w(t, x)
    an, v = analysed([Dt(U1) ~ Dx(Dx(U1)) + W1, 0 ~ Dx(Dx(W1)) + exp(-t) * sin(x)], [U1, W1], [t, x]; t)
    @test isequal(an.mass[2], Num(0))
    @test same(an.remainder[2], -exp(-an.slots.T) * sin(an.slots.X[1]))
    an, v = analysed([Dt(U1) ~ W1, Dt(W1) ~ Dx(Dx(U1))], [U1, W1], [t, x]; t)
    s = an.slots
    @test same(an.mass[1], s.Ud[1]) && same(an.mass[2], s.Ud[2])
    @test same(only(an.flux[2]), -s.G[1])
    @test all(e -> reconstruction_error(an, e) < 1.0e-12, 1:2)

    # Time-dependent coefficient D(t, x)Δu: a timedep flux and a product-rule correction
    an, _ = @test_logs (:warn,) analysed([Dt(U1) ~ (1 + t * x^2) * Dx(Dx(U1))], [U1], [t, x]; t)
    s = an.slots
    @test same(only(an.flux[1]), -(1 + s.T * s.X[1]^2) * s.G[1])
    @test same(an.remainder[1], 2s.T * s.X[1] * s.G[1])

    # ifelse coefficient in divergence form
    an, _ = @test_logs analysed([Dt(U1) ~ Dx(ifelse(x < 0.5, 1.0, 3.0) * Dx(U1))], [U1], [t, x]; t)
    @test same(only(an.flux[1]), -ifelse(an.slots.X[1] < 0.5, 1.0, 3.0) * an.slots.G[1])
end

@testset "Rejections" begin
    U1 = u(t, x)
    reject(eq; kw...) = analysed([eq], [U1], [t, x]; t, kw...)
    @test_throws r"^third-order operator needs a mixed or C¹ formulation" reject(
        Dt(U1) + U1 * Dx(U1) + Dx(Dx(Dx(U1))) ~ 0
    )
    @test_throws r"^fourth-order operator needs a mixed or C¹ formulation" reject(Dt(U1) ~ -Dx(Dx(Dx(Dx(U1)))))
    msg = "second derivative not in a divergence-compatible position"
    @test_throws msg reject(Dt(U1) ~ sin(Dx(Dx(U1))))
    @test_throws msg reject(Dt(U1) ~ Dx(U1) * Dx(Dx(U1)))
    @test_throws msg analysed([Dt(U) ~ Dx(Dx(U)) * Dy(Dy(U))], [U], [t, x, y]; t)
    @test_throws r"is not differentiable; write the equation in divergence form or register the derivative" reject(
        Dt(U1) ~ ifelse(x < 0.5, 1.0, 3.0) * Dx(Dx(U1))
    )
    @test_throws r"unresolved symbols .* register the derivative" reject(Dt(U1) ~ a(x) * Dx(Dx(U1)))

    # The true order wins over the divergence-position message.
    @test_throws r"^third-order operator" reject(Dt(U1) ~ sin(Dx(Dx(Dx(U1)))))
    @test_throws r"^third-order operator" reject(Dt(U1) ~ Dx(Dx(U1)) * Dx(Dx(Dx(U1))))

    # Discontinuous coefficients, listed or caught through a piecewise or NaN derivative.
    for c in (mod(x, 0.5), rem(x, 0.5), clamp(x, 0.2, 0.8))
        @test_throws "is not differentiable" reject(Dt(U1) ~ c * Dx(Dx(U1)))
    end
    # A jump under a derivative that gets expanded would lose its delta; only the flux keeps it.
    jump = ifelse(x < 0.5, 1.0, 3.0)
    for rhs in (Dx(jump * U1), Dx(jump) * Dx(U1), Dx(Dx(jump * U1)))
        @test_throws "is not differentiable; write the equation in divergence form" reject(Dt(U1) ~ rhs)
    end
    @test reject(Dt(U1) ~ Dx(jump * Dx(U1)))[1] isa StrongFormAnalysis
    # The same holds for a mass coefficient.
    @test_throws "is not differentiable" reject(Dx(jump) * Dt(U1) ~ Dx(Dx(U1)))
    @test_throws "is not differentiable" reject(Dt(U1) + Dx(jump) * Dt(U1) ~ 0)
    @test_throws r"^expression .* is not differentiable" reject(Dt(U1) ~ Dx(sign(U1)) + Dx(Dx(U1)))
    # Continuous non-smooth coefficients still cannot take the product rule (M1).
    @test_throws r"^coefficient .* is not differentiable" reject(Dt(U1) ~ abs(x - 0.5) * Dx(Dx(U1)))
    # Continuous non-smooth fluxes of the state may be differentiated; correct almost everywhere.
    for flux in (abs(U1), U1 * abs(U1), max(U1, 0)^2)
        an, _ = @test_logs reject(Dt(U1) ~ Dx(flux) + Dx(Dx(U1)))
        @test reconstruction_error(an, 1) < 1.0e-12
    end
    is_piecewise = PDEBase.StrongForm.is_piecewise
    @test is_piecewise(Symbolics.unwrap(Symbolics.derivative(mod(x, 0.5), x)), Symbolics.unwrap(x))
    @test !is_piecewise(Symbolics.unwrap(Symbolics.derivative(sin(x) * ifelse(y < 0.5, 1, 2), x)), Symbolics.unwrap(x))
end

@testset "Weighted(r²) gives the symmetric spherical form" begin
    @variables r
    @variables q(..)
    Dr = Differential(r)
    Q = q(t, r)
    eq = Dt(Q) ~ 1 / r^2 * Dr(r^2 * Dr(Q))
    an, _ = @test_logs analysed([eq], [Q], [t, r]; t, ibp = Dict(eq => Weighted(r^2)))
    s = an.slots
    @test same(an.mass[1], s.X[1]^2 * s.Ud[1])
    @test same(only(an.flux[1]), -s.X[1]^2 * s.G[1])
    @test isequal(an.remainder[1], Num(0))
    an, _ = @test_logs (:warn, r"ProductRule") analysed([eq], [Q], [t, r]; t)
    @test same(only(an.flux[1]), -an.slots.G[1])
    @test same(an.remainder[1] * an.slots.X[1], -2an.slots.G[1])
end

@testset "Integration-by-parts policies keyed through preprocessing" begin
    U1 = u(t, x)
    eq = Dt(U1) ~ (1 + x^2) * Dx(Dx(U1))
    @test_throws "Reject() integration-by-parts policy forbids" analysed([eq], [U1], [t, x]; t, ibp = Dict(eq => Reject()))
    flipped = (1 + x^2) * Dx(Dx(U1)) ~ Dt(U1)
    @test_throws "Reject()" analysed([eq], [U1], [t, x]; t, ibp = Dict(flipped => Reject()))
    @test_throws "matches no equation of the system" analysed(
        [eq], [U1], [t, x]; t, ibp = Dict((Dt(U1) ~ Dx(Dx(U1))) => Reject())
    )
    @test_throws "matches no equation of the system" analysed([eq], [U1], [t, x]; t, ibp = Dict(2 => Reject()))
    @test_throws "ibp keys must be Equations or integer indices" analysed(
        [eq], [U1], [t, x]; t, ibp = Dict((Dt(U1) - (1 + x^2) * Dx(Dx(U1))) => Reject())
    )

    # A coefficient whose derivative simplifies to zero needs no correction.
    eq = Dt(U1) ~ (sin(x)^2 + cos(x)^2) * Dx(Dx(U1))
    an, _ = @test_logs analysed([eq], [U1], [t, x]; t, ibp = Dict(eq => Reject()))
    @test isequal(an.remainder[1], Num(0))
    G, R, corrected = split_divergence([(sin(x)^2 + cos(x)^2) * Dx(Dx(U1))], [x])
    @test isempty(R) && !corrected
end

@testset "More integration-by-parts cases" begin
    U1 = u(t, x)
    # Mixed partial whose coefficient depends on the outer differentiation variable.
    an, _ = @test_logs (:warn,) analysed([Dt(U) ~ x * Dx(Dy(U))], [U], [t, x, y]; t)
    s = an.slots
    @test all(same.(an.flux[1], [-s.X[1] * s.G[2], Num(0)]))
    @test same(an.remainder[1], s.G[2])
    @test reconstruction_error(an, 1) < 1.0e-12

    # The warning names the policy that generated the correction.
    eq = Dt(U1) ~ (1 + x^2) * Dx(Dx(U1))
    @test_logs (:warn, r"equations 1 \(Weighted") analysed([eq], [U1], [t, x]; t, ibp = Dict(eq => Weighted(x)))

    # Conflicting keys for one equation.
    @test_throws "matched by ibp keys with different policies" analysed(
        [eq], [U1], [t, x]; t, ibp = Dict(eq => Reject(), 1 => ProductRule())
    )
    @test_logs (:warn,) analysed([eq], [U1], [t, x]; t, ibp = Dict(eq => Weighted(x^2), 1 => Weighted(x^2)))

    # Keys written with array-variable fields match after flattening.
    @variables (q(..))[1:2]
    Q1, Q2 = q(t, x)[1], q(t, x)[2]
    eqs = [Dt(Q1) ~ (1 + x^2) * Dx(Dx(Q1)), Dt(Q2) ~ Dx(Dx(Q2))]
    @test_throws "Reject()" analysed(eqs, [q(t, x)], [t, x]; t, ibp = Dict(eqs[1] => Reject()))

    # Derivatives of constants are dropped when distributing.
    @test isequal(PDEBase.StrongForm.distribute_derivatives(Dx(U1 + 2)), Symbolics.unwrap(Dx(U1)))
end

@testset "Opaque coefficients in divergence form" begin
    U1 = u(t, x)
    an, _ = @test_logs analysed([Dt(U1) ~ Dx(opaque_coefficient(x) * Dx(U1))], [U1], [t, x]; t)
    s = an.slots
    @test same(only(an.flux[1]), -opaque_coefficient(s.X[1]) * s.G[1])
    @test isequal(an.remainder[1], Num(0))
    # The strong residual keeps the derivative of the opaque coefficient, in slot coordinates.
    rebuilt = an.mass[1] + total_derivative(only(an.flux[1]), 1, s, 1) + an.remainder[1]
    @test same(an.strong[1], rebuilt)
    @test_throws r"unevaluated derivative Differential\(X₁, 1\)\(.*opaque_coefficient\(X₁\)\)" check_evaluable(
        an.strong[1]
    )
    @test check_evaluable(only(an.flux[1])) === nothing
    @test all(ex -> check_evaluable(ex) === nothing, vcat(an.mass, an.remainder, an.conormal...))
    # In non-divergence form the product rule would need its derivative.
    @test_throws "is not differentiable" analysed([Dt(U1) ~ opaque_coefficient(x) * Dx(Dx(U1))], [U1], [t, x]; t)
end

@testset "Fields follow the declared order, not the order of discovery" begin
    U1, W1 = u(t, x), w(t, x)
    eqs = [Dt(W1) ~ Dx(Dx(W1)) + U1, Dt(U1) ~ Dx(2 * Dx(U1))]
    an, v = analysed(eqs, [U1, W1], [t, x]; t)
    @test isequal(operation.(Symbolics.unwrap.(depvars(v))), operation.(Symbolics.unwrap.([W1, U1])))
    s = an.slots
    @test same(an.mass[1], s.Ud[2]) && same(an.mass[2], s.Ud[1])
    @test same(only(an.flux[1]), -s.G[2]) && same(only(an.flux[2]), -2s.G[1])
    @test same(an.remainder[1], -s.U[1])
    @test all(e -> reconstruction_error(an, e) < 1.0e-12, 1:2)

    # Array variables keep their declared element order after flattening.
    @variables (q(..))[1:2]
    Q1, Q2 = q(t, x)[1], q(t, x)[2]
    an, _ = analysed([Dt(Q2) ~ Dx(Dx(Q2)) + Q1, Dt(Q1) ~ Dx(Dx(Q1))], [q(t, x)], [t, x]; t)
    s = an.slots
    @test same(an.mass[1], s.Ud[2]) && same(an.remainder[1], -s.U[1])
    @test all(e -> reconstruction_error(an, e) < 1.0e-12, 1:2)
end

@testset "Array parameters get one slot per element" begin
    U1 = u(t, x)
    @parameters ap[1:2]
    an, _ = analysed([Dt(U1) ~ ap[1] * Dx(Dx(U1)) + ap[2] * k], [U1], [t, x]; t, ps = [k, ap])
    s = an.slots
    @test length(s.P) == 3
    @test same(only(an.flux[1]), -s.P[2] * s.G[1])
    @test same(an.remainder[1], -s.P[3] * s.P[1])
    @test reconstruction_error(an, 1) < 1.0e-12
end

@testset "split_divergence on its own" begin
    G, R, corrected = split_divergence([Dx(Dx(U)), x * Dy(Dy(U))], [x, y])
    @test isequal(G, [Dx(U), x * Dy(U)])
    @test isempty(R) && !corrected
end

@testset "Boundary conditions are not analysed yet" begin
    U1 = u(t, x)
    sys, v = prepared([Dt(U1) ~ Dx(Dx(U1))], [U1], [t, x]; t, bcs = [u(0.1, x) ~ 0, u(t, 0.1) ~ 0])
    orders = Dict(iv => PDEBase.d_orders(iv, get_bcs(sys)) for iv in PDEBase.all_ivs(v))
    @test_throws "boundary-condition analysis is not implemented yet" analyse(
        sys, v, PDEBase.parse_bcs(get_bcs(sys), v, orders)
    )
    sys, v = prepared([Dt(U1) ~ Dx(Dx(U1))], [U1], [t, x]; t, bcs = [u(0.1, x) ~ 0])
    orders = Dict(iv => PDEBase.d_orders(iv, get_bcs(sys)) for iv in PDEBase.all_ivs(v))
    @test analyse(sys, v, PDEBase.parse_bcs(get_bcs(sys), v, orders)) isa StrongFormAnalysis
end
