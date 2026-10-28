using PDEBase
using ModelingToolkit
using Test

struct XbarOrderDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
end
PDEBase.get_time(disc::XbarOrderDiscretization) = disc.time

@testset "VariableMap x̄ does not depend on equation order" begin
    @parameters t x y
    @variables u(..) q(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)
    eqs = [
        Dt(u(t, x, y)) ~ Dx(Dx(u(t, x, y))) + Dy(Dy(u(t, x, y))),
        Dt(q(t, y)) ~ Dy(Dy(q(t, y))),
    ]
    # `u(t, 0, y)` is a boundary call of `u` with a numeric argument; it must not
    # be able to put `y` ahead of `x`.
    bcs = [
        u(0, x, y) ~ 0, q(0, y) ~ 0, u(t, 0, y) ~ 0, u(t, 1, y) ~ 0,
        u(t, x, 0) ~ 0, u(t, x, 1) ~ 0, q(t, 0) ~ 0, q(t, 1) ~ 0,
    ]
    domains = [t ∈ (0, 1), x ∈ (0, 1), y ∈ (0, 1)]
    # The documented order is the position in `dvs`, full calls first, so `u` first
    # gives [x, y] and `q` first gives [y, x].
    expected = Dict([u(t, x, y), q(t, y)] => [x, y], [q(t, y), u(t, x, y)] => [y, x])

    # The identity, the reverse, every rotation, and the boundary calls of `u` ahead of its full call.
    n = length(bcs)
    bcperms = [collect(1:n), collect(n:-1:1), [circshift(1:n, k) for k in 1:(n - 1)]..., [3, 4, 5, 6, 1, 2, 7, 8]]
    eqperms = [[1, 2], [2, 1]]

    for (dvs, xbar) in expected
        results = map(bcperms) do bcperm
            map(eqperms) do eqperm
                pdesys = PDESystem(
                    eqs[eqperm], bcs[bcperm], domains, [t, x, y], dvs;
                    name = :pdebase_xbar_order_test
                )
                v = VariableMap(pdesys, XbarOrderDiscretization(t))
                (string.(indvars(v)), string.(v.x̄))
            end
        end
        flat = reduce(vcat, results)
        @test all(r -> r[1] == string.(xbar), flat)
        @test all(r -> r[2] == string.(xbar), flat)
    end
end
