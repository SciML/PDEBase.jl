using PDEBase
using ModelingToolkit
using Symbolics
using Test

@testset "Discrete initialization" begin
    @independent_variables t
    @parameters p
    @variables u(t)[1:3] v(t)[1:2]
    D = Differential(t)

    eqs = [
        D(u[1]) ~ 0,
        D(u[2]) ~ v[1],
        D(u[3]) ~ p,
        v[1] ~ u[1],
        v[2] ~ u[2],
    ]
    u0 = [
        u[1] => 1.0,
        u[2] => p,
        u[3] => 3,
        v[1] => p + 1,
        v[2] => 5.0,
    ]

    init_eqs, guesses = PDEBase._discrete_initialization(eqs, t, u0)

    @test length(init_eqs) == 3
    @test Set(PDEBase.safe_unwrap(eq.lhs) for eq in init_eqs) ==
        Set(PDEBase.safe_unwrap(x) for x in collect(u))
    u2_init = only(eq for eq in init_eqs if isequal(eq.lhs, PDEBase.safe_unwrap(u[2])))
    @test isequal(u2_init.rhs, PDEBase.safe_unwrap(p))

    u_guesses = guesses[PDEBase.safe_unwrap(u)]
    v_guesses = guesses[PDEBase.safe_unwrap(v)]
    @test u_guesses[[1, 3]] == [1.0, 3]
    @test isequal(u_guesses[2], PDEBase.safe_unwrap(p))
    @test isequal(v_guesses[1], PDEBase.safe_unwrap(p + 1))
    @test v_guesses[2] == 5.0
    @test_nowarn System(
        eqs, t, [u, v], [p];
        initialization_eqs = init_eqs, guesses = guesses, name = :discrete_initialization_test
    )

    second_order_eqs = [(D^2)(u[1]) ~ 0]
    second_order_u0 = [u[1] => 1.0, D(u[1]) => p]
    second_order_init, derivative_guesses = PDEBase._discrete_initialization(
        second_order_eqs, t,
        [second_order_u0; D(u[2]) => 0.2; D(u[3]) => 0.3]
    )
    @test Set(PDEBase.safe_unwrap(eq.lhs) for eq in second_order_init) ==
        Set(PDEBase.safe_unwrap(first(pair)) for pair in second_order_u0)
    D_u = PDEBase.safe_unwrap(D(u))
    @test derivative_guesses[D_u][[2, 3]] == [0.2, 0.3]
    @test isequal(derivative_guesses[D_u][1], PDEBase.safe_unwrap(p))
end
