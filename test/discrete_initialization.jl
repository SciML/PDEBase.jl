using PDEBase
using ModelingToolkit
using Symbolics
using SymbolicUtils
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

@testset "Discrete initialization with array (slice-form) equations" begin
    @independent_variables t
    @variables u(t)[1:5]
    D = Differential(t)
    bcast(op, args...) = Broadcast.materialize(Broadcast.broadcasted(op, args...))

    # Interior as a single slice-form array equation, as emitted by MethodOfLines'
    # ArrayDiscretization; the differentiated variable is the slice u[2:4].
    inteq = D(u[2:4]) ~ bcast(-, u[1:3], u[2:4])
    eqs = [inteq, u[1] ~ 0.0, u[5] ~ 0.0]
    u0 = [u[i] => Float64(i) for i in 1:5]

    init_eqs, guesses = PDEBase._discrete_initialization(eqs, t, u0)

    # The differentiated slice must be expanded to its elements so that elementwise
    # initial conditions generate initialization equations.
    @test Set(PDEBase.safe_unwrap(eq.lhs) for eq in init_eqs) ==
        Set(PDEBase.safe_unwrap(u[i]) for i in 2:4)
    @test guesses[PDEBase.safe_unwrap(u)] == [1.0, 2.0, 3.0, 4.0, 5.0]
end

@testset "Nonlinear normalization of array equations" begin
    @variables w[1:4]
    bcast(op, args...) = Broadcast.materialize(Broadcast.broadcasted(op, args...))

    iszeroval(x) = iszero(SymbolicUtils.unwrap_const(PDEBase.safe_unwrap(x)))

    arr_eq = bcast(-, w[2:4], w[1:3]) ~ zeros(3)
    neq = PDEBase._normalize_nonlinear_eq(arr_eq)
    @test size(neq.lhs) == (3,)
    @test all(iszeroval, collect(neq.lhs))
    @test PDEBase._is_array_valued(neq.rhs)

    scalar_eq = PDEBase.safe_unwrap(w[1]) ~ 1.0
    sneq = PDEBase._normalize_nonlinear_eq(scalar_eq)
    @test !PDEBase._is_array_valued(sneq.lhs)
    @test iszeroval(sneq.lhs)
end
