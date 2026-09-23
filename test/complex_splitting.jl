using PDEBase
using ModelingToolkit
using Symbolics
using Test

@testset "Complex PDE equations retain real-imaginary coupling" begin
    @parameters t x
    @variables ψ(..) χ(..) Reψ(..) Imψ(..) Reχ(..) Imχ(..)
    Dt = Differential(t)
    Dxx = Differential(x)^2

    pdesys = PDESystem(
        [
            im * Dt(ψ(t, x)) ~ (1 + 2im) * Dxx(χ(t, x)),
            im * Dt(χ(t, x)) ~ ψ(t, x) * χ(t, x),
        ],
        [ψ(t, 0) ~ 0, ψ(t, 1) ~ 0, χ(t, 0) ~ 0, χ(t, 1) ~ 0],
        [t ∈ (0, 1), x ∈ (0, 1)], [t, x], [ψ(t, x), χ(t, x)];
        name = :complex_splitting_test
    )

    split, _ = PDEBase.handle_complex(pdesys)
    equations = ModelingToolkit.equations(split)
    expected = [
        -Dt(Imψ(t, x)) ~ Dxx(Reχ(t, x)) - 2Dxx(Imχ(t, x)),
        Dt(Reψ(t, x)) ~ Dxx(Imχ(t, x)) + 2Dxx(Reχ(t, x)),
        -Dt(Imχ(t, x)) ~ Reψ(t, x) * Reχ(t, x) - Imψ(t, x) * Imχ(t, x),
        Dt(Reχ(t, x)) ~ Reψ(t, x) * Imχ(t, x) + Imψ(t, x) * Reχ(t, x),
    ]

    residual(eq) = Symbolics.simplify(eq.lhs - eq.rhs)
    @test length(equations) == length(expected)
    @test all(
        isequal(Symbolics.value(Symbolics.simplify(residual(equations[i]) - residual(expected[i]))), 0)
            for i in eachindex(expected)
    )
    @test length(PDEBase.get_bcs(split)) == 8

    complex_bc_system = PDESystem(
        [Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
        [ψ(t, 0) ~ cos(t) + im * sin(t)],
        [t ∈ (0, 1), x ∈ (0, 1)], [t, x], [ψ(t, x)];
        name = :complex_dirichlet_test
    )
    complex_bc = only(PDEBase._flatten_bcs(PDEBase.get_bcs(complex_bc_system)))
    redvmaps = Dict(operation(unwrap(ψ(t, 0))) => operation(unwrap(Reψ(t, 0))))
    imdvmaps = Dict(operation(unwrap(ψ(t, 0))) => operation(unwrap(Imψ(t, 0))))
    split_bc = PDEBase.split_complex_bc(complex_bc, redvmaps, imdvmaps)
    expected_bc = [Reψ(t, 0) ~ cos(t), Imψ(t, 0) ~ sin(t)]
    @test all(
        isequal(Symbolics.value(Symbolics.simplify(residual(split_bc[i]) - residual(expected_bc[i]))), 0)
            for i in eachindex(expected_bc)
    )
end

@testset "Mixed complex-typed and untyped dependent variables" begin
    @parameters t x
    @variables ψ(..)::Complex u(..)
    Dt = Differential(t)
    z = ψ(t, x)
    v = u(t, x)
    pdesys = PDESystem(
        [im * Dt(z) ~ conj(z) * v, Dt(v) ~ abs2(z)],
        [ψ(0, x) => 1.0, u(0, x) => 0.0, ψ(t, 0) ~ 0, u(t, 0) ~ 0],
        [t ∈ (0, 1), x ∈ (0, 1)], [t, x], [z, v]; name = :mixed_typed_test
    )

    @test_throws "Complex-typed and real dependent variables cannot be mixed" PDEBase.handle_complex(pdesys)
end

@testset "Complex-typed dependent variable splitting" begin
    @parameters t x
    @variables ψ(..)::Complex Reψ(..) Imψ(..)
    Dt = Differential(t)
    z = ψ(t, x)
    a = Reψ(t, x)
    b = Imψ(t, x)

    residual(eq) = Symbolics.simplify(eq.lhs - eq.rhs)
    same_equation(actual, expected) = isequal(
        Symbolics.value(Symbolics.simplify(residual(actual) - residual(expected))), 0
    )
    stage_orders = (
        (PDEBase.handle_complex, PDEBase.make_pdesys_compatible),
        (PDEBase.make_pdesys_compatible, PDEBase.handle_complex),
    )

    for (nonlinearity, expected_nonlinearity) in (
            (z * conj(z) * z, (a^2 + b^2) * (a + im * b)),
            (abs2(z) * z, (a^2 + b^2) * (a + im * b)),
        )
        pdesys = PDESystem(
            [im * Dt(z) ~ nonlinearity],
            [ψ(0, x) => 1 + 2im, ψ(t, 0) ~ 3 + 4im],
            [t ∈ (0, 1), x ∈ (0, 1)], [t, x], [z]; name = :typed_complex_test
        )
        expected = [
            -Dt(b) ~ real(expected_nonlinearity),
            Dt(a) ~ imag(expected_nonlinearity),
        ]
        expected_bcs = [
            Reψ(0, x) ~ 1, Imψ(0, x) ~ 2,
            Reψ(t, 0) ~ 3, Imψ(t, 0) ~ 4,
        ]
        for (first_stage, second_stage) in stage_orders
            intermediate, _ = first_stage(pdesys)
            split, _ = second_stage(intermediate)
            equations = ModelingToolkit.equations(split)
            @test length(equations) == 2
            @test all(same_equation(equations[i], expected[i]) for i in eachindex(expected))
            @test length(PDEBase.get_dvs(split)) == 2
            bcs = PDEBase._flatten_bcs(PDEBase.get_bcs(split))
            @test length(bcs) == length(expected_bcs)
            @test all(any(same_equation(bc, expected_bc) for bc in bcs) for expected_bc in expected_bcs)
        end
    end

    pdesys = PDESystem(
        [Dt(z) ~ conj(z) + real(z) + im * imag(z) + abs2(z)],
        [ψ(t, 0) ~ 0], [t ∈ (0, 1), x ∈ (0, 1)], [t, x], [z];
        name = :typed_complex_components_test
    )
    expected = [Dt(a) ~ 2a + a^2 + b^2, Dt(b) ~ 0]
    for (first_stage, second_stage) in stage_orders
        intermediate, _ = first_stage(pdesys)
        split, _ = second_stage(intermediate)
        equations = ModelingToolkit.equations(split)
        @test length(equations) == 2
        @test all(same_equation(equations[i], expected[i]) for i in eachindex(expected))
    end
end
