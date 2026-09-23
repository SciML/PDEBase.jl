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
