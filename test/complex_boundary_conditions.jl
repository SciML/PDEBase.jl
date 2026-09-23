using PDEBase
using ModelingToolkit
using Symbolics
using SymbolicUtils
using Test

@testset "Complex boundary condition splitting" begin
    @parameters t x
    @variables ψ(..) Reψ(..) Imψ(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Dx^2
    domain = [t ∈ (0, 1), x ∈ (0, 1)]

    function split(equations, bcs)
        pdesys = PDESystem(equations, bcs, domain, [t, x], [ψ(t, x)]; name = :complex_bc_test)
        return PDEBase.handle_complex(pdesys)[1]
    end

    function same_equations(actual, expected)
        residual(eq) = Symbolics.simplify(eq.lhs - eq.rhs)
        return length(actual) == length(expected) && all(eachindex(expected)) do i
            isequal(Symbolics.value(Symbolics.simplify(residual(actual[i]) - residual(expected[i]))), 0)
        end
    end

    @testset "Complex Neumann condition" begin
        system = split(
            [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [ψ(t, 0) ~ 0, Dx(ψ(t, 1)) ~ im * ψ(t, 1)]
        )
        expected = [
            Reψ(t, 0) ~ 0,
            Imψ(t, 0) ~ 0,
            Dx(Reψ(t, 1)) ~ -Imψ(t, 1),
            Dx(Imψ(t, 1)) ~ Reψ(t, 1),
        ]
        @test same_equations(PDEBase.get_bcs(system), expected)
    end

    @testset "Complex Robin condition" begin
        system = split(
            [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [ψ(t, 0) ~ 0, Dx(ψ(t, 1)) + im * ψ(t, 1) ~ 0]
        )
        expected = [
            Reψ(t, 0) ~ 0,
            Imψ(t, 0) ~ 0,
            Dx(Reψ(t, 1)) ~ Imψ(t, 1),
            Dx(Imψ(t, 1)) ~ -Reψ(t, 1),
        ]
        @test same_equations(PDEBase.get_bcs(system), expected)
    end

    @testset "Complex BC reduced to false is rejected" begin
        @test_throws ArgumentError split(
            [Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [ψ(t, 0) ~ im, ψ(t, 1) ~ 0]
        )
    end

    @testset "Equation forms retain exact real boundary splits" begin
        equations = [
            im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x)) + 0.0 * ψ(t, x),
            (2 + 3im) * Dt(ψ(t, x)) ~ (1 - im) * Dxx(ψ(t, x)),
            Dt(ψ(t, x)) ~ Dxx(ψ(t, x)) / (1 + 2im),
            im * Dt(ψ(t, x)) ~ ψ(t, x)^2,
            im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x)) + abs2(ψ(t, x)) * ψ(t, x),
            im * Dt(ψ(t, x)) ~ conj(ψ(t, x)),
            im * Dt(ψ(t, x)) ~ exp(ψ(t, x)),
            im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x)) + im * sin(x) * ψ(t, x),
        ]
        expected = [Reψ(t, 0) ~ 0, Imψ(t, 0) ~ 0, Reψ(t, 1) ~ 0, Imψ(t, 1) ~ 0]

        for eq in equations
            system = split([eq], [ψ(t, 0) ~ 0, ψ(t, 1) ~ 0])
            @test same_equations(PDEBase.get_bcs(system), expected)
        end
    end

    @testset "Complex initial and boundary values split exactly" begin
        exp0 = SymbolicUtils.term(exp, 0)
        system = split(
            [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [ψ(0, x) => exp(im * x), ψ(t, 0) ~ exp(im * t), ψ(t, 1) ~ (1 + im) * t]
        )
        expected = [
            Reψ(0, x) ~ exp0 * cos(x),
            Imψ(0, x) ~ exp0 * sin(x),
            Reψ(t, 0) ~ exp0 * cos(t),
            Imψ(t, 0) ~ exp0 * sin(t),
            Reψ(t, 1) ~ t,
            Imψ(t, 1) ~ t,
        ]
        @test same_equations(PDEBase.get_bcs(system), expected)
    end

    @testset "Multiple dependent variables retain exact boundary splits" begin
        @variables χ(..) Reχ(..) Imχ(..)
        system = PDESystem(
            [im * Dt(ψ(t, x)) ~ (1 + 2im) * Dxx(χ(t, x)), im * Dt(χ(t, x)) ~ ψ(t, x) * χ(t, x)],
            [ψ(t, 0) ~ 0, ψ(t, 1) ~ 0, χ(t, 0) ~ 0, χ(t, 1) ~ 0],
            domain, [t, x], [ψ(t, x), χ(t, x)]; name = :complex_bc_test
        )
        split_system, _ = PDEBase.handle_complex(system)
        expected = [
            Reψ(t, 0) ~ 0,
            Imψ(t, 0) ~ 0,
            Reψ(t, 1) ~ 0,
            Imψ(t, 1) ~ 0,
            Reχ(t, 0) ~ 0,
            Imχ(t, 0) ~ 0,
            Reχ(t, 1) ~ 0,
            Imχ(t, 1) ~ 0,
        ]
        @test same_equations(PDEBase.get_bcs(split_system), expected)
    end
end
