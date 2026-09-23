using PDEBase
using ModelingToolkit
using Symbolics
using SymbolicUtils
using Test

@testset "Complex boundary condition splitting" begin
    @parameters t x
    @variables ψ(..) Reψ(..) Imψ(..)
    @variables ψc(..)::Complex Reψc(..) Imψc(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Dx^2
    domain = [t ∈ (0, 1), x ∈ (0, 1)]

    function split(equations, bcs)
        pdesys = PDESystem(equations, bcs, domain, [t, x], [ψ(t, x)]; name = :complex_bc_test)
        return PDEBase.handle_complex(pdesys)[1]
    end

    function split_typed(equations, bcs)
        pdesys = PDESystem(equations, bcs, domain, [t, x], [ψc(t, x)]; name = :typed_complex_bc_test)
        return PDEBase.handle_complex(pdesys)[1]
    end

    function same_equations(actual, expected)
        residual(eq) = Symbolics.simplify(eq.lhs - eq.rhs)
        return length(actual) == length(expected) && all(eachindex(expected)) do i
            isequal(Symbolics.value(Symbolics.simplify(residual(actual[i]) - residual(expected[i]))), 0)
        end
    end

    @testset "Complex Neumann condition" begin
        system = split_typed(
            [im * Dt(ψc(t, x)) ~ Dxx(ψc(t, x))],
            [ψc(t, 0) ~ 0, Dx(ψc(t, 1)) ~ im * ψc(t, 1)]
        )
        expected = [
            Reψc(t, 0) ~ 0,
            Imψc(t, 0) ~ 0,
            Dx(Reψc(t, 1)) ~ -Imψc(t, 1),
            Dx(Imψc(t, 1)) ~ Reψc(t, 1),
        ]
        @test same_equations(PDEBase.get_bcs(system), expected)
    end

    @testset "Complex Robin condition" begin
        system = split_typed(
            [im * Dt(ψc(t, x)) ~ Dxx(ψc(t, x))],
            [ψc(t, 0) ~ 0, Dx(ψc(t, 1)) + im * ψc(t, 1) ~ 0]
        )
        expected = [
            Reψc(t, 0) ~ 0,
            Imψc(t, 0) ~ 0,
            Dx(Reψc(t, 1)) ~ Imψc(t, 1),
            Dx(Imψc(t, 1)) ~ -Reψc(t, 1),
        ]
        @test same_equations(PDEBase.get_bcs(system), expected)
    end

    @testset "Complex BC reduced to false is rejected" begin
        @test_throws ArgumentError split(
            [Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [ψ(t, 0) ~ im, ψ(t, 1) ~ 0]
        )
    end

    @testset "Nested real boundary conditions stay independent" begin
        @variables u(..)
        system = PDESystem(
            [Dt(u(t, x)) ~ Dxx(u(t, x))],
            [u(0, x) ~ sin(x), [u(t, 0) ~ 0, u(t, 1) ~ 0]],
            domain, [t, x], [u(t, x)]; name = :nested_real_bc_test
        )
        normalized, complexmap = PDEBase.handle_complex(system)
        @test complexmap === nothing
        @test isequal(PDEBase.get_dvs(normalized), [u(t, x)])
        @test same_equations(PDEBase.get_bcs(normalized), [
            u(0, x) ~ sin(x), u(t, 0) ~ 0, u(t, 1) ~ 0,
        ])
    end

    @testset "Grouped two-field boundary conditions stay independent" begin
        @variables χ(..) Reχ(..) Imχ(..)
        system = PDESystem(
            [im * Dt(ψ(t, x)) ~ Dxx(χ(t, x)), im * Dt(χ(t, x)) ~ Dxx(ψ(t, x))],
            [[ψ(t, 0) ~ 0, χ(t, 0) ~ 0], [ψ(t, 1) ~ 0, χ(t, 1) ~ 0]],
            domain, [t, x], [ψ(t, x), χ(t, x)]; name = :grouped_two_field_bc_test
        )
        normalized, _ = PDEBase.handle_complex(system)
        expected = [
            Reψ(t, 0) ~ 0, Imψ(t, 0) ~ 0, Reχ(t, 0) ~ 0, Imχ(t, 0) ~ 0,
            Reψ(t, 1) ~ 0, Imψ(t, 1) ~ 0, Reχ(t, 1) ~ 0, Imχ(t, 1) ~ 0,
        ]
        @test same_equations(PDEBase.get_bcs(normalized), expected)
    end

    @testset "Real nonlinear equations split after a complex BC" begin
        system = split_typed(
            [Dt(ψc(t, x)) ~ Dxx(ψc(t, x)) + ψc(t, x)^2],
            [ψc(t, 0) ~ 0, Dx(ψc(t, 1)) + im * ψc(t, 1) ~ 0]
        )
        expected = [
            Dt(Reψc(t, x)) ~ Dxx(Reψc(t, x)) + Reψc(t, x)^2 - Imψc(t, x)^2,
            Dt(Imψc(t, x)) ~ Dxx(Imψc(t, x)) + 2 * Reψc(t, x) * Imψc(t, x),
        ]
        @test same_equations(ModelingToolkit.equations(system), expected)
    end

    @testset "Non-holomorphic BC reconstruction is rejected" begin
        @test_throws ArgumentError split(
            [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [ψ(t, 0) ~ 0, Dx(ψ(t, 1)) ~ im * conj(ψ(t, 1))]
        )
    end

    @testset "Truncated complex literal on the left is rejected" begin
        @test_throws ArgumentError split(
            [Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [im ~ ψ(t, 0), ψ(t, 1) ~ 0]
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
        system = split_typed(
            [im * Dt(ψc(t, x)) ~ Dxx(ψc(t, x))],
            [ψc(0, x) => exp(im * x), ψc(t, 0) ~ exp(im * t), ψc(t, 1) ~ (1 + im) * t]
        )
        expected = [
            Reψc(0, x) ~ exp0 * cos(x),
            Imψc(0, x) ~ exp0 * sin(x),
            Reψc(t, 0) ~ exp0 * cos(t),
            Imψc(t, 0) ~ exp0 * sin(t),
            Reψc(t, 1) ~ t,
            Imψc(t, 1) ~ t,
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
