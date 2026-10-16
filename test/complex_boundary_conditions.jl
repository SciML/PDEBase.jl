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

    @testset "Untyped pre-split complex Dirichlet condition splits exactly" begin
        # `ψ(t, 0) ~ cos(t) + im*sin(t)` on an untyped field reaches PDEBase as
        # the pre-split pair `[ψ(t, 0) ~ cos(t), 0 ~ sin(t)]`.
        expected = [
            Reψ(t, 0) ~ cos(t),
            Imψ(t, 0) ~ sin(t),
            Reψ(t, 1) ~ 0,
            Imψ(t, 1) ~ 0,
        ]
        system = split(
            [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [ψ(t, 0) ~ cos(t) + im * sin(t), ψ(t, 1) ~ 0]
        )
        @test same_equations(PDEBase.get_bcs(system), expected)
        system = split(
            [Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [ψ(t, 0) ~ cos(t) + im * sin(t), ψ(t, 1) ~ 0]
        )
        @test same_equations(PDEBase.get_bcs(system), expected)
    end

    @testset "Untyped exp(im*…) Dirichlet conditions with a complex IC split exactly" begin
        # The MethodOfLines "Schroedinger with complex bcs" shape: a complex
        # initial-condition Pair plus untyped `exp(im …)` Dirichlet data.
        exp0 = SymbolicUtils.term(exp, 0.0)
        system = split(
            [(im * 1.0e-2) * Dt(ψ(t, x)) ~ (-0.5 * 1.0e-4) * Dxx(ψ(t, x))],
            [
                ψ(0, x) => exp((im / 1.0e-2) * 0.1 * x),
                ψ(t, 0) ~ exp((im / 1.0e-2) * (0.0 - 0.5e-2 * t)),
                ψ(t, 1) ~ exp((im / 1.0e-2) * (0.1 - 0.5e-2 * t)),
            ]
        )
        expected = [
            Reψ(0, x) ~ exp0 * cos(10.0x),
            Imψ(0, x) ~ exp0 * sin(10.0x),
            Reψ(t, 0) ~ exp0 * cos(-0.5t),
            Imψ(t, 0) ~ exp0 * sin(-0.5t),
            Reψ(t, 1) ~ exp0 * cos(10.0 - 0.5t),
            Imψ(t, 1) ~ exp0 * sin(10.0 - 0.5t),
        ]
        @test same_equations(PDEBase.get_bcs(system), expected)
    end

    @testset "Untyped residual-form complex Dirichlet condition splits exactly" begin
        # `0 ~ ψ(t, 0) - exp(im*t)` on an untyped field reaches PDEBase as
        # `[0 ~ ψ(t, 0) - cos(t)*exp(0), 0 ~ -sin(t)*exp(0)]`: exactly one
        # member is free of dependent variables, so the pair is rebuilt as
        # `0 ~ ψ - exp(im*t)` and split with the real and imaginary parts
        # coupled. Symbolics keeps the `exp(0)` factor unevaluated.
        exp0 = SymbolicUtils.term(exp, 0)
        system = split(
            [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [0 ~ ψ(t, 0) - exp(im * t), ψ(t, 1) ~ 0]
        )
        expected = [
            0 ~ Reψ(t, 0) - exp0 * cos(t),
            0 ~ Imψ(t, 0) - exp0 * sin(t),
            Reψ(t, 1) ~ 0,
            Imψ(t, 1) ~ 0,
        ]
        @test same_equations(PDEBase.get_bcs(system), expected)

        # A `[0 ~ a, 0 ~ b]` pair with neither member mentioning a dependent
        # variable is neither a real boundary grouping nor a complex
        # condition, so it is rejected.
        @test_throws ArgumentError split(
            [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [[0 ~ sin(t), 0 ~ cos(t)], ψ(t, 1) ~ 0]
        )
    end

    @testset "Untyped complex conditions with the field on the right split exactly" begin
        # `exp(im*t) ~ ψ(t, 0)` pre-splits to `[cos(t) ~ ψ(t, 0), sin(t) ~ 0]`,
        # the mirror image of the field-on-the-left pair.
        exp0 = SymbolicUtils.term(exp, 0)
        equations = [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))]
        system = split(equations, [exp(im * t) ~ ψ(t, 0), ψ(t, 1) ~ 0])
        @test same_equations(
            PDEBase.get_bcs(system), [
                exp0 * cos(t) ~ Reψ(t, 0),
                exp0 * sin(t) ~ Imψ(t, 0),
                Reψ(t, 1) ~ 0,
                Imψ(t, 1) ~ 0,
            ]
        )

        system = split(equations, [ψ(t, 0) ~ 0, exp(im * t) ~ Dx(ψ(t, 1))])
        @test same_equations(
            PDEBase.get_bcs(system), [
                Reψ(t, 0) ~ 0,
                Imψ(t, 0) ~ 0,
                exp0 * cos(t) ~ Dx(Reψ(t, 1)),
                exp0 * sin(t) ~ Dx(Imψ(t, 1)),
            ]
        )

        system = split(equations, [ψ(t, 0) - exp(im * t) ~ 0, ψ(t, 1) ~ 0])
        @test same_equations(
            PDEBase.get_bcs(system), [
                Reψ(t, 0) - exp0 * cos(t) ~ 0,
                Imψ(t, 0) - exp0 * sin(t) ~ 0,
                Reψ(t, 1) ~ 0,
                Imψ(t, 1) ~ 0,
            ]
        )

        @test_throws ArgumentError split(
            equations, [[sin(t) ~ 0, cos(t) ~ 0], ψ(t, 1) ~ 0]
        )
    end

    @testset "Pre-split pair with dependent-variable data is rejected" begin
        # `Dx(ψ(t, 1)) ~ im*ψ(t, 1)` pre-splits to `[Dx(ψ) ~ 0, 0 ~ ψ]`: the
        # `0 ~` side carries a dependent variable, so no unique reconstruction
        # exists and the pair is rejected.
        @test_throws ArgumentError split(
            [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [ψ(t, 0) ~ 0, Dx(ψ(t, 1)) ~ im * ψ(t, 1)]
        )
        @test_throws ArgumentError split(
            [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [[0 ~ sin(t), ψ(t, 0) ~ 1], ψ(t, 1) ~ 0]
        )
    end

    @testset "Nested pair of `0 ~` boundary conditions stays real" begin
        @variables u(..)
        system = PDESystem(
            [Dt(u(t, x)) ~ Dxx(u(t, x))],
            [u(0, x) ~ 0, [0 ~ u(t, 0) - sin(t), 0 ~ u(t, 1)]],
            domain, [t, x], [u(t, x)]; name = :zero_lhs_real_pair_test
        )
        normalized, complexmap = PDEBase.handle_complex(system)
        @test complexmap === nothing
        @test same_equations(
            PDEBase.get_bcs(normalized),
            [u(0, x) ~ 0, 0 ~ u(t, 0) - sin(t), 0 ~ u(t, 1)]
        )

        broadcast_system = PDESystem(
            [Dt(u(t, x)) ~ Dxx(u(t, x))],
            [u(0, x) ~ sin(x), 0 .~ Dx.(u.(t, [0, 1]))],
            domain, [t, x], [u(t, x)]; name = :zero_lhs_broadcast_test
        )
        normalized, complexmap = PDEBase.handle_complex(broadcast_system)
        @test complexmap === nothing
        @test same_equations(
            PDEBase.get_bcs(normalized),
            [u(0, x) ~ sin(x), 0 ~ Dx(u(t, 0)), 0 ~ Dx(u(t, 1))]
        )
    end

    @testset "Nested real groups with a leading `0 ~` member stay real" begin
        @variables u(..)
        for group in (
                [0 ~ u(t, 0), u(t, 1) ~ 1],
                [0 ~ u(t, 0), 1 ~ u(t, 1)],
            )
            system = PDESystem(
                [Dt(u(t, x)) ~ Dxx(u(t, x))], [u(0, x) ~ 0, group],
                domain, [t, x], [u(t, x)]; name = :leading_zero_real_group_test
            )
            normalized, complexmap = PDEBase.handle_complex(system)
            @test complexmap === nothing
            @test same_equations(PDEBase.get_bcs(normalized), [u(0, x) ~ 0; group])
        end

        # On a complex system the same shape is also the pre-split of
        # `im*ψ(t, 1) ~ ψ(t, 0) + im`, so it is rejected.
        @test_throws ArgumentError split(
            [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [ψ(0, x) ~ 0, [0 ~ ψ(t, 0), ψ(t, 1) ~ 1]]
        )
    end

    @testset "Untyped complex conditions on a real equation are rejected" begin
        real_pde = [Dt(ψ(t, x)) ~ Dxx(ψ(t, x))]
        for bcs in (
                [ψ(t, 0) ~ 0, Dx(ψ(t, 1)) ~ im * ψ(t, 1)],
                [ψ(t, 0) ~ im * ψ(t, 1), ψ(t, 1) ~ 0],
                [ψ(t, 0) ~ 0, Dx(ψ(t, 1)) ~ (1 + im) * ψ(t, 1)],
                [ψ(t, 0) ~ 0, Dx(ψ(t, 1)) ~ im * conj(ψ(t, 1))],
            )
            @test_throws ArgumentError split(real_pde, bcs)
        end
    end

    @testset "Nested pair with a zero left-hand side is rejected" begin
        @variables u(..)
        system = PDESystem(
            [Dt(u(t, x)) ~ Dxx(u(t, x))],
            [u(0, x) ~ 0, [u(t, 0) ~ 1, 0 ~ Dx(u(t, 0))], u(t, 1) ~ 0],
            domain, [t, x], [u(t, x)]; name = :zero_lhs_pair_test
        )
        @test_throws ArgumentError PDEBase.handle_complex(system)
    end

    @testset "Typed complex Dirichlet conditions split exactly" begin
        system = split_typed(
            [im * Dt(ψc(t, x)) ~ Dxx(ψc(t, x))],
            [ψc(t, 0) ~ cos(t) + im * sin(t), ψc(t, 1) ~ 0]
        )
        expected = [
            Reψc(t, 0) ~ cos(t),
            Imψc(t, 0) ~ sin(t),
            Reψc(t, 1) ~ 0,
            Imψc(t, 1) ~ 0,
        ]
        @test same_equations(PDEBase.get_bcs(system), expected)
    end

    @testset "Typed imaginary literal boundary condition splits exactly" begin
        system = split_typed(
            [Dt(ψc(t, x)) ~ Dxx(ψc(t, x))],
            [ψc(t, 0) ~ im, ψc(t, 1) ~ 0]
        )
        expected = [
            Reψc(t, 0) ~ 0,
            Imψc(t, 0) ~ 1,
            Reψc(t, 1) ~ 0,
            Imψc(t, 1) ~ 0,
        ]
        @test same_equations(PDEBase.get_bcs(system), expected)
    end

    if isdefined(Symbolics, :split_complex_equation)
        @testset "Marked complex boundary conditions preserve coupling" begin
            marked = Symbolics.split_complex_equation(
                ψ(t, 0), cos(t) + im * sin(t)
            )
            @test marked isa Symbolics.SplitComplexEquation

            system = split(
                [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
                [marked, ψ(t, 1) ~ 0]
            )
            expected = [
                Reψ(t, 0) ~ cos(t),
                Imψ(t, 0) ~ sin(t),
                Reψ(t, 1) ~ 0,
                Imψ(t, 1) ~ 0,
            ]
            @test same_equations(PDEBase.get_bcs(system), expected)

            marked_neumann = Symbolics.split_complex_equation(
                Dx(ψ(t, 1)), im * ψ(t, 1)
            )
            @test marked_neumann isa Symbolics.SplitComplexEquation
            system = split(
                [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
                [ψ(t, 0) ~ 0, marked_neumann]
            )
            expected = [
                Reψ(t, 0) ~ 0,
                Imψ(t, 0) ~ 0,
                Dx(Reψ(t, 1)) ~ -Imψ(t, 1),
                Dx(Imψ(t, 1)) ~ Reψ(t, 1),
            ]
            @test same_equations(PDEBase.get_bcs(system), expected)

            system = split(
                [Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
                [marked, ψ(t, 1) ~ 0]
            )
            @test same_equations(
                PDEBase.get_bcs(system), [
                    Reψ(t, 0) ~ cos(t),
                    Imψ(t, 0) ~ sin(t),
                    Reψ(t, 1) ~ 0,
                    Imψ(t, 1) ~ 0,
                ]
            )

            marked_typed = Symbolics.split_complex_equation(
                ψc(t, 0), cos(t) + im * sin(t)
            )
            @test marked_typed isa Symbolics.SplitComplexEquation
            system = split_typed(
                [im * Dt(ψc(t, x)) ~ Dxx(ψc(t, x))],
                [marked_typed, ψc(t, 1) ~ 0]
            )
            @test same_equations(
                PDEBase.get_bcs(system), [
                    Reψc(t, 0) ~ cos(t),
                    Imψc(t, 0) ~ sin(t),
                    Reψc(t, 1) ~ 0,
                    Imψc(t, 1) ~ 0,
                ]
            )
        end
    end

    @testset "Hand-written real/imag pair reads as its complex equivalent" begin
        # On an untyped field `real(ψ)` and `imag(ψ)` simplify to `ψ` and `0`,
        # so a hand-written `real`/`imag` pair arrives in the same shape as a
        # Symbolics pre-split of `ψ(t, 0) ~ cos(t) + im*sin(t)` and splits the
        # same way.
        manual_pair = [
            real(ψ(t, 0)) ~ cos(t),
            imag(ψ(t, 0)) ~ sin(t),
        ]
        system = split(
            [im * Dt(ψ(t, x)) ~ Dxx(ψ(t, x))],
            [manual_pair, ψ(t, 1) ~ 0]
        )
        @test same_equations(
            PDEBase.get_bcs(system),
            [Reψ(t, 0) ~ cos(t), Imψ(t, 0) ~ sin(t), Reψ(t, 1) ~ 0, Imψ(t, 1) ~ 0]
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
        @test same_equations(
            PDEBase.get_bcs(normalized), [
                u(0, x) ~ sin(x), u(t, 0) ~ 0, u(t, 1) ~ 0,
            ]
        )
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

    @testset "Grouped clamped beam boundary conditions stay real" begin
        @variables u(..)
        beam = [Dt(Dt(u(t, x))) ~ -Dx(Dx(Dx(Dx(u(t, x)))))]
        system = PDESystem(
            beam,
            [
                u(0, x) ~ sin(x), Dt(u(0, x)) ~ 0,
                [u(t, 0) ~ 0, Dx(u(t, 0)) ~ 0],
                [u(t, 1) ~ 0, Dx(u(t, 1)) ~ 0],
            ],
            domain, [t, x], [u(t, x)]; name = :grouped_beam_bc_test
        )
        normalized, complexmap = PDEBase.handle_complex(system)
        @test complexmap === nothing
        @test isequal(PDEBase.get_dvs(normalized), [u(t, x)])
        @test same_equations(
            PDEBase.get_bcs(normalized), [
                u(0, x) ~ sin(x), Dt(u(0, x)) ~ 0,
                u(t, 0) ~ 0, Dx(u(t, 0)) ~ 0,
                u(t, 1) ~ 0, Dx(u(t, 1)) ~ 0,
            ]
        )
    end

    @testset "Grouped interface boundary conditions stay real" begin
        @variables u(..) v(..)
        system = PDESystem(
            [Dt(u(t, x)) ~ Dxx(u(t, x)), Dt(v(t, x)) ~ Dxx(v(t, x))],
            [
                u(0, x) ~ 0, v(0, x) ~ 0,
                [u(t, 0) ~ v(t, 0), Dx(u(t, 0)) ~ Dx(v(t, 0))],
                u(t, 1) ~ 1, v(t, 1) ~ 0,
            ],
            domain, [t, x], [u(t, x), v(t, x)]; name = :grouped_interface_bc_test
        )
        normalized, complexmap = PDEBase.handle_complex(system)
        @test complexmap === nothing
        @test isequal(PDEBase.get_dvs(normalized), [u(t, x), v(t, x)])
        @test same_equations(
            PDEBase.get_bcs(normalized), [
                u(0, x) ~ 0, v(0, x) ~ 0,
                u(t, 0) ~ v(t, 0), Dx(u(t, 0)) ~ Dx(v(t, 0)),
                u(t, 1) ~ 1, v(t, 1) ~ 0,
            ]
        )
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

    @testset "Typed non-holomorphic boundary conditions split exactly" begin
        # conj/real/imag of an untyped field simplify before any marker is
        # built, so non-holomorphic conditions require a `::Complex` field.
        system = split_typed(
            [im * Dt(ψc(t, x)) ~ Dxx(ψc(t, x))],
            [ψc(t, 0) ~ 0, Dx(ψc(t, 1)) ~ im * conj(ψc(t, 1))]
        )
        expected = [
            Reψc(t, 0) ~ 0,
            Imψc(t, 0) ~ 0,
            Dx(Reψc(t, 1)) ~ Imψc(t, 1),
            Dx(Imψc(t, 1)) ~ Reψc(t, 1),
        ]
        @test same_equations(PDEBase.get_bcs(system), expected)
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
