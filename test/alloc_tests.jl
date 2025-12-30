using PDEBase
using Symbolics
using SymbolicUtils: operation
using BenchmarkTools
using Test

@testset "Allocation Tests" begin
    # Setup test data
    @variables x y t
    @variables u(..)

    x_sym = Symbolics.unwrap(x)
    y_sym = Symbolics.unwrap(y)
    t_sym = Symbolics.unwrap(t)
    u_term = Symbolics.unwrap(u(x, y, t))
    depvar_ops = [operation(u_term)]

    Dx = Differential(x)
    Dy = Differential(y)
    term = Symbolics.unwrap(Dx(u(x, y, t)))
    term2 = Symbolics.unwrap(Dx(Dy(u(x, y, t))))
    complex_term = Symbolics.unwrap(Dx(u(x, y, t)) + Dy(u(x, y, t)) + u(x, y, t))
    complex_eq = Symbolics.unwrap(Dx(u(x, y, t)) ~ u(x, y, t) + t * y)

    @testset "Performance regression check" begin
        # These benchmarks ensure we don't regress
        # Using median time which is stable across runs
        # Times allow margin for CI variability

        b1 = @benchmark PDEBase.count_differentials($term, $x_sym) samples = 10 evals = 100
        @test median(b1).time < 50_000  # Under 50 microseconds

        b2 = @benchmark PDEBase.differential_order($term, $x_sym) samples = 10 evals = 100
        @test median(b2).time < 50_000  # Under 50 microseconds

        b3 = @benchmark PDEBase.has_derivatives($term) samples = 10 evals = 100
        @test median(b3).time < 1_000  # Under 1 microsecond

        b4 = @benchmark PDEBase.get_depvars($term, $depvar_ops) samples = 10 evals = 100
        @test median(b4).time < 10_000  # Under 10 microseconds

        b5 = @benchmark PDEBase.find_derivative(
            $term, $(operation(u_term))) samples = 10 evals = 100
        @test median(b5).time < 1_000  # Under 1 microsecond

        b6 = @benchmark PDEBase.split_terms($complex_eq) samples = 10 evals = 100
        @test median(b6).time < 50_000  # Under 50 microseconds

        b7 = @benchmark PDEBase.split_additive_terms($complex_eq) samples = 10 evals = 100
        @test median(b7).time < 10_000  # Under 10 microseconds
    end

    @testset "Allocation bounds" begin
        # These tests verify allocations stay bounded
        # We use benchmark minimum allocations for stability

        b1 = @benchmark PDEBase.count_differentials($term, $x_sym) samples = 5 evals = 50
        @test minimum(b1).allocs <= 15  # Bounded allocations

        b2 = @benchmark PDEBase.differential_order($term, $x_sym) samples = 5 evals = 50
        @test minimum(b2).allocs <= 30  # Bounded allocations

        b3 = @benchmark PDEBase.has_derivatives($term) samples = 5 evals = 50
        @test minimum(b3).allocs == 0  # Early return, no allocations

        b4 = @benchmark PDEBase.get_depvars($term, $depvar_ops) samples = 5 evals = 50
        @test minimum(b4).allocs <= 15  # Bounded allocations

        b5 = @benchmark PDEBase.find_derivative(
            $term, $(operation(u_term))) samples = 5 evals = 50
        @test minimum(b5).allocs == 0  # Should be zero for direct match
    end
end
