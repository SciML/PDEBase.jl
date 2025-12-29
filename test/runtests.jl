using SafeTestsets, Test

const GROUP = get(ENV, "GROUP", "All")

const is_APPVEYOR = Sys.iswindows() && haskey(ENV, "APPVEYOR")

const is_TRAVIS = haskey(ENV, "TRAVIS")

const is_CI = haskey(ENV, "CI")

# Currently verified by Downstream tests
@test true

@safetestset "JET static analysis" begin
    using JET
    using PDEBase
    import Symbolics
    using SymbolicUtils: operation

    # Test core symbolic utility functions for type stability
    # Note: JET reports from dependencies (Symbolics, SymbolicUtils) are filtered
    @eval Symbolics.@variables x y t
    @eval Symbolics.@variables u(..)

    x_sym = Symbolics.unwrap(x)
    u_term = Symbolics.unwrap(u(x, y, t))
    depvar_ops = [operation(u_term)]

    Dx = Symbolics.Differential(x)
    term = Symbolics.unwrap(Dx(u(x, y, t)))

    # Test safe_unwrap - simple function that should be type stable
    @testset "safe_unwrap" begin
        rep = JET.report_call(PDEBase.safe_unwrap, (typeof(u(x, y, t)),);
                              target_modules=(PDEBase,))
        @test length(JET.get_reports(rep)) == 0
    end

    # Test has_derivatives with target_modules filter
    @testset "has_derivatives" begin
        rep = JET.report_call(PDEBase.has_derivatives, (typeof(term),);
                              target_modules=(PDEBase,))
        @test length(JET.get_reports(rep)) == 0
    end

    # Test count_differentials with target_modules filter
    @testset "count_differentials" begin
        rep = JET.report_call(PDEBase.count_differentials, (typeof(term), typeof(x_sym));
                              target_modules=(PDEBase,))
        @test length(JET.get_reports(rep)) == 0
    end

    # Test find_derivative with target_modules filter
    @testset "find_derivative" begin
        rep = JET.report_call(PDEBase.find_derivative, (typeof(term), typeof(operation(u_term)));
                              target_modules=(PDEBase,))
        @test length(JET.get_reports(rep)) == 0
    end

    # Test recursive_unwrap with target_modules filter
    @testset "recursive_unwrap" begin
        rep = JET.report_call(PDEBase.recursive_unwrap, (typeof(term),);
                              target_modules=(PDEBase,))
        @test length(JET.get_reports(rep)) == 0
    end
end
