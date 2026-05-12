using JET
using PDEBase
using Test
import Symbolics
using SymbolicUtils: operation

@testset "JET static analysis" begin
    # Test core symbolic utility functions for type stability
    # JET reports from dependencies (Symbolics, SymbolicUtils) are filtered
    # via target_modules to avoid noise from transitive symbolic-stack issues.
    @eval Symbolics.@variables x y t
    @eval Symbolics.@variables u(..)

    x_sym = Symbolics.unwrap(x)
    u_term = Symbolics.unwrap(u(x, y, t))
    depvar_ops = [operation(u_term)]

    Dx = Symbolics.Differential(x)
    term = Symbolics.unwrap(Dx(u(x, y, t)))

    @testset "safe_unwrap" begin
        rep = JET.report_call(
            PDEBase.safe_unwrap, (typeof(u(x, y, t)),);
            target_modules = (PDEBase,)
        )
        @test length(JET.get_reports(rep)) == 0
    end

    @testset "has_derivatives" begin
        rep = JET.report_call(
            PDEBase.has_derivatives, (typeof(term),);
            target_modules = (PDEBase,)
        )
        @test length(JET.get_reports(rep)) == 0
    end

    @testset "count_differentials" begin
        rep = JET.report_call(
            PDEBase.count_differentials, (typeof(term), typeof(x_sym));
            target_modules = (PDEBase,)
        )
        @test length(JET.get_reports(rep)) == 0
    end

    @testset "find_derivative" begin
        rep = JET.report_call(
            PDEBase.find_derivative, (typeof(term), typeof(operation(u_term)));
            target_modules = (PDEBase,)
        )
        @test length(JET.get_reports(rep)) == 0
    end

    @testset "recursive_unwrap" begin
        rep = JET.report_call(
            PDEBase.recursive_unwrap, (typeof(term),);
            target_modules = (PDEBase,)
        )
        @test length(JET.get_reports(rep)) == 0
    end
end
