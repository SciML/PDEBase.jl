using PDEBase
using SciMLBase
using ModelingToolkit
using Test

# A minimal optimization-system discretization: every equation and boundary condition is
# recorded with its kind, and `generate_system` returns the recorded residuals.
struct ResidualDiscretization <: PDEBase.AbstractOptimizationSystemDiscretization end
struct ResidualSpace <: PDEBase.AbstractDiscreteSpace
    varmap::PDEBase.VariableMap
    nbcs::Int
end

mutable struct ResidualState <: PDEBase.AbstractDiscretizationState
    residuals::Vector{Any}
end
struct ResidualMetadata <: SciMLBase.AbstractDiscretizationMetadata{Val(false)}
    pdesys::Any
    metadata::Base.RefValue{Any}
end

PDEBase.construct_disc_state(::ResidualDiscretization) = ResidualState(Any[])
function PDEBase.construct_discrete_space(
        v::PDEBase.VariableMap, pdesys::PDESystem, ::ResidualDiscretization
    )
    return ResidualSpace(v, length(pdesys.bcs))
end
PDEBase.construct_differential_discretizer(pdesys, s::ResidualSpace, ::ResidualDiscretization, orders) = orders
function PDEBase.discretize_equation!(
        state::ResidualState, eq::Equation, kind::Symbol, s::ResidualSpace, orders,
        ::ResidualDiscretization
    )
    push!(state.residuals, (kind, eq.lhs - eq.rhs, orders))
    return nothing
end
function PDEBase.generate_metadata(
        s::ResidualSpace, ::ResidualDiscretization, pdesys, boundarymap, complexmap, u0
    )
    @test boundarymap === nothing
    @test u0 == []
    return ResidualMetadata(pdesys, Ref{Any}(nothing))
end
function PDEBase.generate_system(
        state::ResidualState, s::ResidualSpace, u0, tspan, metadata::ResidualMetadata,
        ::ResidualDiscretization; checks = true
    )
    @test u0 === nothing && tspan === nothing
    return (; residuals = state.residuals, space = s, metadata)
end

@testset "Optimization-system discretization driver" begin
    @parameters t x
    @variables u(..)
    Dt = Differential(t)
    Dxx = Differential(x)^2
    eq = Dt(u(t, x)) ~ Dxx(u(t, x))
    # The interior condition `u(0.5, x)` is not on a boundary of the domain; the
    # optimization driver must accept it because it does not parse a boundary map.
    bcs = [u(0, x) ~ sinpi(x), u(t, 0) ~ 0, u(t, 1) ~ 0, u(0.5, x) ~ 0]
    pdesys = PDESystem(
        [eq], bcs, [t ∈ (0, 1), x ∈ (0, 1)], [t, x], [u(t, x)]; name = :heat
    )
    result = symbolic_discretize(pdesys, ResidualDiscretization())

    kinds = first.(result.residuals)
    @test kinds == [:pde, :bc, :bc, :bc, :bc]
    @test result.space.nbcs == 4
    @test result.metadata isa ResidualMetadata
    @test PDEBase.get_time(ResidualDiscretization()) === nothing
    orders = result.residuals[1][3]
    @test orders[t] == [1] && orders[x] == [2]
    # `cardinalize_eqs!` normalizes the equations to `lhs - rhs ~ 0`.
    @test isequal(result.residuals[1][2], eq.lhs - eq.rhs)
    # Hooks without an override fall back to the equation-system defaults.
    struct BareDiscretization <: PDEBase.AbstractOptimizationSystemDiscretization end
    v = PDEBase.VariableMap(pdesys, BareDiscretization())
    @test PDEBase.construct_discrete_space(v, pdesys, BareDiscretization()) === nothing
    @test PDEBase.discretize_equation!(
        nothing, eq, :pde, nothing, nothing, BareDiscretization()
    ) === nothing
end

@testset "Optimization-system replaced variable map" begin
    @parameters x
    @variables u(..)[1:2]
    eqs = [u(x)[1] ~ u(x)[2]]
    bcs = [u(0)[1] ~ 0, u(1)[2] ~ 1]
    pdesys = PDESystem(eqs, bcs, [x ∈ (0, 1)], [x], [u(x)]; name = :array_system)
    result = symbolic_discretize(pdesys, ResidualDiscretization())

    replacements = PDEBase.replaced_vars(result.space.varmap)
    @test length(replacements) == 2
    @test Set(values(replacements)) == Set([u(x)[1], u(x)[2]])
end
