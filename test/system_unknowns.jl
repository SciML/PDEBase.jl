using PDEBase
using ModelingToolkit
using ModelingToolkitBase: complete, unknowns
using SciMLBase
using Symbolics
using Test

struct UnknownsDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
end
struct ScalarUnknownsSpace <: PDEBase.AbstractDiscreteSpace
    discvars::Any
end
struct ArrayUnknownsSpace <: PDEBase.AbstractDiscreteSpace
    discvars::Any
    arrays::Any
end
struct UnknownsMetadata <: SciMLBase.AbstractDiscretizationMetadata{true}
    pdesys::Any
end

PDEBase.get_time(disc::UnknownsDiscretization) = disc.time
PDEBase.get_discvars(s::ScalarUnknownsSpace) = s.discvars
PDEBase.get_discvars(s::ArrayUnknownsSpace) = s.discvars
PDEBase.get_system_unknowns(s::ArrayUnknownsSpace) = s.arrays

@testset "get_system_unknowns chooses the unknowns of the generated system" begin
    @parameters t x
    @variables u(..)
    Dt = Differential(t)
    Dx = Differential(x)
    pdesys = PDESystem(
        [Dt(u(t, x)) ~ Dx(Dx(u(t, x)))],
        [u(0, x) ~ 0, u(t, 0) ~ 0, u(t, 1) ~ 0],
        [t ∈ (0, 1), x ∈ (0, 1)],
        [t, x],
        [u(t, x)];
        name = :pdebase_system_unknowns_test
    )

    @variables w(t)[1:3]
    elements = collect(w)
    state = PDEBase.EquationState(
        [Dt(w[1]) ~ -w[1], Dt(w[2]) ~ w[1] - w[2], Dt(w[3]) ~ w[2] - w[3]], Equation[]
    )
    u0 = [w[i] => Float64(i) for i in 1:3]
    disc = UnknownsDiscretization(t)
    metadata = UnknownsMetadata(pdesys)

    scalar_space = ScalarUnknownsSpace(Dict(:w => elements))
    @test isequal(PDEBase.get_system_unknowns(scalar_space), vec(elements))
    scalar_sys, _ = PDEBase.generate_system(
        state, scalar_space, u0, (0.0, 1.0), metadata, disc
    )
    @test length(unknowns(scalar_sys)) == 3

    array_space = ArrayUnknownsSpace(Dict(:w => elements), [w])
    @test isequal(PDEBase.get_system_unknowns(array_space), [w])
    array_sys, _ = PDEBase.generate_system(
        state, array_space, u0, (0.0, 1.0), metadata, disc
    )
    @test length(unknowns(array_sys)) == 1
    @test isequal(only(unknowns(array_sys)), unwrap(w))

    # The state vector is the same whichever representation the unknowns take.
    scalar_prob = ODEProblem(complete(scalar_sys), u0, (0.0, 1.0))
    array_prob = ODEProblem(complete(array_sys), u0, (0.0, 1.0))
    @test scalar_prob.u0 == array_prob.u0 == [1.0, 2.0, 3.0]
    @test scalar_prob.f(scalar_prob.u0, scalar_prob.p, 0.0) ==
        array_prob.f(array_prob.u0, array_prob.p, 0.0) == [-1.0, -1.0, -1.0]
end
