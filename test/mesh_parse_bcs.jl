using PDEBase
using PDEBase: PieceBoundary, Essential, Natural, Robin
using ModelingToolkit
using Symbolics: unwrap
using SciMLBase
using Test
import DomainSets
using DomainSets: Interval

# A backend mesh implementing the mesh-domain interface, as in test/mesh_domain.jl.
struct PipelineMesh end
PDEBase.meshdim(::PipelineMesh) = 2
PDEBase.facetset_names(::PipelineMesh) = (:inlet, :wall, :outlet)
PDEBase.cellset_names(::PipelineMesh) = (:all,)
PDEBase.contains_point(::PipelineMesh, p) = all(c -> 0 <= c <= 1, p)

struct DummyMeshDiscretization <: PDEBase.AbstractMeshDiscretization
    time::Any
    captured::Base.RefValue{Any}
end
PDEBase.get_time(disc::DummyMeshDiscretization) = disc.time

struct StopAtBoundaryMap <: Exception end

# Stores the parsed boundary map and stops the pipeline there.
function PDEBase.check_boundarymap(bmap, v::VariableMap, disc::DummyMeshDiscretization)
    disc.captured[] = (bmap, v)
    throw(StopAtBoundaryMap())
end

@parameters t x y k α
@variables u(..)
Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
Ω = MeshDomain(PipelineMesh())
Γ_in, Γ_wall, Γ_out = (boundarypiece(Ω, n) for n in (:inlet, :wall, :outlet))
Dn_wall, Dn_out = NormalDerivative(Γ_wall), NormalDerivative(Γ_out)
heat = Dt(u(t, x, y)) ~ k * (Dx(Dx(u(t, x, y))) + Dy(Dy(u(t, x, y))))
domains = [(x, y) ∈ Ω, t ∈ Interval(0.0, 1.0)]

function run_to_boundarymap(bcs)
    @named pdesys = PDESystem([heat], bcs, domains, [t, x, y], [u(t, x, y)], [k, α])
    disc = DummyMeshDiscretization(t, Ref{Any}(nothing))
    @test_throws StopAtBoundaryMap SciMLBase.symbolic_discretize(pdesys, disc)
    return disc.captured[]..., disc
end

@testset "symbolic_discretize reaches parse_bcs for a mesh system" begin
    bcs = [
        u(t, Γ_in) ~ sin(t),
        Dn_wall(u(t, Γ_wall)) ~ 0,
        k * Dn_out(u(t, Γ_out)) + α * u(t, Γ_out) ~ α * 0.3,
        u(0, x, y) ~ 0,
    ]
    bmap, v, disc = run_to_boundarymap(bcs)
    uop = operation(unwrap(u(t, x, y)))
    @test isequal(collect(keys(bmap)), [uop])
    @test issetequal(keys(bmap[uop]), Any[unwrap(t), Γ_in, Γ_wall, Γ_out])

    for (Γ, kind, order) in ((Γ_in, Essential(), 0), (Γ_wall, Natural(), 1), (Γ_out, Robin(), 1))
        b = only(bmap[uop][Γ])
        @test b isa PieceBoundary
        @test b.piece == Γ
        @test b.kind == kind
        @test b.order == order
        @test isequal(b.u, unwrap(u(t, Γ)))
        @test isequal(b.depvars, [unwrap(u(t, x, y))])
    end
    ic = only(bmap[uop][unwrap(t)])
    @test ic isa LowerBoundary
    @test isequal(ic.u, unwrap(u(0, x, y)))
    @test isequal(ic.x, unwrap(t))

    ics, bcmap = PDEBase.split_ics_bcs(bmap, v, disc)
    @test isequal(ics, [ic])
    @test issetequal(keys(bcmap[uop]), [Γ_in, Γ_wall, Γ_out])
    @test bcmap[uop][Γ_wall] == bmap[uop][Γ_wall]
end

@testset "mesh traces: time derivatives and ODE variables" begin
    @variables w(..)
    eqs = [heat, Dt(w(t)) ~ -w(t)]
    bcs = [u(0, x, y) ~ 0, w(0) ~ 1, Dt(u(t, Γ_in)) ~ 1, u(t, Γ_out) ~ w(t)]
    @named pdesys = PDESystem(eqs, bcs, domains, [t, x, y], [u(t, x, y), w(t)], [k, α])
    disc = DummyMeshDiscretization(t, Ref{Any}(nothing))
    @test_throws StopAtBoundaryMap SciMLBase.symbolic_discretize(pdesys, disc)
    bmap, v = disc.captured[]
    uop = operation(unwrap(u(t, x, y)))
    inlet = only(bmap[uop][Γ_in])
    @test inlet.kind == Essential()
    @test inlet.order == 0
    outlet = only(bmap[uop][Γ_out])
    @test outlet.kind == Essential()
    @test issetequal(outlet.depvars, [unwrap(u(t, x, y)), unwrap(w(t))])
end

@testset "mesh parse_bcs errors" begin
    @named pdesys = PDESystem(
        [heat], [u(0, x, y) ~ 0, Dn_wall(u(t, Γ_out)) ~ 0], domains, [t, x, y], [u(t, x, y)], [k, α]
    )
    disc = DummyMeshDiscretization(t, Ref{Any}(nothing))
    err = try
        SciMLBase.symbolic_discretize(pdesys, disc)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("write the trace on the operator's piece", err.msg)
    @test disc.captured[] === nothing

    # A condition on a mesh must be a trace; `u(t, 0, y)` has no interval to be an endpoint of.
    @named endpoint = PDESystem(
        [heat], [u(0, x, y) ~ 0, u(t, 0, y) ~ 0], domains, [t, x, y], [u(t, x, y)], [k, α]
    )
    @test_throws ArgumentError SciMLBase.symbolic_discretize(endpoint, DummyMeshDiscretization(t, Ref{Any}(nothing)))
end
