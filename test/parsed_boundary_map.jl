using PDEBase
using PDEBase: parsed_boundary_map
using PDEBase: PieceBoundary, LowerBoundary, Essential, Natural, Robin
using ModelingToolkit
using Symbolics: unwrap
using SciMLBase
using Test
import DomainSets
using DomainSets: Interval

struct AccessorMesh end
PDEBase.meshdim(::AccessorMesh) = 2
PDEBase.facetset_names(::AccessorMesh) = (:inlet, :wall, :outlet)
PDEBase.cellset_names(::AccessorMesh) = (:all,)
PDEBase.contains_point(::AccessorMesh, p) = all(c -> 0 <= c <= 1, p)

struct AccessorMeshDiscretization <: PDEBase.AbstractMeshDiscretization
    time::Any
end
struct AccessorBoxDiscretization <: PDEBase.AbstractBoxDiscretization
    time::Any
end
PDEBase.get_time(disc::Union{AccessorMeshDiscretization, AccessorBoxDiscretization}) = disc.time

struct DiscretizationReached <: Exception end

# Stops `symbolic_discretize` right after the boundary map is checked, keeping the map.
const CAPTURED = Ref{Any}(nothing)
function PDEBase.check_boundarymap(bmap, v::VariableMap, disc::AccessorBoxDiscretization)
    CAPTURED[] = bmap
    return nothing
end
PDEBase.construct_disc_state(::AccessorBoxDiscretization) = throw(DiscretizationReached())

function same_map(a, b)
    isequal(collect(keys(a)), collect(keys(b))) || return false
    for u in keys(a)
        issetequal(keys(a[u]), keys(b[u])) || return false
        for x in keys(a[u])
            length(a[u][x]) == length(b[u][x]) || return false
            for (p, q) in zip(a[u][x], b[u][x])
                typeof(p) === typeof(q) || return false
                all(f -> isequal(getfield(p, f), getfield(q, f)), fieldnames(typeof(p))) || return false
            end
        end
    end
    return true
end

@parameters t x y k α
@variables u(..)
Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
heat = Dt(u(t, x, y)) ~ k * (Dx(Dx(u(t, x, y))) + Dy(Dy(u(t, x, y))))

@testset "parsed_boundary_map on a mesh domain, without discretizing" begin
    Ω = MeshDomain(AccessorMesh())
    Γ_in, Γ_wall, Γ_out = (boundarypiece(Ω, n) for n in (:inlet, :wall, :outlet))
    Dn_wall, Dn_out = NormalDerivative(Γ_wall), NormalDerivative(Γ_out)
    bcs = [
        u(t, Γ_in) ~ sin(t),
        Dn_wall(u(t, Γ_wall)) ~ 0,
        k * Dn_out(u(t, Γ_out)) + α * u(t, Γ_out) ~ α * 0.3,
        u(0, x, y) ~ 0,
    ]
    domains = [(x, y) ∈ Ω, t ∈ Interval(0.0, 1.0)]
    @named pdesys = PDESystem([heat], bcs, domains, [t, x, y], [u(t, x, y)], [k, α])
    eqs_before = copy(PDEBase.get_eqs(pdesys))

    parsed = parsed_boundary_map(pdesys, AccessorMeshDiscretization(t))
    (; boundarymap, v) = parsed
    lowered = parsed.pdesys
    uop = operation(unwrap(u(t, x, y)))
    @test issetequal(keys(boundarymap[uop]), Any[unwrap(t), Γ_in, Γ_wall, Γ_out])
    for (Γ, kind, order) in ((Γ_in, Essential(), 0), (Γ_wall, Natural(), 1), (Γ_out, Robin(), 1))
        b = only(boundarymap[uop][Γ])
        @test b isa PieceBoundary
        @test b.piece == Γ
        @test b.kind == kind
        @test b.order == order
    end
    @test only(boundarymap[uop][unwrap(t)]) isa LowerBoundary
    @test v isa VariableMap
    @test lowered isa PDESystem
    @test all(isequal.(PDEBase.get_eqs(pdesys), eqs_before))
end

@testset "parsed_boundary_map agrees with symbolic_discretize on a box" begin
    domains = [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0), y ∈ Interval(0.0, 1.0)]
    bcs = [
        u(0, x, y) ~ 0, u(t, 0, y) ~ 0, u(t, 1, y) ~ 1,
        Dy(u(t, x, 0)) ~ 0, u(t, x, 1) ~ sin(x),
    ]
    @named pdesys = PDESystem([heat], bcs, domains, [t, x, y], [u(t, x, y)], [k])
    disc = AccessorBoxDiscretization(t)

    (; boundarymap) = parsed_boundary_map(pdesys, disc)
    @test_throws DiscretizationReached SciMLBase.symbolic_discretize(pdesys, disc)
    @test CAPTURED[] !== nothing
    @test same_map(boundarymap, CAPTURED[])
end
