using PDEBase
using ModelingToolkit
using Test
import DomainSets
using DomainSets: Interval, Rectangle, (..)

# A backend mesh implementing the mesh-domain interface: the unit square split into
# two named cell sets, with four named facet sets.
struct ToySquareMesh end
PDEBase.meshdim(::ToySquareMesh) = 2
PDEBase.facetset_names(::ToySquareMesh) = (:inlet, :wall, :outlet)
PDEBase.cellset_names(::ToySquareMesh) = (:left_half, :right_half)
PDEBase.contains_point(::ToySquareMesh, p) = all(c -> 0 <= c <= 1, p)

struct TimedDiscretization <: PDEBase.AbstractEquationSystemDiscretization
    time::Any
end
PDEBase.get_time(disc::TimedDiscretization) = disc.time

struct ToyCurveMesh end
PDEBase.meshdim(::ToyCurveMesh) = 3

struct NativeMesh
    points::Vector{Float64}
end
struct WrappedMesh
    native::NativeMesh
end
PDEBase.meshdim(::WrappedMesh) = 1
PDEBase.backend_mesh(m::WrappedMesh) = m.native

@testset "MeshDomain" begin
    Ω = MeshDomain(ToySquareMesh())
    @test Ω isa DomainSets.Domain{DomainSets.SVector{2, Float64}}
    @test DomainSets.dimension(Ω) == 2
    @test [0.5, 0.5] ∈ Ω
    @test [1.5, 0.5] ∉ Ω
    @test DomainSets.SVector(0.0, 1.0) ∈ Ω
    @test PDEBase.backend_mesh(Ω.mesh) === Ω.mesh

    m = WrappedMesh(NativeMesh([0.0, 1.0]))
    @test PDEBase.backend_mesh(MeshDomain(m).mesh) === m.native
    err = try
        [0.5] ∈ MeshDomain(m)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("WrappedMesh", err.msg)
    @test occursin("contains_point", err.msg)

    @testset "Region" begin
        r = Region(Ω, :left_half)
        @test r isa DomainSets.Domain{DomainSets.SVector{2, Float64}}
        @test r.parent === Ω
        @test r.name === :left_half
        @test_throws ArgumentError Region(Ω, :nowhere)
        rs = regions(Ω)
        @test [r.name for r in rs] == [:left_half, :right_half]
        @test all(r -> r.parent === Ω, rs)
    end
end

@testset "VariableMap domains" begin
    @parameters t x y z
    @variables u(..)
    Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
    eq = Dt(u(t, x, y)) ~ Dx(Dx(u(t, x, y))) + Dy(Dy(u(t, x, y)))
    disc = TimedDiscretization(t)

    @testset "tuple pairing on a Rectangle" begin
        rect = Rectangle(0.0 .. 1.0, -1.0 .. 2.0)
        tdomain = t ∈ Interval(0.0, 3.0)
        rectdomain = (x, y) ∈ rect
        bcs = [u(0, x, y) ~ 0, u(t, 0, y) ~ 0, u(t, 1, y) ~ 0, u(t, x, -1) ~ 0, u(t, x, 2) ~ 0]
        @named pdesys = PDESystem([eq], bcs, [rectdomain, tdomain], [t, x, y], [u(t, x, y)])
        v = VariableMap(pdesys, disc)
        @test v.intervals[x] == (0.0, 1.0)
        @test v.intervals[y] == (-1.0, 2.0)
        @test v.intervals[t] == (0.0, 3.0)
        @test v.domains[x] === rectdomain
        @test v.domains[y] === rectdomain
        @test v.domains[t] === tdomain
        @test spatial_domain(v) === rect
    end

    @testset "tuple pairing on a MeshDomain" begin
        Ω = MeshDomain(ToySquareMesh())
        meshdomain = (x, y) ∈ Ω
        @named pdesys = PDESystem(
            [eq], [u(0, x, y) ~ 0], [meshdomain, t ∈ Interval(0.0, 1.0)], [t, x, y], [u(t, x, y)]
        )
        v = VariableMap(pdesys, disc)
        @test isequal(v.x̄, [x, y])
        @test !haskey(v.intervals, x)
        @test !haskey(v.intervals, y)
        @test v.intervals[t] == (0.0, 1.0)
        @test v.domains[x] === meshdomain
        @test v.domains[y] === meshdomain
        @test spatial_domain(v) === Ω

        # A trace on a named subset of the mesh stands in for a boundary piece.
        Γ = Region(Ω, :left_half)
        @named traced = PDESystem(
            [eq], [u(0, x, y) ~ 0, u(t, Γ) ~ sin(t)], [meshdomain, t ∈ Interval(0.0, 1.0)],
            [t, x, y], [u(t, x, y)]
        )
        vt = VariableMap(traced, disc)
        @test isequal(vt.ū, [u(t, x, y)])
        @test isequal(vt.x̄, [x, y])
        @test spatial_domain(vt) === Ω

        @named mismatch = PDESystem(
            [eq], [u(0, x, y) ~ 0], [(x, y) ∈ MeshDomain(ToyCurveMesh()), t ∈ Interval(0.0, 1.0)],
            [t, x, y], [u(t, x, y)]
        )
        err = try
            VariableMap(mismatch, disc)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("has 3 coordinates but the pairing names 2", err.msg)
    end

    @testset "scalar interval pairings" begin
        xdomain, ydomain = x ∈ Interval(0.0, 1.0), y ∈ Interval(2.0, 3.0)
        @named pdesys = PDESystem(
            [eq], [u(0, x, y) ~ 0], [t ∈ Interval(0.0, 1.0), xdomain, ydomain],
            [t, x, y], [u(t, x, y)]
        )
        v = VariableMap(pdesys, disc)
        @test v.intervals[x] == (0.0, 1.0)
        @test v.intervals[y] == (2.0, 3.0)
        @test v.domains[x] === xdomain
        @test spatial_domain(v) == Rectangle(0.0 .. 1.0, 2.0 .. 3.0)

        @variables w(..)
        @named onedim = PDESystem(
            [Dt(w(t, x)) ~ Dx(Dx(w(t, x)))], [w(0, x) ~ 0], [t ∈ Interval(0.0, 1.0), xdomain],
            [t, x], [w(t, x)]
        )
        @test spatial_domain(VariableMap(onedim, disc)) === xdomain.domain
    end

    @testset "spatial variables split across non-interval domains" begin
        @variables w(..)
        eq3 = Dt(w(t, x, y, z)) ~ Dx(Dx(w(t, x, y, z)))
        @named split = PDESystem(
            [eq3], [w(0, x, y, z) ~ 0],
            [(x, y) ∈ MeshDomain(ToySquareMesh()), z ∈ Interval(0.0, 1.0), t ∈ Interval(0.0, 1.0)],
            [t, x, y, z], [w(t, x, y, z)]
        )
        @test_throws ArgumentError spatial_domain(VariableMap(split, disc))
    end

    @testset "missing domain" begin
        @named missing_domain = PDESystem(
            [eq], [u(0, x, y) ~ 0], [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0)],
            [t, x, y], [u(t, x, y)]
        )
        err = try
            VariableMap(missing_domain, disc)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("Could not find domain for variable y", err.msg)
    end

    @testset "positional constructor without domains" begin
        v = VariableMap([], [x], [], t, Dict(x => (0, 1)), Dict(), [], Dict(), Dict(), Dict())
        @test isempty(v.domains)
        @test v.intervals[x] == (0, 1)
    end
end
