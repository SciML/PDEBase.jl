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

struct ToyBareMesh end
PDEBase.meshdim(::ToyBareMesh) = 2
PDEBase.facetset_names(::ToyBareMesh) = ()

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

@testset "MeshDomain boundary pieces" begin
    Ω = MeshDomain(ToySquareMesh())
    Γ = boundarypiece(Ω, :wall)
    @test Γ isa BoundaryPiece
    @test parent(Γ) === Ω
    @test Γ.name === :wall
    @test Γ.geometry === nothing
    @test isequal(Γ, boundarypiece(Ω, :wall))
    err = try
        boundarypiece(Ω, :nowhere)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test all(name -> occursin(string(name), err.msg), (:nowhere, :inlet, :wall, :outlet))

    pieces = boundarypieces(Ω)
    @test [p.name for p in pieces] == [:inlet, :wall, :outlet]
    @test all(p -> parent(p) === Ω && p.geometry === nothing, pieces)
    ∂Ω = DomainSets.boundary(Ω)
    @test ∂Ω isa DomainSets.UnionDomain
    @test isequal(collect(DomainSets.components(∂Ω)), pieces)

    bare = MeshDomain(ToyBareMesh())
    @test isempty(boundarypieces(bare))
    ∂bare = DomainSets.boundary(bare)
    @test ∂bare isa DomainSets.EmptySpace
    @test eltype(∂bare) === eltype(bare)
    @test [0.5, 0.5] ∉ ∂bare
end

@testset "bckind" begin
    @parameters t x y k α
    @variables u(..)
    Dx = Differential(x)
    Ω = MeshDomain(ToySquareMesh())
    Γ_in, Γ_wall, Γ_out = boundarypiece.(Ref(Ω), (:inlet, :wall, :outlet))
    Dn_wall, Dn_out = NormalDerivative(Γ_wall), NormalDerivative(Γ_out)

    @test PDEBase.bckind(u(t, Γ_in) ~ sin(t)) === PDEBase.Essential()
    @test PDEBase.bckind(sin(t) ~ 2 * u(t, Γ_in)) === PDEBase.Essential()
    @test PDEBase.bckind(Dn_wall(u(t, Γ_wall)) ~ 0) === PDEBase.Natural()
    @test PDEBase.bckind(k * Dx(u(t, Γ_wall)) ~ 1) === PDEBase.Natural()
    @test PDEBase.bckind(k * Dn_out(u(t, Γ_out)) + α * u(t, Γ_out) ~ α * 0.3) === PDEBase.Robin()
    @test PDEBase.Essential() isa PDEBase.BCKind
    @test PDEBase.Natural() isa PDEBase.BCKind
    @test PDEBase.Robin() isa PDEBase.BCKind
    @test_throws ArgumentError PDEBase.bckind(u(t, 0, y) ~ 0)
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

        Γ = boundarypiece(Ω, :inlet)
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

    @testset "time inside a tuple pairing" begin
        @variables w(..)
        @named spacetime = PDESystem(
            [Dt(w(t, x)) ~ Dx(Dx(w(t, x)))], [w(0, x) ~ 0],
            [(t, x) ∈ Rectangle(0.0 .. 2.0, 10.0 .. 20.0)], [t, x], [w(t, x)]
        )
        v = VariableMap(spacetime, disc)
        @test isequal(v.x̄, [x])
        @test v.intervals[t] == (0.0, 2.0)
        @test v.intervals[x] == (10.0, 20.0)
        @test spatial_domain(v) == Interval(10.0, 20.0)

        box = DomainSets.ProductDomain(0.0 .. 2.0, 10.0 .. 20.0, -1.0 .. 1.0)
        @named spacetime2 = PDESystem(
            [eq], [u(0, x, y) ~ 0], [(t, x, y) ∈ box], [t, x, y], [u(t, x, y)]
        )
        @test spatial_domain(VariableMap(spacetime2, disc)) ==
            Rectangle(10.0 .. 20.0, -1.0 .. 1.0)

        @named meshtime = PDESystem(
            [Dt(w(t, x)) ~ Dx(Dx(w(t, x)))], [w(0, x) ~ 0],
            [(t, x) ∈ MeshDomain(ToySquareMesh())], [t, x], [w(t, x)]
        )
        err = try
            spatial_domain(VariableMap(meshtime, disc))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("MeshDomain", err.msg)
    end

    @testset "box axis order follows the declaration" begin
        @variables q(..)
        Ix, Iy = Interval(0.0, 1.0), Interval(5.0, 7.0)
        coupled = [
            Dt(q(t, y)) ~ Dy(Dy(q(t, y))),
            Dt(u(t, x, y)) ~ Dx(Dx(u(t, x, y))) + q(t, y),
        ]
        bcs = [q(0, y) ~ 0, u(0, x, y) ~ 0]
        @named split_xy = PDESystem(
            coupled, bcs, [t ∈ Interval(0.0, 1.0), x ∈ Ix, y ∈ Iy], [t, x, y], [q(t, y), u(t, x, y)]
        )
        v = VariableMap(split_xy, disc)
        @test isequal(v.x̄, [y, x])
        @test isequal(PDEBase._box_axes(v), [x, y])
        @test spatial_domain(v) == Rectangle(Ix, Iy)

        @named split_yx = PDESystem(
            coupled, bcs, [t ∈ Interval(0.0, 1.0), y ∈ Iy, x ∈ Ix], [t, x, y], [q(t, y), u(t, x, y)]
        )
        v = VariableMap(split_yx, disc)
        @test isequal(PDEBase._box_axes(v), [y, x])
        @test spatial_domain(v) == Rectangle(Iy, Ix)

        box = DomainSets.ProductDomain(Interval(0.0, 2.0), Ix, Iy)
        @named spacetime = PDESystem(
            [Dt(u(t, y, x)) ~ Dx(Dx(u(t, y, x)))], [u(0, y, x) ~ 0], [(t, x, y) ∈ box],
            [t, x, y], [u(t, y, x)]
        )
        v = VariableMap(spacetime, disc)
        @test isequal(v.x̄, [y, x])
        @test isequal(PDEBase._box_axes(v), [x, y])
        @test spatial_domain(v) == Rectangle(Ix, Iy)
    end

    @testset "no spatial variables" begin
        @variables w(..)
        @named ode = PDESystem(
            [Dt(w(t)) ~ -w(t)], [w(0) ~ 1], [t ∈ Interval(0.0, 1.0)], [t], [w(t)]
        )
        v = VariableMap(ode, disc)
        @test isempty(v.x̄)
        @test spatial_domain(v) === nothing
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
        @test_throws ArgumentError spatial_domain(v)
        @test v.intervals[x] == (0, 1)
    end
end
