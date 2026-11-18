using PDEBase
using PDEBase: expand_vector_operators, parsed_boundary_map
using PDEBase.StrongForm
using ModelingToolkit: @named, @parameters, PDESystem, get_bcs, get_eqs
using Symbolics
const dot = Symbolics.LinearAlgebra.dot
using SymbolicUtils: operation
using SciMLBase
using Test
import DomainSets
using DomainSets: Interval

struct VectorOpsMesh end
PDEBase.meshdim(::VectorOpsMesh) = 2
PDEBase.facetset_names(::VectorOpsMesh) = (:wall,)
PDEBase.cellset_names(::VectorOpsMesh) = (:all,)
PDEBase.contains_point(::VectorOpsMesh, p) = true

struct VectorOpsBoxDiscretization <: PDEBase.AbstractBoxDiscretization
    time::Any
end
PDEBase.get_time(disc::VectorOpsBoxDiscretization) = disc.time

struct VectorOpsOptimizationDiscretization <: PDEBase.AbstractOptimizationSystemDiscretization end
struct OptimizationReached <: Exception
    eqs::Any
end
PDEBase.construct_disc_state(::VectorOpsOptimizationDiscretization) = nothing
function PDEBase.construct_discrete_space(v::VariableMap, pdesys::PDESystem, ::VectorOpsOptimizationDiscretization)
    throw(OptimizationReached(get_eqs(pdesys)))
end

same(a, b) = isequal(Symbolics.expand(Num(a) - Num(b)), Num(0))

if !isdefined(Symbolics, :Grad)
    @testset "vector operators need Symbolics 7.49" begin
        @test_skip false
    end
else
    @parameters k
    @variables t x y z
    @variables u(..) v(..) w(..) A(..)[1:2]
    Dt, Dx, Dy, Dz = Differential(t), Differential(x), Differential(y), Differential(z)
    U = u(t, x, y)
    ∇ = Symbolics.Grad([x, y])
    box = [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0), y ∈ Interval(0.0, 1.0)]
    boxbcs = [u(0, x, y) ~ 0, u(t, 0, y) ~ 0, u(t, 1, y) ~ 0, u(t, x, 0) ~ 0, u(t, x, 1) ~ 0]
    disc = VectorOpsBoxDiscretization(t)

    function expanded(eq; ivs = [t, x, y], dvs = [U], bcs = boxbcs, domain = box)
        @named sys = PDESystem([eq], bcs, domain, ivs, dvs, [k])
        return only(get_eqs(expand_vector_operators(sys)))
    end

    @testset "Laplacian" begin
        eq = expanded(Dt(U) ~ k * Symbolics.Laplacian([x, y])(U))
        @test same(eq.lhs, Dt(U))
        @test same(eq.rhs, k * (Dx(Dx(U)) + Dy(Dy(U))))
    end

    @testset "divergence form is kept" begin
        eq = expanded(Dt(U) ~ Symbolics.Divergence([x, y])(k * ∇(U)))
        @test same(eq.rhs, Dx(k * Dx(U)) + Dy(k * Dy(U)))
        eq = expanded(Dt(U) ~ Symbolics.Divergence([x, y])((1 + x^2) * ∇(U)))
        @test same(eq.rhs, Dx((1 + x^2) * Dx(U)) + Dy((1 + x^2) * Dy(U)))
    end

    @testset "squared gradient" begin
        eq = expanded(Dt(U) ~ ∇(U)[1]^2 + ∇(U)[2]^2)
        @test same(eq.rhs, Dx(U)^2 + Dy(U)^2)
    end

    @testset "2D curl" begin
        V = v(t, x, y)
        eq = expanded(Dt(U) ~ Symbolics.Curl{false}([x, y])([U, V]); dvs = [U, V])
        @test same(eq.rhs, Dx(V) - Dy(U))
    end

    @testset "3D curl" begin
        a, b, c = u(t, x, y, z), v(t, x, y, z), w(t, x, y, z)
        eq = expanded(
            Dt(a) ~ Symbolics.Curl{false}([x, y, z])([a, b, c])[1];
            ivs = [t, x, y, z], dvs = [a, b, c], bcs = Equation[],
            domain = [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0), y ∈ Interval(0.0, 1.0), z ∈ Interval(0.0, 1.0)]
        )
        @test same(eq.rhs, Dy(c) - Dz(b))
    end

    @testset "vector unknown" begin
        Av = A(t, x, y)
        eq = expanded(Dt(U) ~ Symbolics.Divergence([x, y])(Av); dvs = [U, Av])
        @test same(eq.rhs, Dx(Av[1]) + Dy(Av[2]))
    end

    @testset "array-valued equation" begin
        V, W = v(t, x, y), w(t, x, y)
        @named sys = PDESystem([∇(U) ~ [V, W]], boxbcs, box, [t, x, y], [U, V, W], [k])
        eqs = get_eqs(expand_vector_operators(sys))
        @test length(eqs) == 2
        @test isequal(eqs[1], Dx(U) ~ V)
        @test isequal(eqs[2], Dy(U) ~ W)
        @test !PDEBase._has_vector_operator(eqs)

        parsed = parsed_boundary_map(sys, disc)
        @test length(get_eqs(parsed.pdesys)) == 2
        analysis = analyse(parsed.pdesys, parsed.v, parsed.boundarymap)
        @test analysis.nfields == 3
    end

    @testset "array equation with mismatched sides is an error" begin
        @variables B[1:3]
        @named sys = PDESystem([Equation(∇(U), B)], boxbcs, box, [t, x, y], [U], [k])
        @test_throws ArgumentError expand_vector_operators(sys)
    end

    @testset "the description of the system is kept" begin
        @named sys = PDESystem(
            [Dt(U) ~ Symbolics.Laplacian([x, y])(U)], boxbcs, box, [t, x, y], [U], [k];
            description = "Heat equation"
        )
        @test expand_vector_operators(sys).description == "Heat equation"
    end

    @testset "complex-valued sides are left for the complex handling" begin
        eq = expanded(Equation(Dt(U), im * Symbolics.Laplacian([x, y])(U)))
        @test !PDEBase._has_vector_operator(eq)
        @test string(eq.rhs) == string(Symbolics.unwrap(im * (Dx(Dx(U)) + Dy(Dy(U)))))
    end

    @testset "complex-valued boundary pairs" begin
        pair = u(0, x, y) => im * Symbolics.Laplacian([x, y])(U)
        @test PDEBase._has_vector_operator(pair)
        lowered = PDEBase._expand_condition(pair, unwrap.([t, x, y]))
        @test !PDEBase._has_vector_operator(lowered)
    end

    @testset "array operations on operator results" begin
        g = ∇(U)
        eq = expanded(2 * g[1] + dot(g, g) ~ 0)
        @test same(eq.lhs, 2 * Dx(U) + Dx(U)^2 + Dy(U)^2)
        eq = expanded((2 * g)[1] ~ 0)
        @test same(eq.lhs, 2 * Dx(U))
        @test !occursin("array_literal", string(eq))
        eq = expanded(sum(g) ~ 0)
        @test same(eq.lhs, Dx(U) + Dy(U))

        @named sys = PDESystem([(2 * g)[1] + dot(g, g) ~ 0], boxbcs, box, [t, x, y], [U], [k])
        parsed = parsed_boundary_map(sys, disc)
        @test analyse(parsed.pdesys, parsed.v, parsed.boundarymap).nfields == 1
    end

    @testset "broadcast followed by a reduction" begin
        @variables Ac[1:2]
        g = ∇(U)
        gx, gy = Dx(U), Dy(U)
        @test same(expanded(sum(g .^ 2) ~ 0).lhs, gx^2 + gy^2)
        @test same(expanded(sum(g .* g) ~ 0).lhs, gx^2 + gy^2)
        @test same(expanded(sum(sin.(g)) ~ 0).lhs, sin(gx) + sin(gy))
        @test same(expanded(dot(Ac, g) ~ 0).lhs, Ac[1] * gx + Ac[2] * gy)
        @test same(expanded(sum(g .^ 2) ~ 0; dvs = [U, Ac]).lhs, gx^2 + gy^2)
    end

    @testset "divergence of a vector of unknown length" begin
        @variables Vu::Vector{Real}
        eq = expanded(Dt(U) ~ Symbolics.Divergence([x, y])(Vu); dvs = [U])
        @test same(eq.rhs, Dx(Vu[1]) + Dy(Vu[2]))
    end

    @testset "operators nested in a vector of unknown length" begin
        @variables M::Matrix{Real}
        eq = expanded(Dt(U) ~ Symbolics.Divergence([x, y])(M * ∇(U)); dvs = [U])
        @test !PDEBase._has_vector_operator(eq)
        # Both flux components, in order, against the same product built by hand.
        Mflux = M * [Dx(U), Dy(U)]
        @test isequal(eq.rhs, Symbolics.unwrap(Dx(Mflux[1]) + Dy(Mflux[2])))
        @test !isequal(eq.rhs, Symbolics.unwrap(Dx(Mflux[2]) + Dy(Mflux[1])))
    end

    @testset "indexing a flux of unknown shape" begin
        @variables Mi::Matrix{Real}
        eq = expanded((Mi * ∇(U))[1] ~ 0; dvs = [U])
        @test !PDEBase._has_vector_operator(eq)
    end

    @testset "contraction with a vector of unknown length" begin
        @variables Zu::Vector{Real}
        eq = expanded(dot(∇(U), Zu) ~ 0; dvs = [U])
        @test same(eq.lhs, Dx(U) * Zu[1] + Dy(U) * Zu[2])
        eq = expanded(dot(Zu, ∇(U)) ~ 0; dvs = [U])
        @test same(eq.lhs, Dx(U) * Zu[1] + Dy(U) * Zu[2])
    end

    @testset "tensor-weighted gradient equations" begin
        A = [1 2; 3 4]
        g = ∇(U)
        gx, gy = Dx(U), Dy(U)
        reference = gx^2 + 5gx * gy + 4gy^2
        explicit = [reference ~ 0]
        @named refsys = PDESystem(explicit, boxbcs, box, [t, x, y], [U], [k])
        refparsed = parsed_boundary_map(refsys, disc)
        refanalysis = analyse(refparsed.pdesys, refparsed.v, refparsed.boundarymap)
        for lhs in (g' * A * g, dot(g, A * g))
            @test same(expanded(lhs ~ 0).lhs, reference)
            @named sys = PDESystem([lhs ~ 0], boxbcs, box, [t, x, y], [U], [k])
            parsed = parsed_boundary_map(sys, disc)
            @test length(parsed.boundarymap) == length(refparsed.boundarymap)
            analysis = analyse(parsed.pdesys, parsed.v, parsed.boundarymap)
            @test (analysis.nfields, analysis.dim) == (refanalysis.nfields, refanalysis.dim)
        end
    end

    @testset "tensor-weighted gradient equations with a symbolic tensor" begin
        @parameters Mp[1:2, 1:2]
        g = ∇(U)
        gx, gy = Dx(U), Dy(U)
        reference = Mp[1, 1] * gx^2 + (Mp[1, 2] + Mp[2, 1]) * gx * gy + Mp[2, 2] * gy^2
        for lhs in (g' * Mp * g, dot(g, Mp * g))
            @test same(expanded(lhs ~ 0).lhs, reference)
            @named sys = PDESystem([lhs ~ 0], boxbcs, box, [t, x, y], [U], [k, Mp])
            parsed = parsed_boundary_map(sys, disc)
            @test analyse(parsed.pdesys, parsed.v, parsed.boundarymap).nfields == 1
        end
    end

    @testset "divergence of a nonsymmetric tensor flux" begin
        field = x^2 * y + x * y^3
        eq = expanded(Dt(U) ~ Symbolics.Divergence([x, y])([1 2; 3 4] * Symbolics.Grad([x, y])(field)); dvs = [U])
        @test Symbolics.value(Symbolics.substitute(expand_derivatives(eq.rhs), Dict(x => 2, y => 3))) == 305
    end

    @testset "divergence of a concrete tensor flux" begin
        Mc = [1.0 0.5; 0.5 2.0]
        eq = expanded(Dt(U) ~ Symbolics.Divergence([x, y])(Mc * ∇(U)); dvs = [U])
        flux = [1.0 * Dx(U) + 0.5 * Dy(U), 0.5 * Dx(U) + 2.0 * Dy(U)]
        @test same(expand_derivatives(eq.rhs), expand_derivatives(Dx(flux[1]) + Dy(flux[2])))
    end

    @testset "vector of unknown length" begin
        @variables Av2::Vector{Real}
        eq = expanded(Dt(U) ~ Symbolics.Curl{false}([x, y])(Av2); dvs = [U])
        @test same(eq.rhs, Dx(Av2[2]) - Dy(Av2[1]))
    end

    @testset "nested equation groups keep their shape" begin
        V, W = v(t, x, y), w(t, x, y)
        g = ∇(U)
        plain = Dt(V) ~ Dx(W)
        @named sys = PDESystem([[g[1] ~ V, g[2] ~ W], plain], boxbcs, box, [t, x, y], [U, V, W], [k])
        eqs = get_eqs(expand_vector_operators(sys))
        @test eqs[1] isa AbstractVector && length(eqs[1]) == 2
        @test isequal(eqs[1][1], Dx(U) ~ V)
        @test isequal(eqs[1][2], Dy(U) ~ W)
        @test isequal(eqs[2], plain)
    end

    @testset "symbolic index into an operator result is a clear error" begin
        @variables i::Int
        @named sys = PDESystem([∇(U)[i] ~ v(t, x, y)], boxbcs, box, [t, x, y], [U, v(t, x, y)], [k])
        @test_throws ArgumentError expand_vector_operators(sys)
    end

    @testset "d_orders and the pipeline see the expanded equation" begin
        heat = Dt(U) ~ k * Symbolics.Laplacian([x, y])(U)
        @named sys = PDESystem([heat], boxbcs, box, [t, x, y], [U], [k])
        @test isempty(d_orders(x, [heat]))
        lowered = expand_vector_operators(sys)
        @test d_orders(x, get_eqs(lowered)) == [2]
        @test d_orders(y, get_eqs(lowered)) == [2]

        parsed = parsed_boundary_map(sys, disc)
        @test d_orders(x, get_eqs(parsed.pdesys)) == [2]
        @test d_orders(y, get_eqs(parsed.pdesys)) == [2]
        analysis = analyse(parsed.pdesys, parsed.v, parsed.boundarymap)
        @test (analysis.nfields, analysis.dim) == (1, 2)
        @test same(analysis.strong[1], analysis.slots.Ud[1] - analysis.slots.P[1] * (analysis.slots.H[1] + analysis.slots.H[4]))
    end

    @testset "optimisation symbolic_discretize path" begin
        @named sys = PDESystem(
            [k * Symbolics.Laplacian([x, y])(U) ~ 0], boxbcs, box, [t, x, y], [U], [k]
        )
        err = try
            SciMLBase.symbolic_discretize(sys, VectorOpsOptimizationDiscretization())
            nothing
        catch e
            e
        end
        @test err isa OptimizationReached
        @test d_orders(x, err.eqs) == [2]
    end

    @testset "a system without the operators is returned unchanged" begin
        @named sys = PDESystem([Dt(U) ~ k * (Dx(Dx(U)) + Dy(Dy(U)))], boxbcs, box, [t, x, y], [U], [k])
        @test expand_vector_operators(sys) === sys
    end

    @testset "other equations of the system keep their form" begin
        plain = Dt(U) ~ Dx(k * Dx(U))
        @named sys = PDESystem(
            [plain, Dt(v(t, x, y)) ~ Symbolics.Laplacian([x, y])(v(t, x, y))], boxbcs, box,
            [t, x, y], [U, v(t, x, y)], [k]
        )
        eqs = get_eqs(expand_vector_operators(sys))
        @test isequal(eqs[1], plain)
        @test same(eqs[2].rhs, Dx(Dx(v(t, x, y))) + Dy(Dy(v(t, x, y))))
    end

    @testset "variables must be independent variables of the system" begin
        @named sys = PDESystem(
            [Dt(U) ~ Symbolics.Laplacian([x, z])(U)], boxbcs, box, [t, x, y], [U], [k]
        )
        @test_throws ArgumentError expand_vector_operators(sys)
        @test_throws ArgumentError parsed_boundary_map(sys, disc)
    end

    @testset "traces in boundary conditions" begin
        Ω = MeshDomain(VectorOpsMesh())
        Γ = boundarypiece(Ω, :wall)
        domains = [(x, y) ∈ Ω, t ∈ Interval(0.0, 1.0)]
        gradient = Symbolics.Grad([x, y])(u(t, Γ))
        bcs = [gradient[1] ~ 0, Dt(u(t, Γ)) + gradient[2] ~ 0, u(0, x, y) ~ 0]
        @named sys = PDESystem([Dt(U) ~ k * Symbolics.Laplacian([x, y])(U)], bcs, domains, [t, x, y], [U], [k])
        lowered = expand_vector_operators(sys)
        trace = u(t, Γ)
        @test isequal(get_bcs(lowered)[1].lhs, Dx(trace))
        @test same(get_bcs(lowered)[2].lhs, Dt(trace) + Dy(trace))
        @test isequal(get_bcs(lowered)[3], bcs[3])
        @test same(only(get_eqs(lowered)).rhs, k * (Dx(Dx(U)) + Dy(Dy(U))))
    end
end
