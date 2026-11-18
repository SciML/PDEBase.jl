using PDEBase
using PDEBase: expand_vector_operators, parsed_boundary_map
using PDEBase.StrongForm
using ModelingToolkit: @named, @parameters, PDESystem, get_bcs, get_eqs
using Symbolics
using LinearAlgebra: dot
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

# Evaluates the lowered expression itself, so a changed operation or a truncated operand shows in the value.
evaluate(ex, bindings) = Symbolics.unwrap_const(Symbolics.substitute(Symbolics.unwrap(ex), bindings; fold = Val(true)))

sum_both(a::AbstractVector, b::AbstractVector) = sum(a) + sum(b)
@register_symbolic sum_both(a::AbstractVector, b::AbstractVector)

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

    @testset "divergence and curl need a declared length" begin
        @variables Vu::Vector{Real} Mu::Matrix{Real}
        @parameters F2[1:2] F3[1:3]
        # An unsized operand would drop the components beyond the coordinates once it is bound.
        @test_throws ArgumentError expanded(Dt(U) ~ Symbolics.Divergence([x, y])(Vu); dvs = [U])
        @test_throws ArgumentError expanded(Dt(U) ~ Symbolics.Curl{false}([x, y])(Vu); dvs = [U])
        @test_throws ArgumentError expanded(Dt(U) ~ Symbolics.Divergence([x, y])(Mu * ∇(U)); dvs = [U])
        @test_throws ArgumentError expanded(Dt(U) ~ Symbolics.Divergence([x, y])(∇(U) + Vu); dvs = [U])
        @test_throws ArgumentError expanded(Dt(U) ~ Symbolics.Curl{false}([x, y])(∇(U) + Vu); dvs = [U])
        # A declared length is checked by the operator, so a matching one lowers.
        @test same(expanded(Dt(U) ~ Symbolics.Divergence([x, y])(F2); dvs = [U]).rhs, Dx(F2[1]) + Dy(F2[2]))
        @test same(expanded(Dt(U) ~ Symbolics.Curl{false}([x, y])(F2); dvs = [U]).rhs, Dx(F2[2]) - Dy(F2[1]))
        # Constructor controls: Symbolics itself rejects the oversized operator calls before PDEBase runs.
        @test_throws ArgumentError Symbolics.Divergence([x, y])(F3)
        @test_throws ArgumentError Symbolics.Curl{false}([x, y])(F3)
    end

    @testset "indexing a flux of unknown shape" begin
        @variables Mi::Matrix{Real}
        eq = expanded((Mi * ∇(U))[1] ~ 0; dvs = [U])
        @test !PDEBase._has_vector_operator(eq)
    end

    @testset "dot with an unknown-length vector stays a contraction until it is bound" begin
        @parameters weights::Vector{Real}
        gradient = [19, -3]
        bindings(w) = Dict(
            Symbolics.unwrap(Dx(U)) => gradient[1], Symbolics.unwrap(Dy(U)) => gradient[2],
            Symbolics.unwrap(weights) => w
        )
        residual(eq, w) = evaluate(eq.lhs, bindings(w))
        for lhs in (dot(∇(U), weights), dot(weights, ∇(U)))
            @named sys = PDESystem([lhs ~ 0], boxbcs, box, [t, x, y], [U], [k, weights])
            for eq in (only(get_eqs(expand_vector_operators(sys))), only(get_eqs(parsed_boundary_map(sys, disc).pdesys)))
                @test residual(eq, [7, 13]) == 94
                @test_throws DimensionMismatch residual(eq, [7, 13, 5])
            end
        end
        @variables Z3[1:3]
        @test_throws ArgumentError expanded(dot(∇(U), Z3) ~ 0; dvs = [U])
    end

    @testset "a function named dot in a package named LinearAlgebra is an ordinary function" begin
        # Satisfies a name-and-path check for the standard library: a package called LinearAlgebra that lives
        # under a `stdlib` directory, with its own `dot` whose operands have independent lengths.
        pkg = joinpath(mktempdir(), "stdlib", "LinearAlgebra")
        mkpath(joinpath(pkg, "src"))
        uuid = Base.UUID("11111111-2222-3333-4444-555555555555")
        write(joinpath(pkg, "Project.toml"), "name = \"LinearAlgebra\"\nuuid = \"$uuid\"\nversion = \"0.0.1\"\n")
        write(
            joinpath(pkg, "src", "LinearAlgebra.jl"),
            "module LinearAlgebra\ndot(a::AbstractVector, b::AbstractVector) = sum(a) + sum(b)\nend\n"
        )
        saved_load_path = copy(LOAD_PATH)
        pushfirst!(LOAD_PATH, pkg)
        try
            loader = Module(:LookAlikeLoader)
            Core.eval(loader, :(import LinearAlgebra))
            fake = Base.invokelatest(getproperty, loader, :LinearAlgebra)
            @test fake.dot !== dot && occursin("stdlib", pathof(fake))
            Core.eval(@__MODULE__, :(@register_symbolic $fake.dot(a::AbstractVector, b::AbstractVector)))
            @parameters zq[1:3]
            @named sys = PDESystem([Base.invokelatest(fake.dot, ∇(U), zq) ~ 0], boxbcs, box, [t, x, y], [U], [k, zq])
            eq = only(get_eqs(Base.invokelatest(expand_vector_operators, sys)))
            bindings = Dict(
                Symbolics.unwrap(Dx(U)) => 19, Symbolics.unwrap(Dy(U)) => -3,
                [Symbolics.unwrap(zq[i]) => c for (i, c) in enumerate([7, 13, 5])]...
            )
            @test evaluate(eq.lhs, bindings) == 41
        finally
            copy!(LOAD_PATH, saved_load_path)
        end
    end

    @testset "a registered function named dot is an ordinary function" begin
        mod = Module(:UserDot)
        Core.eval(mod, :(using Symbolics))
        Core.eval(mod, :(dot(a::AbstractVector, b::AbstractVector) = sum(a) + sum(b)))
        Core.eval(mod, :(@register_symbolic dot(a::AbstractVector, b::AbstractVector)))
        @parameters zw[1:3]
        lhs = Base.invokelatest(mod.dot, ∇(U), zw)
        @named sys = PDESystem([lhs ~ 0], boxbcs, box, [t, x, y], [U], [k, zw])
        eq = only(get_eqs(expand_vector_operators(sys)))
        bindings = Dict(
            Symbolics.unwrap(Dx(U)) => 19, Symbolics.unwrap(Dy(U)) => -3,
            [Symbolics.unwrap(zw[i]) => c for (i, c) in enumerate([7, 13, 5])]...
        )
        @test Base.invokelatest(evaluate, eq.lhs, bindings) == 41
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

    @testset "divergence of a spatially varying nonsymmetric tensor flux" begin
        field = x^2 * y + x * y^3
        tensor = [x y; x * y 1 + 2 * x]
        at(ex) = Symbolics.value(Symbolics.substitute(expand_derivatives(ex), Dict(x => 2, y => 3)))
        grad = Symbolics.Grad([x, y])(field)
        flux1 = expanded(Dt(U) ~ (tensor * grad)[1]; dvs = [U])
        flux2 = expanded(Dt(U) ~ (tensor * grad)[2]; dvs = [U])
        divergence = expanded(Dt(U) ~ Symbolics.Divergence([x, y])(tensor * grad); dvs = [U])
        @test (at(flux1.rhs), at(flux2.rhs), at(divergence.rhs)) == (252, 524, 588)
    end

    @testset "elementwise arithmetic with an unknown-length vector stays symbolic" begin
        @parameters pv::Vector{Real}
        gradient = [19, -3]
        bindings(p) = Dict(
            Symbolics.unwrap(Dx(U)) => gradient[1], Symbolics.unwrap(Dy(U)) => gradient[2], Symbolics.unwrap(pv) => p
        )
        @named sys = PDESystem([sum(∇(U) + pv) ~ 0], boxbcs, box, [t, x, y], [U], [k, pv])
        for eq in (only(get_eqs(expand_vector_operators(sys))), only(get_eqs(parsed_boundary_map(sys, disc).pdesys)))
            @test evaluate(eq.lhs, bindings([7, 13])) == 36
            @test_throws DimensionMismatch evaluate(eq.lhs, bindings([7, 13, 5]))
        end
    end

    @testset "tuple-valued boundary conditions are expanded" begin
        gradient = Symbolics.Grad([x])(u(t, 0))
        bcs = (gradient[1] ~ 0, u(t, 1) ~ 0, u(0, x) ~ 0)
        @named sys = PDESystem([Dt(u(t, x)) ~ Dx(Dx(u(t, x)))], bcs, [t ∈ Interval(0.0, 1.0), x ∈ Interval(0.0, 1.0)], [t, x], [u(t, x)], [])
        lowered = get_bcs(expand_vector_operators(sys))
        @test lowered isa Tuple && length(lowered) == 3
        @test isequal(lowered[1].lhs, Dx(u(t, 0)))
        @test isequal(lowered[2], bcs[2])
    end

    @testset "a derivative of an operator result is taken per component" begin
        g = ∇(U)
        eq = expanded(dot(Dt(g), g) ~ 0)
        @test same(eq.lhs, Dt(Dx(U)) * Dx(U) + Dt(Dy(U)) * Dy(U))
    end

    @testset "an unsized vector used only through a scalar index does not defer" begin
        @parameters ps::Vector{Real}
        field = t^2 * x^2 * y + t * x * y^3
        gf = Symbolics.Grad([x, y])(field)
        eq = expanded(dot(Dt(ps[1] * gf), gf) ~ 0; dvs = [U])
        point = Dict(t => 2, x => 3, y => 4)
        # The lowered residual is linear in p[1]: its coefficient is the value for p[1] = 1.
        lowered = Symbolics.substitute(expand_derivatives(eq.lhs), point)
        @test same(lowered, 94160 * ps[1])
        by_hand = Dt(7 * Dx(field)) * Dx(field) + Dt(7 * Dy(field)) * Dy(field)
        @test Symbolics.unwrap_const(Symbolics.substitute(expand_derivatives(by_hand), point)) == 659120 == 7 * 94160
    end

    @testset "a derivative of an array expression with an unsized vector is rejected" begin
        @parameters pq::Vector{Real}
        @test_throws ArgumentError expanded(dot(Dt(∇(U) + pq), ∇(U)) ~ 0; dvs = [U])
    end

    @testset "a scalar-weighted declared vector still contracts after binding" begin
        @parameters pw::Vector{Real} q2[1:2]
        @named sys = PDESystem([dot(∇(U), pw[1] * q2) ~ 0], boxbcs, box, [t, x, y], [U], [k, pw, q2])
        bindings = Dict(
            Symbolics.unwrap(Dx(U)) => 19, Symbolics.unwrap(Dy(U)) => -3, Symbolics.unwrap(pw) => [7, 13],
            Symbolics.unwrap(q2) => [2, -3]
        )
        for eq in (only(get_eqs(expand_vector_operators(sys))), only(get_eqs(parsed_boundary_map(sys, disc).pdesys)))
            @test evaluate(eq.lhs, bindings) == 7 * (19 * 2 + 3 * 3)
        end
    end

    @testset "an index into a composite does not escape the derivative rule" begin
        @parameters pc::Vector{Real}
        g = ∇(U)
        @test_throws ArgumentError expanded(dot(Dt((g + pc)[1] * g), g) ~ 0; dvs = [U])
    end

    @testset "an integer symbol indexes a bare unsized parameter like a literal" begin
        @parameters pi_::Vector{Real} ii::Int
        field = t^2 * x^2 * y + t * x * y^3
        gf = Symbolics.Grad([x, y])(field)
        eq = expanded(dot(Dt(pi_[ii] * gf), gf) ~ 0; dvs = [U])
        lowered = Symbolics.substitute(expand_derivatives(eq.lhs), Dict(t => 2, x => 3, y => 4))
        @test same(lowered, 94160 * pi_[ii])
    end

    @testset "a scalar coefficient does not delay a known dimension error" begin
        @parameters ps2::Vector{Real} q3[1:3]
        @test_throws ArgumentError expanded(dot(∇(U), ps2[1] * q3) ~ 0; dvs = [U])
    end

    @testset "an unsized vector defers every enclosing operation" begin
        @parameters pd::Vector{Real}
        g = ∇(U)
        bindings(p) = Dict(Symbolics.unwrap(Dx(U)) => 19, Symbolics.unwrap(Dy(U)) => -3, Symbolics.unwrap(pd) => p)
        cases = [
            (sum(2 * (g + pd)), 72), (dot(g, g + pd), 464), (sum(abs2.(2 * (g + pd))), 4 * (26^2 + 10^2)),
        ]
        for (lhs, value) in cases
            @named sys = PDESystem([lhs ~ 0], boxbcs, box, [t, x, y], [U], [k, pd])
            for eq in (only(get_eqs(expand_vector_operators(sys))), only(get_eqs(parsed_boundary_map(sys, disc).pdesys)))
                @test evaluate(eq.lhs, bindings([7, 13])) == value
                @test_throws DimensionMismatch evaluate(eq.lhs, bindings([7, 13, 5]))
            end
        end
    end

    @testset "independent-length operands keep their size" begin
        @parameters zv[1:3] zu::Vector{Real}
        for weights in (zv, zu)
            bindings = Dict(
                Symbolics.unwrap(Dx(U)) => 19, Symbolics.unwrap(Dy(U)) => -3,
                (weights === zv ? [Symbolics.unwrap(zv[i]) => c for (i, c) in enumerate([7, 13, 5])] : [Symbolics.unwrap(zu) => [7, 13, 5]])...
            )
            @named sys = PDESystem([sum_both(∇(U), weights) ~ 0], boxbcs, box, [t, x, y], [U], [k, weights])
            @test evaluate(only(get_eqs(expand_vector_operators(sys))).lhs, bindings) == 41
            parsed = parsed_boundary_map(sys, disc)
            @test evaluate(only(get_eqs(parsed.pdesys)).lhs, bindings) == 41
        end
    end

    @testset "divergence of a concrete tensor flux" begin
        Mc = [1.0 0.5; 0.5 2.0]
        eq = expanded(Dt(U) ~ Symbolics.Divergence([x, y])(Mc * ∇(U)); dvs = [U])
        flux = [1.0 * Dx(U) + 0.5 * Dy(U), 0.5 * Dx(U) + 2.0 * Dy(U)]
        @test same(expand_derivatives(eq.rhs), expand_derivatives(Dx(flux[1]) + Dy(flux[2])))
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

    @testset "concrete slices of an operator result" begin
        @test same(expanded(sum(∇(U)[:]) ~ 0).lhs, Dx(U) + Dy(U))
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
