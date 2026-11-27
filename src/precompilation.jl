using PrecompileTools

@setup_workload begin
    # Minimal imports needed for precompilation
    @compile_workload begin
        # Create simple symbolic variables for precompilation
        @variables x y t
        @variables u(..)

        # Get unwrapped symbolic types
        x_sym = SymbolicUtils.unwrap(x)
        y_sym = SymbolicUtils.unwrap(y)
        t_sym = SymbolicUtils.unwrap(t)
        u_term = SymbolicUtils.unwrap(u(x, y, t))

        # Precompile safe_unwrap for common types
        safe_unwrap(u(x, y, t))
        safe_unwrap(x)
        safe_unwrap(1.0)

        # Create a differential term for testing
        Dx = Differential(x)
        term = Dx(u(x, y, t))
        term_sym = SymbolicUtils.unwrap(term)

        # Precompile has_derivatives
        has_derivatives(term_sym)
        has_derivatives(u_term)
        has_derivatives(x_sym)

        # Precompile count_differentials
        count_differentials(term_sym, x_sym)
        count_differentials(u_term, x_sym)

        # Precompile differential_order
        differential_order(term_sym, x_sym)
        differential_order(u_term, x_sym)

        # Precompile get_depvars
        depvar_ops = [operation(u_term)]
        get_depvars(term_sym, depvar_ops)
        get_depvars(u_term, depvar_ops)

        # Precompile split_terms
        eq = SymbolicUtils.unwrap(Dx(u(x, y, t)) ~ u(x, y, t))
        split_terms(eq)

        # Precompile split_additive_terms
        split_additive_terms(eq)

        # Precompile recursive_unwrap
        recursive_unwrap(term_sym)
    end
end

# Private stand-ins for a discretizer package's types, so the entry points below run on
# the concrete discretization subtypes that downstream packages define.
struct _WarmupBoxDiscretization <: AbstractBoxDiscretization
    time::Any
end
get_time(disc::_WarmupBoxDiscretization) = disc.time

struct _WarmupMesh end
meshdim(::_WarmupMesh) = 2
facetset_names(::_WarmupMesh) = (:inlet, :wall)
cellset_names(::_WarmupMesh) = (:all,)
contains_point(::_WarmupMesh, p) = all(c -> 0 <= c <= 1, p)

struct _WarmupMeshDiscretization <: AbstractMeshDiscretization
    time::Any
end
get_time(disc::_WarmupMeshDiscretization) = disc.time

@setup_workload begin
    @variables t x y
    @variables u(..)
    Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
    unit_interval = DomainSets.Interval(0.0, 1.0)
    box_disc = _WarmupBoxDiscretization(t)
    U = u(t, x, y)
    # Lazy vector operators exist only in newer Symbolics releases.
    rhs = isdefined(Symbolics, :Laplacian) ? Symbolics.Laplacian([x, y])(U) : Dx(Dx(U)) + Dy(Dy(U))
    bcs = [u(0, x, y) ~ 0, u(t, 0, y) ~ 0, u(t, 1, y) ~ 0, u(t, x, 0) ~ 0, u(t, x, 1) ~ 0]
    heat2d = PDESystem(
        [Dt(U) ~ rhs], bcs, [iv ∈ unit_interval for iv in (t, x, y)], [t, x, y], [U], [], name = :warmup_heat2d
    )

    @compile_workload begin
        parsed = parsed_boundary_map(heat2d, box_disc)
        StrongForm.analyse(parsed.pdesys, parsed.v, parsed.boundarymap)
    end
end

@setup_workload begin
    @variables t x y k
    @variables u(..)
    Dt, Dx, Dy = Differential(t), Differential(x), Differential(y)
    Ω = MeshDomain(_WarmupMesh())
    Γin, Γwall = boundarypiece(Ω, :inlet), boundarypiece(Ω, :wall)
    Dn = NormalDerivative(Γwall)
    U = u(t, x, y)
    bcs = [u(t, Γin) ~ sin(t), Dn(u(t, Γwall)) ~ 0, u(0, x, y) ~ 0]
    meshsys = PDESystem(
        [Dt(U) ~ k * (Dx(Dx(U)) + Dy(Dy(U)))], bcs,
        [(x, y) ∈ Ω, t ∈ DomainSets.Interval(0.0, 1.0)], [t, x, y], [U], [k], name = :warmup_mesh
    )
    mesh_disc = _WarmupMeshDiscretization(t)

    @compile_workload begin
        parsed = parsed_boundary_map(meshsys, mesh_disc)
        StrongForm.analyse(parsed.pdesys, parsed.v, parsed.boundarymap)
    end
end
