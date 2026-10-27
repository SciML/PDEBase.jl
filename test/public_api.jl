using Test
using PDEBase

# `Base.ispublic` and the `public` keyword only exist on Julia 1.11+. On the LTS
# `SciMLPublic.@public` expands to nothing, so there is no publicness to query there.
const CAN_QUERY_PUBLIC = isdefined(Base, :ispublic)

# `Base.Docs.doc(mod, name)` reads `name` as a type signature and returns the *module's*
# docstring, so it is true even for names that do not exist. Query the binding instead.
function hasdocstring(mod::Module, name::Symbol)
    return haskey(Base.Docs.meta(mod), Base.Docs.Binding(mod, name))
end

@testset "VariableMap accessor public API" begin
    @test isdefined(PDEBase, :replaced_vars)
    @test hasdocstring(PDEBase, :replaced_vars)
    if CAN_QUERY_PUBLIC
        @test Base.ispublic(PDEBase, :replaced_vars)
    end
end

@testset "Mesh-domain API" begin
    for name in (:MeshDomain, :Region, :regions, :spatial_domain)
        @test hasdocstring(PDEBase, name)
        @test Base.isexported(PDEBase, name)
    end
    for name in (
            :meshdim, :facetset_names, :cellset_names, :contains_point, :backend_mesh,
            :BCKind, :Essential, :Natural, :Robin, :bckind,
        )
        @test hasdocstring(PDEBase, name)
        if CAN_QUERY_PUBLIC
            @test Base.ispublic(PDEBase, name)
        end
    end
end

@testset "Argument classifier API" begin
    for name in (
            :classify_argument, :classify_call, :classify_condition, :substitute_traces, :assert_no_traces,
            :to_box_bcs,
        )
        @test hasdocstring(PDEBase, name)
        @test Base.isexported(PDEBase, name)
    end
    for name in (:ArgumentClass, :Free, :BoxFace, :Piece, :Interior, :ConditionClass, :location, :boxface, :normal)
        @test hasdocstring(PDEBase, name)
        @test !Base.isexported(PDEBase, name)
        if CAN_QUERY_PUBLIC
            @test Base.ispublic(PDEBase, name)
        end
    end
end

@testset "Boundary types" begin
    for name in (:LowerBoundary, :UpperBoundary, :InterfaceBoundary, :HigherOrderInterfaceBoundary, :PieceBoundary)
        @test hasdocstring(PDEBase, name)
        @test Base.isexported(PDEBase, name)
    end
end

@testset "Developer extension API" begin
    developer_api = (
        :parse_bcs,
        :lower_boundaries,
        :interface_errors,
        :check_boundarymap,
        :should_transform,
        :transform_pde_system!,
        :split_ics_bcs,
        :construct_disc_state,
        :construct_discrete_space,
        :construct_var_equation_mapping,
        :construct_differential_discretizer,
        :discretize_equation!,
        :generate_ic_defaults,
        :generate_metadata,
        :generate_system,
        :get_discvars,
        :get_system_unknowns,
        :get_eqvar,
        :add_metadata!,
        :make_pdesys_compatible,
        :supported_conditions,
        :SupportedConditions,
    )

    for name in developer_api
        @test isdefined(PDEBase, name)
        @test hasdocstring(PDEBase, name)
        if CAN_QUERY_PUBLIC
            @test Base.ispublic(PDEBase, name)
        end
    end

    @test !Base.isexported(PDEBase, :parse_bcs)

    # The names a mesh discretization's `parse_bcs` map is built from are public already.
    for name in (:PieceBoundary, :BoundaryPiece, :boundarypieces, :bckind, :BCKind, :Essential, :Natural, :Robin)
        @test hasdocstring(PDEBase, name)
    end

    @test !hasdocstring(PDEBase, :not_a_pdebase_name)

    # A module's docstring lives in the module itself.
    @test hasdocstring(PDEBase.StrongForm, :StrongForm)
    if CAN_QUERY_PUBLIC
        @test Base.ispublic(PDEBase, :StrongForm)
    end
end
