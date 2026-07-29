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

@testset "Developer extension API" begin
    developer_api = (
        :interface_errors,
        :check_boundarymap,
        :should_transform,
        :transform_pde_system!,
        :construct_disc_state,
        :construct_discrete_space,
        :construct_var_equation_mapping,
        :construct_differential_discretizer,
        :discretize_equation!,
        :generate_ic_defaults,
        :generate_metadata,
        :generate_system,
        :get_discvars,
        :get_eqvar,
        :add_metadata!,
    )

    for name in developer_api
        @test isdefined(PDEBase, name)
        @test hasdocstring(PDEBase, name)
        if CAN_QUERY_PUBLIC
            @test Base.ispublic(PDEBase, name)
        end
    end

    @test !hasdocstring(PDEBase, :not_a_pdebase_name)
end
