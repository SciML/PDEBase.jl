using Test
using PDEBase

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
        @test Base.ispublic(PDEBase, name)
        @test Base.Docs.doc(PDEBase, name) !== nothing
    end
end
