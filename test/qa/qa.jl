using SciMLTesting, PDEBase, Test
using JET

run_qa(
    PDEBase;
    explicit_imports = true,
    # symbolic_discretize(::PDESystem, ::AbstractDiscretization) extends
    # SciMLBase.symbolic_discretize as PDEBase's discretization entry point; the
    # function is conceptually owned by the discretization stack PDEBase implements.
    aqua_kwargs = (;
        piracies = (; treat_as_own = [PDEBase.symbolic_discretize]),
        deps_compat = (; check_extras = false),
    ),
    ei_kwargs = (;
        all_explicit_imports_via_owners = (;
            ignore = (
                # re-exported by ModelingToolkit; owner ModelingToolkitBase
                :ProblemTypeCtx, :get_bcs, :get_connector_type, :get_domain,
                :get_dvs, :get_eqs, :get_gui_metadata, :get_iv, :get_ivs,
                :get_metadata, :get_ps, :get_systems, :get_unknowns,
                # re-exported by SymbolicUtils; owner TermInterface
                :maketerm, :metadata,
                # re-exported by Symbolics; owner SymbolicUtils
                :unwrap,
            ),
        ),
        all_qualified_accesses_via_owners = (;
            ignore = (
                :ExtraVariablesSystemException,  # ModelingToolkit re-export; owner ModelingToolkitBase
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                :ExtraVariablesSystemException,  # ModelingToolkit non-public
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                :AbstractDiscretization, :AbstractDiscretizationMetadata,  # SciMLBase non-public
                :Chain, :Prewalk, :maketerm, :metadata,                    # SymbolicUtils non-public
                :ProblemTypeCtx,                                           # ModelingToolkit non-public
                :rename, :setname,                                         # Symbolics non-public
            ),
        ),
    ),
    # PDEBase's public surface rests on heavy `using ModelingToolkit / Symbolics /
    # SymbolicUtils / SciMLBase / DomainSets`; making each implicitly used name
    # explicit is a large, risky refactor tracked separately. See the QA PR.
    ei_broken = (:no_implicit_imports,),
)
