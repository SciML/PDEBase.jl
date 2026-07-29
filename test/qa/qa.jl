using SciMLTesting, PDEBase

run_qa(
    PDEBase;
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            # `ProblemTypeCtx` is the metadata key SciMLBase's `wrap_sol` reads back off
            # `prob.problem_type` to build a `PDETimeSeriesSolution`, so PDEBase cannot
            # substitute a key of its own. Drop once the upstream `public` declaration
            # ships: https://github.com/SciML/ModelingToolkit.jl/pull/4842
            ignore = (:ProblemTypeCtx,),
        ),
    ),
)
