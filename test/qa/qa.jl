using SciMLTesting, PDEBase

run_qa(
    PDEBase;
    # DomainSets declares no `public` names, and `HyperRectangle` is the supertype of its boxes.
    ei_kwargs = (; all_qualified_accesses_are_public = (; ignore = (:HyperRectangle,))),
    # Upstream syntax PDEBase defines only until Symbolics or DomainSets does (src/fallbacks.jl).
    reexports_allow = (:BoundaryPiece, :boundarypiece, :boundarypieces, :NormalDerivative),
)
