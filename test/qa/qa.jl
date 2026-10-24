using SciMLTesting, PDEBase

run_qa(
    PDEBase;
    # DomainSets declares no `public` names: `HyperRectangle` is the supertype of its boxes,
    # and `indomain` is the method it documents for new domain types to implement membership.
    ei_kwargs = (; all_qualified_accesses_are_public = (; ignore = (:HyperRectangle, :indomain))),
    # Upstream syntax PDEBase defines only until Symbolics or DomainSets does (src/fallbacks.jl).
    reexports_allow = (:BoundaryPiece, :boundarypiece, :boundarypieces, :NormalDerivative),
)
