using PDEBase
using Test

@testset "Documentation project" begin
    version = pkgversion(PDEBase)
    docs_project = read(joinpath(pkgdir(PDEBase), "docs", "Project.toml"), String)
    @test contains(docs_project, "PDEBase = \"$(version.major).$(version.minor)\"")
end
