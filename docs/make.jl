using Documenter, PDEBase

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "PDEBase.jl",
    authors = "Alex Jones, Chris Rackauckas et al.",
    clean = true,
    doctest = false,
    warnonly = [:missing_docs],
    modules = [PDEBase],
    format = Documenter.HTML(analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/PDEBase/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/PDEBase.jl"; push_preview = true)
