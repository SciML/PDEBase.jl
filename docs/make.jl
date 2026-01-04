using Documenter, PDEBase

# Create assets directory if needed
mkpath("./src/assets")

# Copy manifest/project files if they exist
isfile("./Manifest.toml") && cp("./Manifest.toml", "./src/assets/Manifest.toml", force = true)
isfile("./Project.toml") && cp("./Project.toml", "./src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(
    sitename = "PDEBase.jl",
    authors = "Alex Jones, Chris Rackauckas et al.",
    clean = true,
    doctest = false,
    warnonly = [:missing_docs],
    modules = [PDEBase],
    format = Documenter.HTML(
        analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/PDEBase/stable/"
    ),
    pages = pages
)

deploydocs(repo = "github.com/SciML/PDEBase.jl"; push_preview = true)
