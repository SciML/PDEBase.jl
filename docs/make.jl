using Documenter, PDEBase

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

# Make sure that plots don't throw a bunch of warnings / errors!
ENV["GKSwstype"] = "100"

include("pages.jl")

makedocs(sitename = "PDEBase.jl",
    authors = "Alex Jones, Chris Rackauckas et al.",
    clean = true,
    doctest = false,
    strict = [
        :doctest,
        :linkcheck,
        :parse_error,
        :example_block        # Other available options are        # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
    ],
    modules = [MethodOfLines],
    format = Documenter.HTML(analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/PDEBase/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/PDEBase.jl"; push_preview = true)
