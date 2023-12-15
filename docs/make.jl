using Documenter
using DormandPrince
using DocThemeIndigo

indigo = DocThemeIndigo.install(Configurations)

makedocs(;
    modules = [DormandPrince],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://queracomputing.github.io/DormandPrince.jl",
        assets=String[indigo],
    ),
    pages = [
        "Home" => "index.md",
    ],
    repo = "https://github.com/QuEraComputing/DormandPrince.jl",
    sitename = "DormandPrince.jl",
)

deploydocs(; repo = "https://github.com/QuEraComputing/DormandPrince.jl")
