using Documenter
using CSV
include("../src/EpistemicNetworkAnalysis.jl")

# BUGFIX
# https://github.com/JuliaDocs/Documenter.jl/issues/1414#issuecomment-761891250
ENV["GKSwstype"] = "100"

makedocs(
    sitename="EpistemicNetworkAnalysis.jl",
    authors="Mariah A. Knowles <snotskie@gmail.com> and contributors",
    format=Documenter.HTML(assets=[
        "assets/script.js"
    ]),
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "models.md",
            "rotations.md",
            "plots.md",
            "functions.md"
        ],
        "Examples" => "examples/" .* sort(readdir("src/examples")),
        "ICQE23 Workshop" => "icqe23.md"
    ]
)

mv("build", "latest", force=true)
cp("src/index.md", "../README.md", force=true)