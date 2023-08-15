using Documenter
using CSV
include("../src/EpistemicNetworkAnalysis.jl")

makedocs(
    sitename="EpistemicNetworkAnalysis.jl",
    authors="Mariah A. Knowles <snotskie@gmail.com> and contributors",
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "models.md",
            "rotations.md",
            "plots.md",
            "functions.md"
        ],
        "ICQE23 Workshop" => "icqe23.md"
    ]
)

mv("build", "latest", force=true)
cp("src/index.md", "../README.md", force=true)