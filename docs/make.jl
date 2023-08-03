using Documenter
include("../src/EpistemicNetworkAnalysis.jl")

```@meta
CurrentModule = EpistemicNetworkAnalysis
```

makedocs(
    sitename="EpistemicNetworkAnalysis.jl",
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "models.md",
            "rotations.md",
            "plots.md",
            "functions.md"
        ]
    ]
)

mv("build", "latest", force=true)