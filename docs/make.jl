cd("/home/jovyan/docs/")
push!(LOAD_PATH,"../src/")
using EpistemicNetworkAnalysis
using Pkg
Pkg.add("Documenter")
using Documenter

makedocs(
    sitename="EpistemicNetworkAnalysis.jl",
    modules=[EpistemicNetworkAnalysis],
    pages=[
        "Home" => "index.md"
    ]
)