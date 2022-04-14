module EpistemicNetworkAnalysis

# Imports
## Plotting
using Plots
using Plots.PlotMeasures
import Plots.plot
import Plots.plot!
import Plots.Plot
using Colors
using Dierckx

## Data
using DataFrames
using CSV

## Math
using Statistics
using LinearAlgebra
using MultivariateStats
using HypothesisTests
using GLM

## Nonlinear
using Random
using Distances
using UMAP

# Includes
include("./helpers.jl")
include("./typetree.jl")
include("./plotbase.jl")
include("./SVDRotation.jl")
include("./LDARotation.jl")
include("./MulticlassRotation.jl")
include("./FormulaRotation.jl")
include("./Formula2Rotation.jl")
include("./MeansRotation.jl")
include("./Means2Rotation.jl")
include("./ThematicRotation.jl")
include("./UMAPRotation.jl")
include("./ENAModel.jl")
include("./DigraphENAModel.jl")
include("./NonlinearENAModel.jl")
# include("./CopyRotation.jl")
include("./ena_dataset.jl")

# Exports
export ENAModel
export DigraphENAModel
export NonlinearENAModel
export plot
export SVDRotation
export LDARotation
export MulticlassRotation
export MeansRotation
export Means2Rotation
export FormulaRotation
export Formula2Rotation
export ThematicRotation
export UMAPRotation
export ena_dataset
export derivedAnyCode!
export derivedAllCode!

# @warn "Running EpistemicNetworkAnalysis.jl as main. Performing kitchen sink operation."
# let

# data = ena_dataset("shakespeare.data")
# conversations = [:Play, :Act]
# units = [:Play, :Speaker]
# codes = [
#     :Love,
#     :Death,
#     :Honor,
#     :Men,
#     :Women,
#     # :Beauty,
#     # :Fear,
#     # :Friendship,
#     # :Hate,
#     # :Pride
# ]

# myENA = DirectedENAModel(
#     data, codes, conversations, units,
#     windowSize=4,
#     rotateBy=MeansRotation(:Play, "Romeo and Juliet", "Hamlet"),
#     # fitNodesToCircle=true
# )

# display(myENA)
# p = plot(myENA, flipY=true)
# savefig(p, "~/Downloads/temp.svg")
# display(p)

# end # let


end # module
