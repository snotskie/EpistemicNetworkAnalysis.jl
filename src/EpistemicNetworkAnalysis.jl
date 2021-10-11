module EpistemicNetworkAnalysis

# Imports
## Plotting
using Plots
using Plots.PlotMeasures
import Plots.plot
import Plots.plot!
import Plots.Plot
using Colors

## Data
using DataFrames
using CSV

## Math
using Statistics
using LinearAlgebra
using MultivariateStats
using HypothesisTests
using GLM

# Includes
include("./helpers.jl")
include("./typetree.jl")
include("./plotbase.jl")
include("./SVDRotation.jl")
include("./LDARotation.jl")
include("./FormulaRotation.jl")
include("./Formula2Rotation.jl")
include("./MeansRotation.jl")
include("./Means2Rotation.jl")
include("./ThematicRotation.jl")
include("./ENAModel.jl")
include("./BiplotModel.jl")
# include("./CopyRotation.jl")
# include("./CodeNetworkRotation.jl")
# include("./DifferenceRotation.jl")
# include("./DirectionRotation.jl")
include("./ena_dataset.jl")

# Exports
export ENAModel
export BiplotModel
export plot
export SVDRotation
export LDARotation
export MeansRotation
export Means2Rotation
export FormulaRotation
export Formula2Rotation
export ThematicRotation
export ena_dataset

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

# # data = ena_dataset("RS.data")
# # conversations = [:Condition, :GameHalf, :GroupName]
# # units = [:Condition, :GameHalf, :UserName]
# # codes = [
# #     :Data,
# #     :Technical_Constraints,
# #     :Performance_Parameters,
# #     :Client_and_Consultant_Requests,
# #     :Design_Reasoning,
# #     :Collaboration,
# # ]

# myENA = ENAModel(
#     data, codes, conversations, units,
#     windowSize=4,
#     rotateBy=MeansRotation(:Play, "Romeo and Juliet", "Hamlet"),
#     recenterEmpty=true,
#     # dropEmpty=true,
#     # meanCenter=false
# )

# display(myENA)
# savefig(plot(myENA, flipY=true, showCodeLabels=false), "~/Downloads/temp.png")
# end # let


end # module