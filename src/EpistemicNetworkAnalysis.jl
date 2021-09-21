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

# # changes should be adding an elseif at line 167 of ENAModel.jl and a param recenterEmpty::Bool=true, also need to swap the comments at lines 332 and 337
# @warn "Running EpistemicNetworkAnalysis.jl as main. Performing kitchen sink operation."
# let

# data = ena_dataset("shakespeare.data")
# conversations = [:Play, :Act, :Scene]
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
#     rotateBy=MeansRotation(:Play, "Hamlet", "Romeo and Juliet"),
#     recenterEmpty=true,
#     # dropEmpty=true,
#     # meanCenter=false
# )

# display(myENA)
# savefig(plot(myENA, flipX=true, flipY=true), "~/Downloads/temp.png")
# end # let


end # module