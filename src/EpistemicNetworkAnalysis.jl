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
include("./ENAModel.jl")
include("./RSData.jl")

# Exports
export ENAModel
export plot
export SVDRotation
export LDARotation
export MeansRotation
export Means2Rotation
export FormulaRotation
export Formula2Rotation
export ena_dataset


# @warn "Running EpistemicNetworkAnalysis.jl as main. Performing kitchen sink operation."
# let
# data = ena_dataset("RS.data")
# codes = [
#     :Data,
#     :Technical_Constraints,
#     :Performance_Parameters,
#     :Client_and_Consultant_Requests,
#     :Design_Reasoning,
#     :Collaboration
# ]

# conversations = [:Condition, :GameHalf, :GroupName]
# units = [:Condition, :GameHalf, :UserName]

# using Random
# Random.seed!(4321)
# data[!, :RND] = rand(nrow(data))
# data[!, :RNDGroup] = map(data[!, :RND]) do x
#     if x <= 0.5
#         return "First"
#     else
#         return "Second"
#     end
# end

# # rotation = SVDRotation()
# # rotation = LDARotation(:Condition)
# # rotation = LDARotation(:GameHalf)
# # rotation = LDARotation(:RNDGroup)
# # rotation = LDARotation(:GroupName)
# # rotation = MeansRotation(:Condition, "FirstGame", "SecondGame")
# # rotation = MeansRotation(:Condition, "SecondGame", "FirstGame")
# # rotation = MeansRotation(:GameHalf, "First", "Second")
# # rotation = MeansRotation(:RNDGroup, "First", "Second")
# # using Lasso
# # rotation = FormulaRotation(
# #     LassoModel, 2, @formula(col ~ 0 + RND), nothing
# # )
# rotation = FormulaRotation(
#     LinearModel, 2, @formula(col ~ 1 + RND), nothing
# )
# myENA = ENAModel(data, codes, conversations, units, rotateBy=rotation, deflateEmpty=true)
# # myENA = ENAModel(data, codes, conversations, units, rotateBy=rotation, subsetFilter=x->x[:GameHalf]=="First"&&x[:Condition]=="FirstGame")
# display(myENA)
# savefig(plot(myENA), "~/Downloads/temp.png")
# # savefig(plot(myENA, groupBy=:Condition), "~/Downloads/temp.png")
# run(`firefox "~/Downloads/temp.png"`)
# end # let


end # module