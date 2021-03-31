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
include("./CopyRotation.jl")
include("./SVDRotation.jl")
include("./LDARotation.jl")
include("./FormulaRotation.jl")
include("./Formula2Rotation.jl")
include("./MeansRotation.jl")
include("./Means2Rotation.jl")
include("./ENAModel.jl")
include("./BiplotModel.jl")
include("./RSData.jl")

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
export ena_dataset


# @warn "Running EpistemicNetworkAnalysis.jl as main. Performing kitchen sink operation."
# let
# data = ena_dataset("RS.data")
# data[!, :_] = ones(nrow(data))
# codes = [
#     :Data,
#     :Technical_Constraints,
#     :Performance_Parameters,
#     :Client_and_Consultant_Requests,
#     :Design_Reasoning,
#     :Collaboration,
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
# # rotation = SVDRotation(1)
# # rotation = LDARotation(:Condition)
# # rotation = LDARotation(:GameHalf)
# # rotation = LDARotation(:RNDGroup)
# # rotation = LDARotation(:GroupName)
# rotation = LDARotation(:GroupName, 2)
# # rotation = MeansRotation(:Condition, "FirstGame", "SecondGame")
# # rotation = MeansRotation(:Condition, "SecondGame", "FirstGame")
# # rotation = MeansRotation(:GameHalf, "First", "Second")
# # rotation = MeansRotation(:RNDGroup, "First", "Second")
# # rotation = Means2Rotation(:Condition, "SecondGame", "FirstGame",
# #                           :GameHalf, "First", "Second")
# # rotation = FormulaRotation(
# #     LinearModel, 2, @formula(col ~ 1 + RND), nothing
# # )
# # using Lasso
# # rotation = FormulaRotation(
# #     LassoModel, 2, @formula(col ~ 0 + RND), nothing
# # )

# # someData = data[data[!, :Condition] .== "FirstGame", :]
# # someData = someData[someData[!, :GameHalf] .== "Second", :]
# myENA = ENAModel(
# # myENA = BiplotModel(
#     data, codes, conversations, units,
#     rotateBy=rotation,
#     # rotateOn=:accumModel,
#     # rotateOn=:codeModel,
#     # rotateOn=:centroidModel,
#     deflateEmpty=true,
#     # meanCenter=false
# )

# # println(std(myENA.accumModel[!, :pos_x]))
# # println(median(myENA.accumModel[!, :pos_x]))
# # println(median(myENA.accumModel[myENA.metadata[!, :Condition] .== "FirstGame", :pos_x]))
# # println(median(myENA.accumModel[myENA.metadata[!, :Condition] .== "SecondGame", :pos_x]))

# # # myENA = ENAModel(
# # # # myENA = BiplotModel(
# # #     data, codes, conversations, units,
# # #     rotateBy=CopyRotation(myENA),
# # #     rotateOn=:accumModel,
# # #     # rotateOn=:codeModel,
# # #     # deflateEmpty=true,
# # #     # meanCenter=false
# # # )

# display(myENA)
# savefig(plot(myENA), "~/Downloads/temp.png")
# end # let


end # module