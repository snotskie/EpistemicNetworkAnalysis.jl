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
# # data = ena_dataset("RS.data")
# # data[!, :_] = ones(nrow(data))
# # data[!, :_1] = zeros(nrow(data))
# # data[!, :_2] = zeros(nrow(data))
# # data[!, :_3] = [row[:Condition] == "FirstGame" ? 1 : 0 for row in eachrow(data)]
# # codes = [
# #     :Data,
# #     :Technical_Constraints,
# #     :Performance_Parameters,
# #     :Client_and_Consultant_Requests,
# #     :Design_Reasoning,
# #     :Collaboration,
# #     # :_1,
# #     # :_2,
# #     # :_3,
# #     # :_
# # ]

# # conversations = [:Condition, :GameHalf, :GroupName]
# # units = [:Condition, :GameHalf, :UserName]

# # Data
# data = ena_dataset("shakespeare.data")

# # Config
# codes = [
#     :Love,
#     :Beauty,
#     :Death,
#     :Fear,
#     :Friendship,
#     :Hate,
#     :Honor,
#     :Men,
#     :Women,
#     :Pride
# ]

# conversations = [:Play, :Act, :Scene]
# units = [:Play, :Act, :Speaker]

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

# data[!, :RND3Group] = map(data[!, :RND]) do x
#     if x <= 1/3
#         return "First"
#     elseif x <= 2/3
#         return "Second"
#     else
#         return "Third"
#     end
# end

# # rotation = SVDRotation()
# # rotation = SVDRotation(5)
# # rotation = LDARotation(:Play)
# # rotation = LDARotation(:GameHalf)
# # rotation = LDARotation(:RNDGroup)
# # rotation = LDARotation(:RND3Group)
# # rotation = LDARotation(:GroupName)
# # rotation = LDARotation(:GroupName, 2)
# # rotation = MeansRotation(:Condition, "FirstGame", "SecondGame")
# # rotation = MeansRotation(:Collaboration, 0, 1)
# # rotation = MeansRotation(:Condition, "SecondGame", "FirstGame")
# # rotation = MeansRotation(:GameHalf, "First", "Second")
# # rotation = MeansRotation(:RNDGroup, "First", "Second")
# # rotation = Means2Rotation(:Condition, "SecondGame", "FirstGame",
# #                           :GameHalf, "First", "Second")
# # rotation = FormulaRotation(
# #     LinearModel, 2, @formula(col ~ 1 + RND), nothing
# # )
# # using NetworkLayout
# # rotation = CodeNetworkRotation(NetworkLayout.Stress)
# # using Lasso
# # rotation = FormulaRotation(
# #     LassoModel, 2, @formula(col ~ 0 + RND), nothing
# # )

# # someData = data[data[!, :Condition] .== "FirstGame", :]
# # someData = someData[someData[!, :GameHalf] .== "Second", :]
# # rotation = DifferenceRotation(1, 4)
# # rotation = DirectionRotation(6)
# # rotation = ThematicRotation([:Technical_Constraints], [:Performance_Parameters])
# # rotation = ThematicRotation([:Client_and_Consultant_Requests, :Collaboration], [:Design_Reasoning, :Performance_Parameters, :Data])
# # rotation = ThematicRotation([:Men], [:Women])
# # rotation = ThematicRotation([:Honor], [:Women])
# rotation = ThematicRotation([:Men, :Honor], [:Women, :Hate, :Love])
# myENA = ENAModel(
# # myENA = BiplotModel(
#     data, codes, conversations, units,
#     rotateBy=rotation,
#     # subspace = 6,
#     # subsetFilter=(x->x[:RND] < 1),
#     # subsetFilter=(x->x[:GroupName] in ["1", "2", "3"] && x[:Condition] == "FirstGame"),
#     # dimensionNormalize=true,
#     # rotateOn=:accumModel,
#     # rotateOn=:codeModel,
#     rotateOn=:centroidModel,
#     # deflateEmpty=true,
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
# savefig(plot(myENA, showWarps=false, weakLinks=true, groupBy=:Play), "~/Downloads/temp.png")
# end # let


end # module