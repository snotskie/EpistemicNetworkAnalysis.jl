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
# TODO using Lasso

# Includes
include("./helpers.jl")
include("./typetree.jl")
include("./plotbase.jl")
include("./SVDRotation.jl")
include("./FormulaRotation.jl")
include("./Formula2Rotation.jl")
include("./MeansRotation.jl")
include("./Means2Rotation.jl")
include("./ENAModel.jl")
include("./RSData.jl")

# Exports
export ENAModel
export plot
export MeansRotation
export Means2Rotation
export SVDRotation
export FormulaRotation
export Formula2Rotation
export ena_dataset


# @warn "Running EpistemicNetworkAnalysis.jl as main. Performing kitchen sink operation."
# RSdata = ena_dataset("RS.data")
# using Random
# Random.seed!(4321)
# RSdata[!, :RND] = rand(nrow(RSdata))
# println(names(RSdata))
# codes = [:Data,
#     :Technical_Constraints,
#     :Performance_Parameters,
#     :Client_and_Consultant_Requests,
#     :Design_Reasoning,
#     :Collaboration]

# conversations = [:Condition, :GameHalf, :GroupName]
# units = [:Condition, :GameHalf, :UserName]

# # myENA = ENAModel(RSdata, codes, conversations, units)

# # rotation = SVDRotation()

# # rotation = MeansRotation(:Condition, "FirstGame", "SecondGame")

# rotation = Means2Rotation(
#     :Condition, "FirstGame", "SecondGame",
#     :GameHalf, "First", "Second"
# )

# # rotation = FormulaRotation(
# #     LinearModel, 2, @formula(y ~ 1 + CONFIDENCE_Pre), nothing
# # )

# # rotation = FormulaRotation(
# #     LinearModel, 2, @formula(y ~ 1 + RND), nothing
# # )

# # rotation = Formula2Rotation(
# #     LinearModel, 2, @formula(y ~ 1 + CONFIDENCE_Pre), nothing,
# #     LinearModel, 1, @formula(y ~ 1 + 0), nothing
# #     # LinearModel, 1, @formula(y ~ 1 + 0), nothing,
# #     # LinearModel, 2, @formula(y ~ 1 + CONFIDENCE_Pre), nothing
# # )

# myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation, dropEmpty=false, deflateEmpty=false, meanCenter=true, sphereNormalize=true)
# display(myENA)
# # p = plot(myENA, showExtras=false, showUnits=true, showNetworks=true, showCIs=true, flipY=false, leg=true)
# # p = plot(myENA, groupBy=:GroupName)
# p = plot(myENA)
# display(p)


end # module