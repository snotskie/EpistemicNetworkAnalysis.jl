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
include("./ENAModel.jl")
include("./RSData.jl")

# Exports
export ENAModel
export plot
export MeansRotation
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

# rotation = MeansRotation(:Condition, "FirstGame", "SecondGame")
# myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation)

# # rotation = FormulaRotation(
# #     LinearModel, 2, @formula(y ~ 1 + CONFIDENCE_Pre), nothing
# # )
# # myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation)

# # rotation = FormulaRotation(
# #     LinearModel, 2, @formula(y ~ 1 + RND), nothing
# # )
# # myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation)

# display(myENA)
# p = plot(myENA, groupBy=:Condition)
# # p = plot(myENA, groupBy=:GroupName)
# # p = plot(myENA)
# display(p)


end # module