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
include("./SVDRotation.jl")
include("./FormulaRotation.jl")
include("./Formula2Rotation.jl")
include("./MeansRotation.jl")
include("./ENAModel.jl")
include("./RSData.jl")

include("./CarlsFormula2Rotation.jl") # TEMP

# # Exports
# export ENAModel
# export ENARotation, SVDRotation, MeansRotation
# export ENAArtist, DefaultArtist, MeansArtist
# export plot
# export ena_dataset



# @warn "Running EpistemicNetworkAnalysis.jl as main. Performing kitchen sink operation."

# # rotation = SVDRotation()
# # rotation = FormulaRotation(
# #     LinearModel, 2, @formula(y ~ 1 + Condition + GameHalf), Dict(:Condition => EffectsCoding(), :GameHalf => EffectsCoding())
# # )

# # rotation = FormulaRotation(
# #     LinearModel, 2, @formula(y ~ 1 + RND), nothing
# # )

# # rotation = MeansRotation(
# #     :Condition, "FirstGame", "SecondGame"
# # )

# rotation = Formula2Rotation(
#     LinearModel, 2, @formula(y ~ 1 + Condition + GameHalf), Dict(:Condition => EffectsCoding(), :GameHalf => EffectsCoding()),
#     LinearModel, 2, @formula(y ~ 1 + GameHalf + Condition), Dict(:Condition => EffectsCoding(), :GameHalf => EffectsCoding())
# )

# # rotation = CarlsFormula2Rotation(
# #     LinearModel, 2, @formula(y ~ 1 + Condition + GameHalf), Dict(:Condition => EffectsCoding(), :GameHalf => EffectsCoding()),
# #     LinearModel, 2, @formula(y ~ 1 + GameHalf + Condition), Dict(:Condition => EffectsCoding(), :GameHalf => EffectsCoding())
# # )

# # rotation = Formula2Rotation(
# #     LinearModel, 2, @formula(y ~ 1 + CONFIDENCE_Pre + GameHalf), Dict(:GameHalf => EffectsCoding()),
# #     LinearModel, 2, @formula(y ~ 1 + GameHalf + CONFIDENCE_Pre), Dict(:GameHalf => EffectsCoding())
# # )

# RSdata = ena_dataset("RS.data")
# using Random
# Random.seed!(4321)
# RSdata[!, :RND] = rand(nrow(RSdata))
# codes = [:Data,
#     :Technical_Constraints,
#     :Performance_Parameters,
#     :Client_and_Consultant_Requests,
#     :Design_Reasoning,
#     :Collaboration]

# conversations = [:Condition, :GameHalf, :GroupName]
# units = [:Condition, :GameHalf, :UserName]
# for t in 1:3
#     RScopy = copy(RSdata)
#     rnd = rand()
#     for copyRow in eachrow(RScopy)
#         for col in [:GroupName, :UserName]
#             copyRow[col] = string(copyRow[col], rnd)
#         end

#         for col in codes
#             if rand() < 0.05
#                 copyRow[col] = 1 - copyRow[col]
#             end
#         end
#     end

#     global RSdata = vcat(RSdata, RScopy)
# end

# # rotation = Formula2Rotation(
# #     LinearModel, 2, @formula(y ~ 1 + RND + GameHalf), Dict(:GameHalf => EffectsCoding()),
# #     LinearModel, 2, @formula(y ~ 1 + GameHalf + RND), Dict(:GameHalf => EffectsCoding())
# # )

# myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation)
# display(myENA)

# p = plot(myENA, title="Bryans", xlabel="Condition", ylabel="GameHalf", minColor=colorant"blue", maxColor=colorant"red")#, display_filter=unitRow->unitRow[:Condition]=="FirstGame")
# display(p)



end # module