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

# # Includes
# include("./helpers.jl")
# include("./ENARotations.jl")
# include("./ENAArtists.jl")
# include("./ENAModel.jl")
# include("./ENADisplay.jl")
# include("./RSData.jl")

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

# # rotation = Formula2Rotation(
# #     LinearModel, 2, @formula(y ~ 1 + Condition + GameHalf), Dict(:Condition => EffectsCoding(), :GameHalf => EffectsCoding()),
# #     LinearModel, 2, @formula(y ~ 1 + GameHalf + Condition), Dict(:Condition => EffectsCoding(), :GameHalf => EffectsCoding())
# # )



# # rotation = Formula2Rotation(
# #     LinearModel, 2, @formula(y ~ 1 + CONFIDENCE_Pre + GameHalf), Dict(:GameHalf => EffectsCoding()),
# #     LinearModel, 2, @formula(y ~ 1 + GameHalf + CONFIDENCE_Pre), Dict(:GameHalf => EffectsCoding())
# # )

# RSdata = ena_dataset("RS.data")
# RSdata[!, :RND] = rand(nrow(RSdata))
# codes = [:Data,
#     :Technical_Constraints,
#     :Performance_Parameters,
#     :Client_and_Consultant_Requests,
#     :Design_Reasoning,
#     :Collaboration]

# conversations = [:Condition, :GameHalf, :GroupName]
# units = [:Condition, :GameHalf, :UserName]

# rotation = Formula2Rotation(
#     LinearModel, 2, @formula(y ~ 1 + RND + GameHalf), Dict(:GameHalf => EffectsCoding()),
#     LinearModel, 2, @formula(y ~ 1 + GameHalf + RND), Dict(:GameHalf => EffectsCoding())
# )

# myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation)
# display(myENA)

# p = plot(myENA, title="test", xlabel="RND", ylabel="GameHalf")
# display(p)



end # module