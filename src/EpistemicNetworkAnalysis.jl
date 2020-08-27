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

# Exports
export ENAModel
export plot
export MeansRotation
export SVDRotation
export FormulaRotation
export Formula2Rotation
export ena_dataset


# @warn "Running EpistemicNetworkAnalysis.jl as main. Performing kitchen sink operation."
# # include("./CarlsFormula2Rotation.jl")

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
# for t in 1:4
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

# ## Bryan's
# rotation = Formula2Rotation(
#     LinearModel, 2, @formula(y ~ 1 + RND + Condition + GameHalf), Dict(:Condition => EffectsCoding(), :GameHalf => EffectsCoding()),
#     LinearModel, 3, @formula(y ~ 1 + RND + Condition + GameHalf), Dict(:Condition => EffectsCoding(), :GameHalf => EffectsCoding())
# )

# myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation)
# display(myENA)

# p = plot(myENA, title="Bryans", ylabel="Condition", xlabel="RND", minColor=colorant"blue", maxColor=colorant"red")#, display_filter=unitRow->unitRow[:Condition]=="FirstGame")
# display(p)

# # ## Carl's
# # rotation = CarlsFormula2Rotation(
# #     LinearModel, 2, @formula(y ~ 1 + Condition + GameHalf), Dict(:Condition => EffectsCoding(), :GameHalf => EffectsCoding()),
# #     LinearModel, 3, @formula(y ~ 1 + Condition + GameHalf), Dict(:Condition => EffectsCoding(), :GameHalf => EffectsCoding())
# # )

# # myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation)
# # display(myENA)

# # p = plot(myENA, title="Carls", ylabel="GameHalf", xlabel="Condition", minColor=colorant"blue", maxColor=colorant"red")#, display_filter=unitRow->unitRow[:Condition]=="FirstGame")
# # display(p)

# # ## MR1
# # rotation = MeansRotation(:Condition, "FirstGame", "SecondGame")

# # myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation, windowSize=4)
# # display(myENA)

# # p = plot(myENA, title="MR1", ylabel="SVD", xlabel="Condition", minColor=colorant"blue", maxColor=colorant"red")#, display_filter=unitRow->unitRow[:Condition]=="FirstGame")
# # display(p)



end # module