module EpistemicNetworkAnalysis

# Imports
## Plotting
using Plots
using Plots.PlotMeasures
import Plots.plot
import Plots.plot!
import Plots.Plot

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
# RSdata = ena_dataset("RS.data")
# codes = [:Data,
#     :Technical_Constraints,
#     :Performance_Parameters,
#     :Client_and_Consultant_Requests,
#     :Design_Reasoning,
#     :Collaboration]

# conversations = [:Condition, :GameHalf, :GroupName]
# units = [:Condition, :GameHalf, :UserName]
# # rotation = SVDRotation()
# rotation = Formula2Rotation(
#     LinearModel, 2, @formula(y ~ 1 + Condition), Dict(:Condition => DummyCoding()),
#     LinearModel, 2, @formula(y ~ 1 + GameHalf), Dict(:GameHalf => DummyCoding())
# )

# myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation)
# display(myENA)
# p = plot(myENA, title="test", flipX=true, groupVar=:Condition)
# display(p)



end # module