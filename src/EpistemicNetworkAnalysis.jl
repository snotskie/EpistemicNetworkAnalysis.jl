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
include("./DirectedENAModel.jl")
# include("./CopyRotation.jl")
include("./ena_dataset.jl")

# Exports
export ENAModel
export DirectedENAModel
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

# data = ena_dataset("shakespeare.data")
# conversations = [:Play, :Act]
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

# myENA = DirectedENAModel(
#     data, codes, conversations, units,
#     windowSize=4,
#     rotateBy=MeansRotation(:Play, "Romeo and Juliet", "Hamlet"),
#     # fitNodesToCircle=true
# )

# display(myENA)
# p = plot(myENA, flipY=true)
# savefig(p, "~/Downloads/temp.svg")
# display(p)

# x_men = Matrix{Float64}(myENA.accumModel[myENA.accumModel[!, :pos_x] .< 0, [:Death_Men, :Honor_Men, :Women_Men, :Love_Men]])
# men_x = Matrix{Float64}(myENA.accumModel[myENA.accumModel[!, :pos_x] .< 0, [:Men_Death, :Men_Honor, :Men_Women, :Men_Love]])
# display(x_men - men_x)
# println(mean(x_men - men_x))
# end # let


end # module