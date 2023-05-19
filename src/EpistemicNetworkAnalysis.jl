module EpistemicNetworkAnalysis

# Imports
## Plotting
using Plots
using Plots.PlotMeasures
import Plots.plot
import Plots.plot!
import Plots.Plot
using Colors
using Dierckx
using DataStructures

## Data
using DataFrames
using CSV

## Math
using Statistics
using LinearAlgebra
using MultivariateStats
using HypothesisTests
using GLM

## Nonlinear
# using Random
# using Distances
# using UMAP

# Dataflow Logic
using Pipe

# Helpers
include("./enadevtools.jl")
include("./linearmodeling.jl")
include("./linearplotting.jl")
include("./utils.jl")

# Rotations
include("./rotations/ManualRotation.jl")
include("./rotations/SVDRotation.jl")
include("./rotations/FormulaRotation.jl")
include("./rotations/MeansRotation.jl")
# include("./rotations/LDARotation.jl")
# include("./rotations/MulticlassRotation.jl")
include("./rotations/TopicRotation.jl")

# Models
include("./models/ENAModel.jl")
# include("./models/DigraphENAModel.jl")
include("./models/BiplotENAModel.jl")

# Exports
export ENAModel
export BiplotENAModel
# export DigraphENAModel
# export NonlinearENAModel
export plot
export SVDRotation
export ManualRotation
# export LDARotation
# export MulticlassRotation
export MeansRotation
export FormulaRotation
export TopicRotation
# export UMAPRotation
export loadExample
# export derivedAnyCode!
# export derivedAllCode!

end # module
