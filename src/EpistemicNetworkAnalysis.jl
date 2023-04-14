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

# Helpers
include("./enadevtools.jl")
include("./lineardefaults.jl")
include("./utils.jl")

# Rotations
include("./rotations/SVDRotation.jl")
include("./rotations/ModeratedRotation.jl")
# include("./rotations/LDARotation.jl")
# include("./rotations/MulticlassRotation.jl")
# include("./rotations/FormulaRotation.jl")
# include("./rotations/Formula2Rotation.jl")
# include("./rotations/MeansRotation.jl")
# include("./rotations/Means2Rotation.jl")
# include("./rotations/ThematicRotation.jl")

# Models
include("./models/ENAModel.jl")
# include("./models/DigraphENAModel.jl")
include("./models/BiplotENAModel.jl")

# Exports
export ENAModel
export BiplotENAModel
# export DigraphENAModel
# export NonlinearENAModel
# export plot
export SVDRotation
export ModeratedRotation
# export LDARotation
# export MulticlassRotation
# export MeansRotation
# export Means2Rotation
# export FormulaRotation
# export Formula2Rotation
# export ThematicRotation
# export UMAPRotation
export loadExample
# export derivedAnyCode!
# export derivedAllCode!

end # module
