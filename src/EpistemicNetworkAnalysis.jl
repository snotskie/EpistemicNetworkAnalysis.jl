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
include("./devtools.jl")
using devtools
include("./defaults.jl")
include("./utils.jl")

# Rotations
include("./SVDRotation.jl")
include("./LDARotation.jl")
include("./MulticlassRotation.jl")
include("./FormulaRotation.jl")
include("./Formula2Rotation.jl")
include("./MeansRotation.jl")
include("./Means2Rotation.jl")
include("./ThematicRotation.jl")

# Models
include("./ENAModel.jl")
include("./DigraphENAModel.jl")
include("./Biplot.jl")

# Example Data
include("./loadExample.jl")

# Exports
export ENAModel
export DigraphENAModel
export NonlinearENAModel
export plot
export SVDRotation
export LDARotation
export MulticlassRotation
export MeansRotation
export Means2Rotation
export FormulaRotation
export Formula2Rotation
export ThematicRotation
export UMAPRotation
export loadExample
export derivedAnyCode!
export derivedAllCode!

end # module
