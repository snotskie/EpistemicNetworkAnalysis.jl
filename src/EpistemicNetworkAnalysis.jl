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
using StatsBase
using Statistics
using LinearAlgebra
using MultivariateStats
using HypothesisTests
using GLM

## Printing
using PkgVersion
using NamedTupleTools
using PrettyPrinting
using Tables
import Base.show

## Serialization
using XLSX
using Dates

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
include("./lineartesting.jl")
include("./utils.jl")

# Rotations
include("./rotations/TrainedRotation.jl")
include("./rotations/SVDRotation.jl")
include("./rotations/FormulaRotation.jl")
include("./rotations/MeansRotation.jl")
include("./rotations/AbstractGroupDifferenceRotation.jl")
include("./rotations/LDARotation.jl")
include("./rotations/MulticlassRotation.jl")
include("./rotations/TopicRotation.jl")

# Models
include("./models/ENAModel.jl")
include("./models/DigraphENAModel.jl")
include("./models/BiplotENAModel.jl")

# Exports
export ENAModel
export BiplotENAModel
export DigraphENAModel
# export NonlinearENAModel
export plot
export SVDRotation
export TrainedRotation
export LDARotation
export MulticlassRotation
export MeansRotation
export FormulaRotation
export TopicRotation
# export UMAPRotation
export loadExample
export to_xlsx
export deriveAnyCode!
export deriveAllCode!
export statistics
export pointcloud

end # module
