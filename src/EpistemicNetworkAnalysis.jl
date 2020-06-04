module EpistemicNetworkAnalysis

# Imports
## Math
using Statistics
using LinearAlgebra
using MultivariateStats
using HypothesisTests
# TODO using Lasso

## Data
using DataFrames
using CSV

## Plotting
using ImageMagick # fix for cairo bug, have to have this first
using CairoMakie # fallback for when no gpu present
using Makie

# Includes
include("./helpers.jl")
include("./ENARotations.jl")
include("./ENAArtists.jl")
include("./ENAModel.jl")
include("./ENADisplay.jl")
include("./RSData.jl")

# include("./examples.jl")
# temp_example()

# Exports
export ENAModel
export ENARotation, SVDRotation, MeansRotation
export ENAArtist, DefaultArtist, MeansArtist
export Makie.plot
export ena_dataset

end # module