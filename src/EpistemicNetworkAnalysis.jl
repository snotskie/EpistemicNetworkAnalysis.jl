module EpistemicNetworkAnalysis

# Imports
## Plotting
using ImageMagick # fix for cairo bug, have to have this first
using CairoMakie # fallback for when no gpu present
using Makie

## Data
using DataFrames
using CSV

## Math
using Statistics
using LinearAlgebra
using MultivariateStats
using HypothesisTests
# TODO using Lasso

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
export plot
export ena_dataset

end # module