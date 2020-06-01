module EpistemicNetworkAnalysis

using DataFrames
using Statistics
using LinearAlgebra
using MultivariateStats
using HypothesisTests
using CSV
# TODO using Lasso
using ImageMagick # fix for cairo bug, have this first
using CairoMakie # using CairoMakie instead of GLMakie because the latter requires a gpu
using Makie

include("./helpers.jl")
include("./ENARotations.jl")
include("./ENAArtists.jl")
include("./ENAModel.jl")
include("./ENADisplay.jl")
include("./RSData.jl")
include("./examples.jl") # TEMP

# TODO exports

temp_example() # TEMP

end # module