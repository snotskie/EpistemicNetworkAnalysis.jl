module EpistemicNetworkAnalysis

using DataFrames
using Statistics
using LinearAlgebra
using GLM
using CSV
# TODO using Lasso
using Makie

include("./ENAModel.jl")
include("./RSData.jl")
include("./examples.jl")

temp_example()

end # module
