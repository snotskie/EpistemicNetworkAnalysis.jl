module EpistemicNetworkAnalysis

# Imports
## Plotting
using Plots

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

# Debugging behavior for when run as main
if abspath(PROGRAM_FILE) == @__FILE__
    @warn "Running EpistemicNetworkAnalysis.jl as main."
    RSdata = ena_dataset("RS.data")
    codes = [:Data,
        :Technical_Constraints,
        :Performance_Parameters,
        :Client_and_Consultant_Requests,
        :Design_Reasoning,
        :Collaboration]

    conversations = [:Condition, :GameHalf, :GroupName]
    units = [:Condition, :GameHalf, :UserName]
    myRotation = TwoGroupRotation(:Condition, "FirstGame", "SecondGame",
                                    :GameHalf, "First", "Second",
                                    [])

    myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=myRotation)
    display(myENA)

    # myArtist = MeansArtist(:Condition, "FirstGame", "SecondGame")
    # myArtist = MeansArtist(:GameHalf, "First", "Second")
    # myArtist = WindowsArtist(:Condition, "FirstGame", "SecondGame",
    #                          :GameHalf, "First", "Second")
    myArtist = TVRemoteArtist(:Condition, "FirstGame", "SecondGame",
                                :GameHalf, "First", "Second")
    p = plot(myENA,
        showprojection=true,
        artist=myArtist
    )

    display(p)
end

# Exports
export ENAModel
export ENARotation, SVDRotation, MeansRotation
export ENAArtist, DefaultArtist, MeansArtist
export plot
export ena_dataset

end # module