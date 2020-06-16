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
using GLM
# TODO using Lasso

# Includes
include("./helpers.jl")
include("./ENARotations.jl")
include("./ENAArtists.jl")
include("./ENAModel.jl")
include("./ENADisplay.jl")
include("./RSData.jl")

# Exports
export ENAModel
export ENARotation, SVDRotation, MeansRotation
export ENAArtist, DefaultArtist, MeansArtist
export plot
export ena_dataset

# Debugging behavior for when run as main
if abspath(PROGRAM_FILE) == @__FILE__
    @warn "Running EpistemicNetworkAnalysis.jl as main. Performing kitchen sink operation."
    RSdata = ena_dataset("RS.data")
    RSdata[!, :FactoredCondition] = map(RSdata[!, :Condition]) do rowValue
        if rowValue == "FirstGame"
            return 0
        elseif rowValue == "SecondGame"
            return 1
        else
            return missing
        end
    end

    RSdata[!, :FactoredGameHalf] = map(RSdata[!, :GameHalf]) do rowValue
        if rowValue == "First"
            return 0
        elseif rowValue == "Second"
            return 1
        else
            return missing
        end
    end

    RSdata[!, :FactoredCondition] = RSdata[!, :FactoredCondition] .- mean(RSdata[!, :FactoredCondition])
    RSdata[!, :FactoredGameHalf] = RSdata[!, :FactoredGameHalf] .- mean(RSdata[!, :FactoredGameHalf])

    codes = [:Data,
        :Technical_Constraints,
        :Performance_Parameters,
        :Client_and_Consultant_Requests,
        :Design_Reasoning,
        :Collaboration]

    conversations = [:Condition, :GameHalf, :GroupName]
    units = [:Condition, :GameHalf, :UserName]

    myRotations = ENARotation[]
    push!(myRotations, SVDRotation())
    push!(myRotations, MeansRotation(:Condition, "FirstGame", "SecondGame"))
    push!(myRotations, MeansRotation(:GameHalf, "First", "Second"))
    push!(myRotations, Formula2Rotation(
        LinearModel, @formula(y ~ 1 + FactoredCondition),
        LinearModel, @formula(y ~ 1 + FactoredCondition)
    ))

    push!(myRotations, Formula2Rotation(
        LinearModel, @formula(y ~ 1 + FactoredCondition + FactoredGameHalf + FactoredGameHalf&FactoredCondition),
        LinearModel, @formula(y ~ 1 + FactoredGameHalf + FactoredCondition + FactoredGameHalf&FactoredCondition)
    ))

    push!(myRotations, Formula2Rotation(
        LinearModel, @formula(y ~ 1 + FactoredGameHalf + FactoredCondition + FactoredGameHalf&FactoredCondition),
        LinearModel, @formula(y ~ 1 + FactoredCondition + FactoredGameHalf&FactoredCondition)
    ))

    myArtists = ENAArtist[]
    push!(myArtists, DefaultArtist())
    push!(myArtists, MeansArtist(:Condition, "FirstGame", "SecondGame"))
    push!(myArtists, MeansArtist(:GameHalf, "First", "Second"))
    push!(myArtists, WindowsArtist(
        :Condition, "FirstGame", "SecondGame",
        :GameHalf, "First", "Second"
    ))

    push!(myArtists, TVRemoteArtist(
        :Condition, "FirstGame", "SecondGame",
        :GameHalf, "First", "Second"
    ))

    push!(myArtists, TVRemoteArtist(
        :GameHalf, "First", "Second",
        :Condition, "FirstGame", "SecondGame"
    ))

    counter = 0
    for rotation in myRotations
        println(typeof(rotation))
        myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation)
        display(myENA)
        for artist in myArtists
            println(typeof(artist))
            p1 = plot(myENA,
                showprojection=true,
                artist=artist
            )

            p2 = plot(myENA,
                showunits=false,
                showconfidence=true,
                artist=artist
            )

            global counter
            xticks!(p1, [0])
            yticks!(p1, [0])
            savefig(p1, "kitchen-sink-units-$(counter)-$(typeof(rotation))-$(typeof(artist)).svg")
            savefig(p1, "kitchen-sink-units-$(counter)-$(typeof(rotation))-$(typeof(artist)).png")

            xticks!(p2, [0])
            yticks!(p2, [0])
            savefig(p2, "kitchen-sink-means-$(counter)-$(typeof(rotation))-$(typeof(artist)).svg")
            savefig(p2, "kitchen-sink-means-$(counter)-$(typeof(rotation))-$(typeof(artist)).png")
            counter += 1
        end
    end
end

end # module