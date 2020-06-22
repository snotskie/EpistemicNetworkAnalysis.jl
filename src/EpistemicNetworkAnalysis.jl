module EpistemicNetworkAnalysis

# Imports
## Plotting
using Plots
using Plots.PlotMeasures

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

    myRotations = Dict{String,ENARotation}()
    myRotations["svd"] = SVDRotation()
    myRotations["means-condition"] = MeansRotation(:Condition, "FirstGame", "SecondGame")
    myRotations["formula2-repeat"] =  Formula2Rotation(
        LinearModel, @formula(y ~ 1 + FactoredCondition + FactoredGameHalf + FactoredGameHalf&FactoredCondition),
        LinearModel, @formula(y ~ 1 + FactoredGameHalf + FactoredCondition + FactoredGameHalf&FactoredCondition)
    )

    # myRotations["formula2-norepeat"] =  Formula2Rotation(
    #     LinearModel, @formula(y ~ 1 + FactoredCondition + FactoredGameHalf + FactoredGameHalf&FactoredCondition),
    #     LinearModel, @formula(y ~ 1 + FactoredGameHalf + FactoredGameHalf&FactoredCondition)
    # )

    myArtists = Dict{String,ENAArtist}()
    # myArtists["black"] = DefaultArtist()
    myArtists["color-condition"] = MeansArtist(:Condition, "FirstGame", "SecondGame")
    # myArtists["color-gamehalf"] = MeansArtist(:GameHalf, "First", "Second")
    myArtists["windows"] = WindowsArtist(
        :Condition, "FirstGame", "SecondGame",
        :GameHalf, "First", "Second"
    )

    myArtists["tv-condition-x-gamehalf-y"] = TVRemoteArtist(
        :Condition, "FirstGame", "SecondGame",
        :GameHalf, "First", "Second"
    )

    counter = 0
    for (label1, rotation) in myRotations
        println(typeof(rotation))
        myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=rotation)
        display(myENA)
        for (label2, artist) in myArtists
            println(typeof(artist))
            p1 = plot(myENA,
                artist=artist
            )

            p2 = plot(myENA,
                showunits=false,
                showconfidence=true,
                artist=artist
            )

            global counter
            savefig(p1, "kitchensink-units-$(counter)-$(label1)-$(label2).svg")
            savefig(p2, "kitchensink-means-$(counter)-$(label1)-$(label2).svg")
            counter += 1
        end
    end
end

end # module