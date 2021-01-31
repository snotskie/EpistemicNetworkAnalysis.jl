# Rotations
abstract type AbstractENARotation
    # fields: (none)
    # plot accepts: flipX, flipY, xlabel, ylabel
    # test reports: variance_x, variance_y
end

abstract type AbstractSVDRotation <: AbstractENARotation
    # fields: (inherit)
    # plot accepts: (inherit), groupVar
end

abstract type AbstractFormulaRotation <: AbstractENARotation
    # fields: (inherit), regression_model, coefindex, f1, contrasts
    # plot accepts: (inherit), minLabel, maxLabel, minColor, maxColor
end

abstract type AbstractFormula2Rotation <: AbstractFormulaRotation
    # fields: (inherit), regression_model2, coefindex2, f2, contrasts
    # plot accepts: (inherit)
end

abstract type AbstractMeansRotation <: AbstractFormulaRotation
    # fields: (inherit), groupVar, controlGroup, treatmentGroup
    # plot accepts: (inherit)
end

abstract type AbstractMeans2Rotation <: AbstractFormula2Rotation
    # fields: (inherit), groupVar, controlGroup, treatmentGroup, groupVar2, controlGroup2, treatmentGroup2
    # plot accepts: (inherit)
end

# Accumulation Models
abstract type AbstractENAModel{T<:AbstractENARotation}
    # fields: units, conversations, codes, rotation, accumModel, centroidModel, metadata, codeModel, networkModel, relationshipMap
    # test reports: coregistration
end

# Default Functions
## Rotations
function rotate!(rotation::AbstractENARotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)
    error("Unimplemented")
end

## Tests
function test(ena::AbstractENAModel)

    # For reference:
    # p = pvalue(OneSampleTTest(fitDiffs, dimDiffs))
    # pvalue(EqualVarianceTTest(x, y))
    # pvalue(UnequalVarianceTTest(x, y))
    # pvalue(MannWhitneyUTest(x, y))
    # pvalue(SignedRankTest(x, y))

    ### Find difference between each pair of points, in accum and centroid
    centroidDiffsX = Real[]
    centroidDiffsY = Real[]
    accumDiffsX = Real[]
    accumDiffsY = Real[]
    for (i, unitRowA) in enumerate(eachrow(ena.accumModel))
        for (j, unitRowB) in enumerate(eachrow(ena.accumModel))
            if i < j
                push!(centroidDiffsX, ena.centroidModel[i, :pos_x] - ena.centroidModel[j, :pos_x])
                push!(centroidDiffsY, ena.centroidModel[i, :pos_y] - ena.centroidModel[j, :pos_y])
                push!(accumDiffsX, unitRowA[:pos_x] - unitRowB[:pos_x])
                push!(accumDiffsY, unitRowA[:pos_y] - unitRowB[:pos_y])
            end
        end
    end

    corr_x = cor(ena.centroidModel[!, :pos_x], ena.accumModel[!, :pos_x])
    corr_y = cor(ena.centroidModel[!, :pos_y], ena.accumModel[!, :pos_y])

    ### Do those differences correlate?
    pearson_x = cor(centroidDiffsX, accumDiffsX)
    pearson_y = cor(centroidDiffsY, accumDiffsY)

    ### Find the percent variance explained by the x and y axis of the entire high dimensional space
    total_variance = sum(var.(eachcol(ena.centroidModel[!, ena.networkModel[!, :relationship]])))
    variance_x = var(ena.centroidModel[!, :pos_x]) / total_variance
    variance_y = var(ena.centroidModel[!, :pos_y]) / total_variance

    ### Package and return
    return Dict(
        :correlation_x => corr_x,
        :correlation_y => corr_y,
        :coregistration_x => pearson_x,
        :coregistration_y => pearson_y,
        :variance_x => variance_x,
        :variance_y => variance_y
    )
end

## Text display
function Base.display(ena::AbstractENAModel) # TODO should this be print, display, or show?

    ### Show centroids
    println("Units (centroids):")
    show(ena.centroidModel, allrows=true)
    println()

    ### Show codes
    println("Codes:")
    show(ena.codeModel, allrows=true)
    println()

    ### Show network
    println("Network:")
    show(ena.networkModel, allrows=true)
    println()

    ### Show every test result we have
    results = test(ena)
    for key in keys(results)
        println("$key:")
        println(results[key])
        println()
    end
end