# Rotations
abstract type AbstractENARotation
    # fields: (none)
    # plot accepts: flipX, flipY, xlabel, ylabel
    # test reports: variance_x, variance_y
end

abstract type AbstractSVDRotation <: AbstractENARotation
    # fields: (inherit)
    # plot accepts: (inherit), groupVar
    # test reports: (inherit)
end

abstract type AbstractFormulaRotation <: AbstractENARotation
    # fields: (inherit), regression_model, coefindex, f1, contrasts
    # plot accepts: (inherit), minLabel, maxLabel, minColor, maxColor
    # test reports: (inherit), pvalue, effect_size
end

abstract type AbstractFormula2Rotation <: AbstractFormulaRotation
    # fields: (inherit), regression_model2, coefindex2, f2, contrasts
    # plot accepts: (inherit)
    # test reports: (inherit), pvalue2, effect_size2
end

abstract type AbstractMeansRotation <: AbstractFormulaRotation
    # fields: (inherit), groupVar, controlGroup, treatmentGroup
    # plot accepts: (inherit)
    # test reports: (inherit)
end

abstract type AbstractDoubleMeansRotation <: AbstractFormula2Rotation
    # fields: (inherit), groupVar, controlGroup, treatmentGroup, groupVar2, controlGroup2, treatmentGroup2
    # plot accepts: (inherit)
    # test reports: (inherit)
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


# ## Plotting
# ### Top-level wrapper - Do not override
# function plot(ena::AbstractENAModel;
#     margin=10mm, size=500, lims=1, ticks=[-1, 0, 1], title="", leg=:bottomleft,
#     kwargs...)

#     #### Create and empty plot
#     p = plot(leg=leg, margin=margin, size=(size, size))

#     #### Call mutating wrapper
#     plot!(p, ena; kwargs...)

#     #### Set ticks, limits, and title
#     xticks!(p, ticks)
#     yticks!(p, ticks)
#     xlims!(p, -lims, lims)
#     ylims!(p, -lims, lims)
#     title!(p, title)

#     #### Done, return so user can modify from there
#     return p
# end

# ### Mutating wrapper - Do not override
# function plot!(p::Plot, ena::AbstractENAModel;
#     showUnits::Bool=true, showNetwork::Bool=true, showIntervals::Bool=true, showExtras::Bool=true,
#     unitFilter=x->true,
#     kwargs...)

#     #### Run the filter just this once, pass to helpers
#     displayRows = map(unitFilter, eachrow(ena.metadata))

#     #### Call each helper, unless the user asked us not to
#     if showUnits
#         plot_units!(p, ena, displayRows; kwargs...)
#     end

#     if showNetwork
#         plot_network!(p, ena, displayRows; kwargs...)
#     end

#     if showIntervals
#         plot_intervals!(p, ena, displayRows; kwargs...)
#     end

#     if showExtras
#         plot_extras!(p, ena, displayRows; kwargs...)
#     end

#     #### Always call the helper to label the axes
#     plot_labels!(p, ena; kwargs...)
# end

# ### Unit-level helper
# function plot_units!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
#     flipX::Bool=false, flipY::Bool=false,
#     kwargs...)

#     #### Get the x/y positions
#     displayCentroids = ena.centroidModel[displayRows, :]
#     x = displayCentroids[!, :pos_x] * (flipX ? -1 : 1)
#     y = displayCentroids[!, :pos_y] * (flipY ? -1 : 1)

#     #### Draw them in black
#     plot!(p, x, y,
#         label="Units",
#         seriestype=:scatter,
#         markershape=:circle,
#         markersize=2,
#         markercolor=:black,
#         markerstrokecolor=:black)
# end

# ### Network-level helper
# function plot_network!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
#     flipX::Bool=false, flipY::Bool=false,
#     kwargs...)

#     #### Find the true weight on each line
#     displayAccums = ena.accumModel[displayRows, :]
#     lineWidths = map(eachrow(ena.networkModel)) do networkRow
#         return sum(displayAccums[!, networkRow[:relationship]])
#     end

#     #### Rescale the lines
#     lineWidths *= 2 / maximum(lineWidths)

#     #### Initialize code widths, compute while we visit each line
#     codeWidths = zeros(nrow(ena.codeModel))

#     #### For each line...
#     for (i, networkRow) in enumerate(eachrow(ena.networkModel))
#         j, k = ena.relationshipMap[networkRow[:relationship]]

#         #### ...add to its code weights
#         codeWidths[j] += lineWidths[i]
#         codeWidths[k] += lineWidths[i]

#         #### and plot that line
#         x = ena.codeModel[[j, k], :pos_x] * (flipX ? -1 : 1)
#         y = ena.codeModel[[j, k], :pos_y] * (flipY ? -1 : 1)
#         plot!(p, x, y,
#             label=nothing,
#             seriestype=:line,
#             linewidth=lineWidths[i],
#             linecolor=:black)
#     end

#     #### Rescale and draw the codes
#     codeWidths *= 8 / maximum(codeWidths)
#     x = ena.codeModel[!, :pos_x] * (flipX ? -1 : 1)
#     y = ena.codeModel[!, :pos_y] * (flipY ? -1 : 1)
#     labels = map(label->text(label, :top, 8), ena.codeModel[!, :code])
#     plot!(p, x, y,
#         label=nothing,
#         seriestype=:scatter,
#         series_annotations=labels,
#         markershape=:circle,
#         markersize=codeWidths,
#         markercolor=:black,
#         markerstrokecolor=:black)
# end

# ### CI-level helper
# function plot_intervals!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
#     flipX::Bool=false, flipY::Bool=false,
#     kwargs...)
#     #### do nothing
# end

# ### Extras helper
# function plot_extras!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
#     flipX::Bool=false, flipY::Bool=false,
#     kwargs...)
#     #### do nothing
# end

# ### Labels helper
# function plot_labels!(p::Plot, ena::AbstractENAModel;
#     xlabel="X", ylabel="Y",
#     kwargs...)

#     #### Run tests and report the variances explained in the axis labels
#     results = test(ena)
#     xlabel!(p, "$xlabel ($(round(Int, results[:variance_x]*100))%)")
#     ylabel!(p, "$ylabel ($(round(Int, results[:variance_y]*100))%)")
# end