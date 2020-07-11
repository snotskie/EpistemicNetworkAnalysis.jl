"""
TODO document
"""
struct MeansRotation <: AbstractMeansRotation
    regression_model::Type{LinearModel}
    coefindex::Int
    f1::FormulaTerm
    contrasts::Union{Nothing,Dict}
    groupVar::Symbol
    controlGroup::Any
    treatmentGroup::Any
end

# Simplified constructor
function MeansRotation(groupVar::Symbol, controlGroup::Any, treatmentGroup::Any)
    regression_model = LinearModel
    coefindex = 2
    f1 = @formula(y ~ 1 + FactoredGroupVar)
    contrasts = nothing
    return MeansRotation(regression_model, coefindex, f1, contrasts, groupVar, controlGroup, treatmentGroup)
end

# Implement rotate
function rotate!(rotation::AbstractMeansRotation, networkModel::DataFrame, centroidModel::DataFrame)

    ## Manually factor the grouping variable to 0/1/missing
    centroidModel[!, :FactoredGroupVar] = map(eachrow(centroidModel)) do unitRow
        if unitRow[rotation.groupVar] == rotation.controlGroup
            return 0.0
        elseif unitRow[rotation.groupVar] == rotation.treatmentGroup
            return 1.0
        else
            return missing
        end
    end

    ## Use a FormulaRotation to do the rest of the work
    invoke(rotate!, Tuple{AbstractFormulaRotation, DataFrame, DataFrame}, rotation, networkModel, centroidModel)
end

# Override plotting pieces
# ## Units - we can color them into two groups
# function plot_units!(p::Plot, ena::AbstractENAModel{<:AbstractMeansRotation}, displayCentroids::DataFrame, displayCounts::DataFrame;
#     flipX::Bool=false, flipY::Bool=false,
#     kwargs...)


#     unitColors = map(eachrow(displayCentroids)) do unitRow
#         if unitRow[ena.rotation.groupVar] == ena.rotation.controlGroup
#             return :purple
#         elseif unitRow[ena.rotation.groupVar] == ena.rotation.treatmentGroup
#             return :orange
#         else
#             return :black
#         end
#     end

#     x = displayCentroids[!, :pos_x] * (flipX ? -1 : 1)
#     y = displayCentroids[!, :pos_y] * (flipY ? -1 : 1)
#     plot!(p, x, y,
#         seriestype=:scatter,
#         markershape=:circle,
#         markersize=1.5,
#         markercolor=unitColors,
#         markerstrokecolor=unitColors)
# end

## CIs - we can color them into two groups
function plot_intervals!(p::Plot, ena::AbstractENAModel{<:AbstractMeansRotation}, displayCentroids::DataFrame, displayCounts::DataFrame;
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    controlUnits = filter(displayCentroids) do unitRow
        return unitRow[ena.rotation.groupVar] == ena.rotation.controlGroup
    end

    xs = controlUnits[!, :pos_x] * (flipX ? -1 : 1)
    ys = controlUnits[!, :pos_y] * (flipY ? -1 : 1)
    help_plot_ci(p, xs, ys, :purple, :square)

    treatmentUnits = filter(displayCentroids) do unitRow
        return unitRow[ena.rotation.groupVar] == ena.rotation.treatmentGroup
    end

    xs = treatmentUnits[!, :pos_x] * (flipX ? -1 : 1)
    ys = treatmentUnits[!, :pos_y] * (flipY ? -1 : 1)
    help_plot_ci(p, xs, ys, :orange, :square)
end

## Network - we can do a subtraction plot
function plot_network!(p::Plot, ena::AbstractENAModel{<:AbstractMeansRotation}, displayCentroids::DataFrame, displayCounts::DataFrame;
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    controlRows = map(x->x[ena.rotation.groupVar] == ena.rotation.controlGroup, eachrow(displayCentroids))
    treatmentRows = map(x->x[ena.rotation.groupVar] == ena.rotation.treatmentGroup, eachrow(displayCentroids))

    lineWidths = map(eachrow(network(ena))) do networkRow
        controlMean = mean(displayCounts[controlRows, networkRow[:relationship]])
        treatmentMean = mean(displayCounts[treatmentRows, networkRow[:relationship]])
        if !isnan(controlMean) && !isnan(treatmentMean)
            return treatmentMean - controlMean
        elseif !isnan(controlMean)
            return -controlMean
        elseif !isnan(treatmentMean)
            return treatmentMean
        else
            return 0
        end
    end

    lineColors = map(lineWidths) do width
        if width < 0
            return :purple
        else
            return :orange
        end
    end

    lineWidths = abs.(lineWidths)
    lineWidths *= 2 / maximum(lineWidths)
    codeWidths = zeros(nrow(nodes(ena)))
    for (i, networkRow) in enumerate(eachrow(network(ena)))
        j, k = relationships(ena)[networkRow[:relationship]]
        codeWidths[j] += lineWidths[i]
        codeWidths[k] += lineWidths[i]

        x = nodes(ena)[[j, k], :pos_x] * (flipX ? -1 : 1)
        y = nodes(ena)[[j, k], :pos_y] * (flipY ? -1 : 1)
        plot!(p, x, y,
            seriestype=:line,
            linewidth=lineWidths[i],
            linecolor=lineColors[i])
    end

    codeWidths *= 8 / maximum(codeWidths)
    x = nodes(ena)[!, :pos_x] * (flipX ? -1 : 1)
    y = nodes(ena)[!, :pos_y] * (flipY ? -1 : 1)
    labels = map(label->text(label, :top, 8), nodes(ena)[!, :code])
    plot!(p, x, y,
        seriestype=:scatter,
        series_annotations=labels,
        markershape=:circle,
        markersize=codeWidths,
        markercolor=:black,
        markerstrokecolor=:black)
end