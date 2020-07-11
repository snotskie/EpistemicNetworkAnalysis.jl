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
function rotate!(rotation::AbstractMeansRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)

    ## Manually factor the grouping variable to 0/1/missing
    metadata[!, :FactoredGroupVar] = map(eachrow(metadata)) do unitRow
        if unitRow[rotation.groupVar] == rotation.controlGroup
            return 0.0
        elseif unitRow[rotation.groupVar] == rotation.treatmentGroup
            return 1.0
        else
            return missing
        end
    end

    ## Use a FormulaRotation to do the rest of the work
    invoke(rotate!, Tuple{AbstractFormulaRotation, DataFrame, DataFrame, DataFrame}, rotation, networkModel, unitModel, metadata)
end

# Override plotting pieces
# ## Units - different default for labels
function plot_units!(p::Plot, ena::AbstractENAModel{<:AbstractMeansRotation}, displayRows::Array{Bool,1};
    flipX::Bool=false, flipY::Bool=false, minLabel::Union{Nothing,String}=nothing, maxLabel::Union{Nothing,String}=nothing,
    kwargs...)

    if isnothing(minLabel)
        minLabel = "$(ena.rotation.controlGroup) Units"
    end

    if isnothing(maxLabel)
        maxLabel = "$(ena.rotation.treatmentGroup) Units"
    end

    invoke(plot_units!, Tuple{Plot, AbstractENAModel{<:AbstractFormulaRotation}, Array{Bool,1}},
        p, ena, displayRows; flipX=flipX, flipY=flipY, minLabel=minLabel, maxLabel=maxLabel, kwargs...)
end

## CIs - we can color them into two groups
function plot_intervals!(p::Plot, ena::AbstractENAModel{<:AbstractMeansRotation}, displayRows::Array{Bool,1};
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    displayCentroids = ena.centroidModel[displayRows, :]
    displayMetadata = ena.metadata[displayRows, :]
    controlRows = map(x->x[ena.rotation.groupVar] == ena.rotation.controlGroup, eachrow(displayMetadata))
    treatmentRows = map(x->x[ena.rotation.groupVar] == ena.rotation.treatmentGroup, eachrow(displayMetadata))

    controlUnits = displayCentroids[controlRows, :]
    treatmentUnits = displayCentroids[treatmentRows, :]

    xs = controlUnits[!, :pos_x] * (flipX ? -1 : 1)
    ys = controlUnits[!, :pos_y] * (flipY ? -1 : 1)
    help_plot_ci(p, xs, ys, :purple, :square, "$(ena.rotation.controlGroup) Mean")

    xs = treatmentUnits[!, :pos_x] * (flipX ? -1 : 1)
    ys = treatmentUnits[!, :pos_y] * (flipY ? -1 : 1)
    help_plot_ci(p, xs, ys, :orange, :square, "$(ena.rotation.treatmentGroup) Mean")
end

## Network - we can do a subtraction plot
function plot_network!(p::Plot, ena::AbstractENAModel{<:AbstractMeansRotation}, displayRows::Array{Bool,1};
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    displayAccum = ena.accumModel[displayRows, :]
    displayMetadata = ena.metadata[displayRows, :]
    controlRows = map(x->x[ena.rotation.groupVar] == ena.rotation.controlGroup, eachrow(displayMetadata))
    treatmentRows = map(x->x[ena.rotation.groupVar] == ena.rotation.treatmentGroup, eachrow(displayMetadata))

    lineWidths = map(eachrow(ena.networkModel)) do networkRow
        controlMean = mean(displayAccum[controlRows, networkRow[:relationship]])
        treatmentMean = mean(displayAccum[treatmentRows, networkRow[:relationship]])
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
    codeWidths = zeros(nrow(ena.codeModel))
    for (i, networkRow) in enumerate(eachrow(ena.networkModel))
        j, k = ena.relationshipMap[networkRow[:relationship]]
        codeWidths[j] += lineWidths[i]
        codeWidths[k] += lineWidths[i]

        x = ena.codeModel[[j, k], :pos_x] * (flipX ? -1 : 1)
        y = ena.codeModel[[j, k], :pos_y] * (flipY ? -1 : 1)
        plot!(p, x, y,
            label=nothing,
            seriestype=:line,
            linewidth=lineWidths[i],
            linecolor=lineColors[i])
    end

    codeWidths *= 8 / maximum(codeWidths)
    x = ena.codeModel[!, :pos_x] * (flipX ? -1 : 1)
    y = ena.codeModel[!, :pos_y] * (flipY ? -1 : 1)
    labels = map(label->text(label, :top, 8), ena.codeModel[!, :code])
    plot!(p, x, y,
        label=nothing,
        seriestype=:scatter,
        series_annotations=labels,
        markershape=:circle,
        markersize=codeWidths,
        markercolor=:black,
        markerstrokecolor=:black)
end