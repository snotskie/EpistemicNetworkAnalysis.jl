"""
TODO document
"""
struct SVDRotation <: AbstractSVDRotation
end

# Implement rotation
function rotate!(rotation::AbstractSVDRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)
    pcaModel = projection(help_deflating_svd(networkModel, unitModel))
    networkModel[!, :weight_x] = pcaModel[:, 1]
    networkModel[!, :weight_y] = pcaModel[:, 2]
end

# Override plotting pieces
## Units - we can color them into n groups
function plot_units!(p::Plot, ena::AbstractENAModel{<:AbstractSVDRotation}, displayRows::Array{Bool,1};
    flipX::Bool=false, flipY::Bool=false, groupVar::Union{Symbol,Nothing}=nothing,
    kwargs...)

    displayCentroids = ena.centroidModel[displayRows, :]
    displayMetadata = ena.metadata[displayRows, :]
    unitColors = :black
    x = displayCentroids[!, :pos_x] * (flipX ? -1 : 1)
    y = displayCentroids[!, :pos_y] * (flipY ? -1 : 1)
    lastLabel = "Units"
    if !isnothing(groupVar)
        colors = [:purple, :orange, :green, :blue, :pink, :cyan]
        groups = sort(unique(ena.metadata[!, groupVar]))
        if length(groups) <= length(colors)
            colorMap = Dict(g => colors[i] for (i, g) in enumerate(groups))
            unitColors = map(eachrow(displayMetadata)) do unitRow
                return colorMap[unitRow[groupVar]]
            end

            lastLabel = nothing
            for g in keys(colorMap)
                plot!(p, x, y,
                    label="$(g) Units",
                    seriestype=:scatter,
                    markershape=:circle,
                    markersize=2,
                    markercolor=colorMap[g],
                    markerstrokecolor=colorMap[g])
            end
        else
            @warn "Too many groups to color nicely, falling back to all black"
        end
    end

    plot!(p, x, y,
        label=lastLabel,
        seriestype=:scatter,
        markershape=:circle,
        markersize=2,
        markercolor=unitColors,
        markerstrokecolor=unitColors)
end

## CIs - we can color them into n groups
function plot_intervals!(p::Plot, ena::AbstractENAModel{<:AbstractSVDRotation}, displayRows::Array{Bool,1};
    flipX::Bool=false, flipY::Bool=false, groupVar::Union{Symbol,Nothing}=nothing,
    kwargs...)

    displayCentroids = ena.centroidModel[displayRows, :]
    displayMetadata = ena.metadata[displayRows, :]
    if !isnothing(groupVar)
        colors = [:purple, :orange, :green, :blue, :pink, :cyan]
        groups = sort(unique(ena.metadata[!, groupVar]))
        if length(groups) <= length(colors)
            for (i, g) in enumerate(groups)
                color = colors[i]
                groupRows = map(eachrow(displayMetadata)) do unitRow
                    return unitRow[groupVar] == g
                end

                groupCentroids = displayCentroids[groupRows, :]
                xs = groupCentroids[!, :pos_x] * (flipX ? -1 : 1)
                ys = groupCentroids[!, :pos_y] * (flipY ? -1 : 1)
                help_plot_ci(p, xs, ys, color, :square, "$(g) Mean")
            end
        end
    end
end