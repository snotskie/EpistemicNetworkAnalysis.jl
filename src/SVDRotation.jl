"""
TODO document
"""
struct SVDRotation <: AbstractSVDRotation
end

# Implement rotation
function rotate!(rotation::SVDRotation, networkModel::DataFrame, centroidModel::DataFrame)
    pcaModel = projection(help_deflating_svd(networkModel, centroidModel))
    networkModel[!, :weight_x] = pcaModel[:, 1]
    networkModel[!, :weight_y] = pcaModel[:, 2]
end

# Override plotting pieces
## Units - we can color them into n groups
function plot_units!(p::Plots.Plot, ena::AbstractENAModel{SVDRotation}, displayCentroids::DataFrame, displayCounts::DataFrame;
    flipX::Bool=false, flipY::Bool=false, groupVar::Union{Symbol,Nothing}=nothing,
    kwargs...)

    unitColors = :black
    if !isnothing(groupVar)
        colors = [:purple, :orange, :green, :blue, :pink, :cyan]
        groups = sort(unique(centroids(ena)[!, groupVar]))
        if length(groups) <= length(colors)
            colorMap = Dict(g => colors[i] for (i, g) in enumerate(groups))
            unitColors = map(eachrow(centroids(ena))) do unitRow
                return colorMap[unitRow[groupVar]]
            end
        else
            @warn "Too many groups to color nicely, falling back to all black"
        end
    end

    x = displayCentroids[!, :pos_x] * (flipX ? -1 : 1)
    y = displayCentroids[!, :pos_y] * (flipY ? -1 : 1)
    Plots.plot!(p, x, y,
        seriestype=:scatter,
        markershape=:circle,
        markersize=1.5,
        markercolor=unitColors,
        markerstrokecolor=unitColors)
end

## CIs - we can color them into n groups
function plot_intervals!(p::Plots.Plot, ena::AbstractENAModel{SVDRotation}, displayCentroids::DataFrame, displayCounts::DataFrame;
    flipX::Bool=false, flipY::Bool=false, groupVar::Union{Symbol,Nothing}=nothing,
    kwargs...)

    if !isnothing(groupVar)
        colors = [:purple, :orange, :green, :blue, :pink, :cyan]
        groups = sort(unique(centroids(ena)[!, groupVar]))
        if length(groups) <= length(colors)
            for (i, g) in enumerate(groups)
                color = colors[i]
                groupUnits = filter(displayCentroids) do unitRow
                    return unitRow[groupVar] == g
                end

                xs = groupUnits[!, :pos_x] * (flipX ? -1 : 1)
                ys = groupUnits[!, :pos_y] * (flipY ? -1 : 1)
                help_plot_ci(p, xs, ys, color, :square)
            end
        end
    end
end