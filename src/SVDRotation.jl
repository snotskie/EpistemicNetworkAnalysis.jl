"""
TODO document
"""
struct SVDRotation <: AbstractSVDRotation
    # no fields
end

# Implement rotation
function rotate!(rotation::AbstractSVDRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)

    ## Run an ortho svd and use those values as the axis weights
    pcaModel = projection(help_deflating_svd(networkModel, unitModel))
    networkModel[!, :weight_x] = pcaModel[:, 1]
    networkModel[!, :weight_y] = pcaModel[:, 2]
end

# Override plotting pieces
## Units - we can color them into n groups
function plot_units!(p::Plot, ena::AbstractENAModel{<:AbstractSVDRotation}, displayRows::Array{Bool,1};
    flipX::Bool=false, flipY::Bool=false, groupVar::Union{Symbol,Nothing}=nothing,
    kwargs...)

    ### Grab filtered values
    displayCentroids = ena.centroidModel[displayRows, :]
    displayMetadata = ena.metadata[displayRows, :]

    ### Default, all units are black
    unitColors = :black

    ### Find positions of units to show
    x = displayCentroids[!, :pos_x] * (flipX ? -1 : 1)
    y = displayCentroids[!, :pos_y] * (flipY ? -1 : 1)

    ### Default, show "Units" in the legend
    lastLabel = "Units"

    ### If we have a groupVar from the user...
    if !isnothing(groupVar)

        ### ...and we don't have more groups than we have colors...
        colors = [:purple, :orange, :green, :blue, :pink, :cyan]
        groups = sort(unique(ena.metadata[!, groupVar]))
        if length(groups) <= length(colors)

            ### ...then map each unit to its corresponding color...
            colorMap = Dict(g => colors[i] for (i, g) in enumerate(groups))
            unitColors = map(eachrow(displayMetadata)) do unitRow
                return colorMap[unitRow[groupVar]]
            end

            ### ...and add that color to the legend
            lastLabel = nothing
            for g in keys(colorMap)
                plot!(p, [-999], [-999],
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

    ### Draw the units with the colors found above
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

    ### Grab filtered values
    displayCentroids = ena.centroidModel[displayRows, :]
    displayMetadata = ena.metadata[displayRows, :]

    ### Only if we have a groupVar from the user...
    if !isnothing(groupVar)

        ### ...and we don't have more groups than colors...
        colors = [:purple, :orange, :green, :blue, :pink, :cyan]
        groups = sort(unique(ena.metadata[!, groupVar]))
        if length(groups) <= length(colors)

            ### ...then for each group...
            for (i, g) in enumerate(groups)

                ### ...grab the units in that group...
                groupRows = map(eachrow(displayMetadata)) do unitRow
                    return unitRow[groupVar] == g
                end

                groupCentroids = displayCentroids[groupRows, :]

                ### ...and show their CI
                xs = groupCentroids[!, :pos_x] * (flipX ? -1 : 1)
                ys = groupCentroids[!, :pos_y] * (flipY ? -1 : 1)
                help_plot_ci(p, xs, ys, colors[i], :square, "$(g) Mean")
            end
        end
    end
end