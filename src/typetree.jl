# Rotations
abstract type AbstractENARotation end
abstract type AbstractSVDRotation <: AbstractENARotation end
abstract type AbstractFormulaRotation <: AbstractENARotation end
abstract type AbstractFormula2Rotation <: AbstractFormulaRotation end
abstract type AbstractMeansRotation <: AbstractFormulaRotation end
abstract type AbstractDoubleMeansRotation <: AbstractFormula2Rotation end

# Models
abstract type AbstractENAModel{T<:AbstractENARotation} end

# Default Functions
## Rotations
function rotate!(rotation::AbstractENARotation, networkModel::DataFrame, centroidModel::DataFrame)
    error("Unimplemented")
end

## Model Access
function centroids(ena::AbstractENAModel)
    return ena.centroidModel #! assumes this field by default
end

function counts(ena::AbstractENAModel)
    return ena.countModel #! assumes this field by default
end

function network(ena::AbstractENAModel)
    return ena.networkModel #! assumes this field by default
end

function relationships(ena::AbstractENAModel)
    return ena.relationshipMap #! assumes this field by default
end

function nodes(ena::AbstractENAModel)
    return ena.codeModel #! assumes this field by default
end

function pearson(ena::AbstractENAModel)
    return ena.pearson #! assumes this field by default
end

## Text display
function Base.display(ena::AbstractENAModel) # TODO should this be print, display, or show?
    println("Units (centroids):")
    show(centroids(ena)[!, [:ENA_UNIT, :pos_x, :pos_y]], allrows=true)
    println()
    println("Codes:")
    show(nodes(ena), allrows=true)
    println()
    println("Network:")
    show(network(ena), allrows=true)
    println()
    println("Model fit (Pearson):")
    println(round(pearson(ena), digits=4))
    println()
end

## Plotting
### Top-level wrapper
function plot(ena::AbstractENAModel;
    margin=10mm, size=500, lims=1, ticks=[-1, 0, 1], title="", 
    kwargs...)

    p = plot(leg=false, margin=margin, size=(size, size))
    plot!(p, ena; kwargs...)
    xticks!(p, ticks)
    yticks!(p, ticks)
    xlims!(p, -lims, lims)
    ylims!(p, -lims, lims)
    title!(p, title)
    return p
end

### Mutating wrapper
function plot!(p::Plot, ena::AbstractENAModel;
    showUnits::Bool=true, showNetwork::Bool=true, showIntervals::Bool=true, showExtras::Bool=true,
    display_filter=x->true,
    kwargs...)

    goodRows = map(display_filter, eachrow(centroids(ena)))
    displayCentroids = centroids(ena)[goodRows, :]
    displayCounts = counts(ena)[goodRows, :]
    if showUnits
        plot_units!(p, ena, displayCentroids, displayCounts; kwargs...)
    end

    if showNetwork
        plot_network!(p, ena, displayCentroids, displayCounts; kwargs...)
    end

    if showIntervals
        plot_intervals!(p, ena, displayCentroids, displayCounts; kwargs...)
    end

    if showExtras
        plot_extras!(p, ena, displayCentroids, displayCounts; kwargs...)
    end

    plot_labels!(p, ena; kwargs...)
end

### Unit-level helper
function plot_units!(p::Plot, ena::AbstractENAModel, displayCentroids::DataFrame, displayCounts::DataFrame;
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    x = displayCentroids[!, :pos_x] * (flipX ? -1 : 1)
    y = displayCentroids[!, :pos_y] * (flipY ? -1 : 1)
    plot!(p, x, y,
        seriestype=:scatter,
        markershape=:circle,
        markersize=1.5,
        markercolor=:black,
        markerstrokecolor=:black)
end

### Network-level helper
function plot_network!(p::Plot, ena::AbstractENAModel, displayCentroids::DataFrame, displayCounts::DataFrame;
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    lineWidths = map(eachrow(network(ena))) do networkRow
        return sum(displayCounts[!, networkRow[:relationship]])
    end

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
            linecolor=:black)
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

### CI-level helper
function plot_intervals!(p::Plot, ena::AbstractENAModel, displayCentroids::DataFrame, displayCounts::DataFrame;
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)
    # do nothing
end

### Extras helper
function plot_extras!(p::Plot, ena::AbstractENAModel, displayCentroids::DataFrame, displayCounts::DataFrame;
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)
    # do nothing
end

### Labels helper
function plot_labels!(p::Plot, ena::AbstractENAModel;
    xlabel="X", ylabel="Y",
    kwargs...)

    total_variance = sum(var.(eachcol(centroids(ena)[!, network(ena)[!, :relationship]])))
    variance_x = var(centroids(ena)[!, :pos_x]) / total_variance
    variance_y = var(centroids(ena)[!, :pos_y]) / total_variance
    xlabel!(p, "$xlabel ($(round(Int, variance_x*100))%)")
    ylabel!(p, "$ylabel ($(round(Int, variance_y*100))%)")
end