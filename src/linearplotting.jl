# Defaults
const DEFAULT_NEG_COLOR = colorant"#cc423a"
const DEFAULT_POS_COLOR = colorant"#218ebf"
const DEFAULT_GROUP_COLORS = [
    colorant"#56BD7C",
    colorant"#EF691B",
    colorant"#9d5dbb",
    colorant"#FBC848",
    colorant"#D0386C",
    colorant"#f18e9f",
    colorant"#9A9EAB",
    colorant"#ff8c39",
    colorant"#346B88"
]

const DEFAULT_ALPHABET = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZθωερτψυπασδφγηςκλζχξωβνμΘΩΨΥΠΣΔΦΓΛΞΩ"
const GLOBAL_MAX_EDGE_SIZE = 8
const GLOBAL_MAX_NODE_SIZE = 8
const GLOBAL_UNIT_SIZE = 3
const DEFAULT_INNER_MEAN_MARKERS = [
    :utriangle, :dtriangle, :rtriangle, :ltriangle,
    :pentagon, :heptagon, :octagon, :star4,
    :star6, :star7, :star8
]

function defaultplotkwargs(
        ::Type{M},
        model::AbstractLinearENAModel;
        x::Int=1,
        y::Int=2,
        margin::PlotMeasures.AbsoluteLength=10mm,
        size::Real=700,
        meanCenter::Bool=model.config.sphereNormalize,
        origin::Array{<:Real}=(meanCenter ?  [mean(model.points[x, :]), mean(model.points[y, :])] : [0,0]),
        zoom::Real=1,
        lims::Real=1/zoom,
        flipX::Bool=false,
        flipY::Bool=false,
        xticks::Array{<:Real}=(
            round.([origin[1]-lims, origin[1], origin[1]+lims], digits=4) |>
            (xticks) -> (flipX ? -reverse(xticks) : xticks)
        ),
        yticks::Array{<:Real}=(
            round.([origin[2]-lims, origin[2], origin[2]+lims], digits=4) |>
            (yticks) -> (flipY ? -reverse(yticks) : yticks)
        ),
        xlims::Array{<:Real}=xticks[[1, end]],
        ylims::Array{<:Real}=yticks[[1, end]],
        titles::Array{<:AbstractString}=String[],
        xlabel::AbstractString=model.embedding[x, :label],
        ylabel::AbstractString=model.embedding[y, :label],
        unitLabel::AbstractString="Unit",
        leg::Union{Symbol,Bool}=:topleft,
        negColor::Colorant=DEFAULT_NEG_COLOR,
        posColor::Colorant=DEFAULT_POS_COLOR,
        groupColors::Array{<:Colorant,1}=DEFAULT_GROUP_COLORS,
        alphabet::String=DEFAULT_ALPHABET,
        groupBy::Union{Symbol,Nothing}=nothing,
        innerGroupBy::Union{Symbol,Nothing}=nothing,
        spectralColorBy::Union{Symbol,Nothing}=nothing,
        trajectoryBy::Union{Symbol,Nothing}=nothing,
        trajectoryBins::Int=5,
        spectoryBy::Union{Symbol,Nothing}=nothing,
        spectoryPercent::Real=1/3,
        showExtras::Bool=true,
        showNetworks::Bool=true,
        showUnits::Bool=true,
        showMeans::Bool=true,
        showWarps::Bool=false,
        confidenceShape::Symbol=:rect,
        fitNodesToCircle::Bool=false,
        showWeakEdges::Bool=true,
        colorbar::Bool=false,
        kwargs...
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    kwargs = NamedTuple(kwargs)
    defaults = (
        x=x,
        y=y,
        margin=margin,
        size=size,
        meanCenter=meanCenter,
        origin=origin,
        zoom=zoom,
        lims=lims,
        xticks=xticks,
        yticks=yticks,
        xlims=xlims,
        ylims=ylims,
        titles=titles,
        xlabel=xlabel,
        ylabel=ylabel,
        unitLabel=unitLabel,
        leg=leg,
        negColor=negColor,
        posColor=posColor,
        groupColors=groupColors,
        alphabet=alphabet,
        flipX=flipX,
        flipY=flipY,
        groupBy=groupBy,
        innerGroupBy=innerGroupBy,
        spectralColorBy=spectralColorBy,
        trajectoryBy=trajectoryBy,
        trajectoryBins=trajectoryBins,
        spectoryBy=spectoryBy,
        spectoryPercent=spectoryPercent,
        showExtras=showExtras,
        showNetworks=showNetworks,
        showUnits=showUnits,
        showMeans=showMeans,
        showWarps=showWarps,
        confidenceShape=confidenceShape,
        fitNodesToCircle=fitNodesToCircle,
        showWeakEdges=showWeakEdges,
        colorbar=colorbar,
        kwargs...
    )

    return merge(defaults, kwargs)
end

# Top-level
function plot(
        ::Type{M}, model::AbstractLinearENAModel, plotconfig::NamedTuple
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    # Produce the plots by calling helpers
    ps = Plot[]
    plot_omnibus!(M, ps, model, plotconfig)
    plot_trends!(M, ps, model, plotconfig)
    plot_groups!(M, ps, model, plotconfig)
    plot_subtractions!(M, ps, model, plotconfig)
    plot_extras!(M, ps, model, plotconfig)
    return layoutSubplots(ps, plotconfig)
end

# Documentation
"""
    plot(
        model::AbstractLinearENAModel,
        x::Int=1,
        y::Int=2,
        margin::PlotMeasures.AbsoluteLength=10mm,
        size::Real=600,
        meanCenter::Bool=model.config.sphereNormalize,
        origin::Array{<:Real}=(meanCenter ?  [mean(model.points[x, :]), mean(model.points[y, :])] : [0,0]),
        zoom::Real=1,
        lims::Real=1/zoom,
        flipX::Bool=false,
        flipY::Bool=false,
        xticks::Array{<:Real}=(
            round.([origin[1]-lims, origin[1], origin[1]+lims], digits=4) |>
            (xticks) -> (flipX ? -reverse(xticks) : xticks)
        ),
        yticks::Array{<:Real}=(
            round.([origin[2]-lims, origin[2], origin[2]+lims], digits=4) |>
            (yticks) -> (flipY ? -reverse(yticks) : yticks)
        ),
        xlims::Array{<:Real}=xticks[[1, end]],
        ylims::Array{<:Real}=yticks[[1, end]],
        titles::Array{<:AbstractString}=String[],
        xlabel::AbstractString=model.embedding[x, :label],
        ylabel::AbstractString=model.embedding[y, :label],
        unitLabel::AbstractString="Unit",
        leg::Union{Symbol,Bool}=:topleft,
        negColor::Colorant=DEFAULT_NEG_COLOR,
        posColor::Colorant=DEFAULT_POS_COLOR,
        groupColors::Array{<:Colorant,1}=DEFAULT_GROUP_COLORS,
        alphabet::String=DEFAULT_ALPHABET,
        groupBy::Union{Symbol,Nothing}=nothing,
        innerGroupBy::Union{Symbol,Nothing}=nothing,
        spectralColorBy::Union{Symbol,Nothing}=nothing,
        trajectoryBy::Union{Symbol,Nothing}=nothing,
        trajectoryBins::Int=5,
        spectoryBy::Union{Symbol,Nothing}=nothing,
        spectoryPercent::Real=1/3,
        showExtras::Bool=true,
        showNetworks::Bool=true,
        showUnits::Bool=true,
        showMeans::Bool=true,
        showWarps::Bool=false,
        confidenceShape::Symbol=:rect,
        fitNodesToCircle::Bool=false,
        showWeakEdges::Bool=true,
        colorbar::Bool=false
    )

Plot an ENA model using the [GR backend](https://docs.juliaplots.org/latest/gallery/gr/)

See also [savefig](https://docs.juliaplots.org/latest/output/#savefig-/-format)

## Arguments

At minimum, the only required argument is the ENA model itself.

Several optional arguments are available:

- `x` and `y` control which dimension to show on the x- and y-axis respectively
- `margin`, `size`, `meanCenter`, `origin`, `zoom`, `lims`, `flipX`, `flipY`, `xticks`, `yticks`, `xlims`, and `ylims` together control aspects of the plot size and axes
- `titles`, `xlabel`, `ylabel`, `unitLabel`, `leg`, and `alphabet` together control the text that labels the plot
- `negColor`, `posColor`, and `groupColors` together control the colors used in the plot
- `groupBy` and `innerGroupBy` define which metadata columns to use as grouping variables for the sake of color coding and confidence intervals
- `spectralColorBy` defines which metadata column to use to color-code units as a spectrum, to show how networks relate to the variable of interest
- `trajectoryBy` and `trajectoryBins` together define and control how a trajectory path should be overlaid on the plot, to show how the mean network changes along the variable of interest. Similarly, `spectoryBy` and `spectoryPercent` define how a sequence of 1-level densities should be overlaid on the plot, to show how the distribution of networks changes along the variable of interest
- `showExtras`, `showNetworks`, `showUnits`, and `showMeans` control which plot elements to show or hide. Additionally, `confidenceShape` can be set to `:rect` (default) or `:density` to choose which shape to use around the means
- `showWarps` controls if edges should be drawn straight (`false`) or "warped" to show their true location in the space (`true`)
- `fitNodesToCircle` controls if nodes should be shown in their optimized positions for goodness of fit, or at a circular position around the origin
- `showWeakEdges` controls if edges with weak correlations to trends should be shown
- `colorbar` controls if subplots that use gradient color-coding for their edges should have an explicit colorbar added

## Example

```julia
model = ENAModel(data, codes, conversations, units)
p = plot(model)

# Move the legends to an outer position
p = plot(model, leg=:outertopright)

# Grab one subplot
sp = plot(p.subplots[1], size=(700,700))

# Save
savefig(p, "example.png")
```

## Interpretation

`plot(model)` produces a plot with the following subplots:

- `(a)` an overall mean, which tells us the baseline everything compares against
- `(b)` and `(c)` rates of change for each connection across the x- and y-axes, which helps identify what is being modeled by each axis
- Subsequent subplots show each subgroup on its own. It's good to compare these to the overall mean
- And the last subplots show how each pair of subgroups compare. Similar to the trend plots, these show you what is being modeled by the difference of the two groups

Some differences from WebENA and rENA:

- Saturation shows *correlation strength*. In WebENA and rENA saturation and line thickness both show magnitude of an effect
- Plots are mean centered by moving the origin of the plot, not by changing the underlying data. This preserves information that may or may not be useful for downstream analyses
- Plots are opinionated. Based on the model config, the plot's default settings to change to what I believed was the best way to plot that kind of model. This gives you the "right" plot without having to specify what "right" means each time
- A [known issue](https://github.com/snotskie/EpistemicNetworkAnalysis.jl/issues/11) is that the y-axis label can get cutoff when there are a lot of subplots
"""
plot(model::AbstractLinearENAModel)

# Helpers
function fixX(x::Real, plotconfig::NamedTuple)
    return x * (plotconfig.flipX ? -1 : 1)
end

function fixY(y::Real, plotconfig::NamedTuple)
    return y * (plotconfig.flipY ? -1 : 1)
end

function fixPoint(point::Array{<:Real}, plotconfig::NamedTuple)
    return [
        fixX(point[1], plotconfig),
        fixY(point[2], plotconfig)
    ]
end

function nonlinearGradientMap(lo, mid, hi; grains=100, curve=1.5)
    return vcat(
        [weighted_color_mean((grains-i)^curve/grains^curve, lo, mid) for i in 1:grains],
        [mid],
        [weighted_color_mean(1-i^curve/grains^curve, mid, hi) for i in 1:grains]
    )
end

function layoutSubplots(ps::Array{Plot}, plotconfig::NamedTuple)
    for p in ps
        xticks!(p, plotconfig.xticks)
        yticks!(p, plotconfig.yticks)
        if plotconfig.lims > 0
            xlims!(p, plotconfig.xlims...)
            ylims!(p, plotconfig.ylims...)
        end
    end

    N = ceil(Int, sqrt(length(ps)))
    M = ceil(Int, length(ps)/N)
    layout = grid(N, M)
    while length(ps) < N*M
        push!(ps, plot(legend=false,grid=false,foreground_color_subplot=:white))
    end

    # BUGFIX: divide the default arrow size by ceil(sqrt(number of plots in the multi plot)) = number of plots on the first row
    GR.setarrowsize(1/N)

    return plot(ps..., size=(plotconfig.size*M, plotconfig.size*N), layout=layout)
end

function computeNetworkDensities(model, rows=!; normalize=false)

    # sum densities for each edge
    edgeIDs = model.edges.edgeID
    edgeDensities = Dict(
        edgeID => sum(model.accum[rows, edgeID])
        for edgeID in edgeIDs
    )

    # find the density of each code row dot, by "splitting" the density of each line between its two codes
    nodeIDs = model.nodes.nodeID
    nodeDensities = Dict(
        nodeID => 0.0
        for nodeID in nodeIDs
    )

    for edge in eachrow(model.edges)
        nodeDensities[edge.ground] += edgeDensities[edge.edgeID]
        nodeDensities[edge.response] += edgeDensities[edge.edgeID]
    end

    # optionally normalize each
    if normalize
        s = maximum(values(edgeDensities))
        if s != 0
            for edgeID in edgeIDs
                edgeDensities[edgeID] /= s
            end
        end

        s = maximum(values(nodeDensities))
        if s != 0
            for nodeID in nodeIDs
                nodeDensities[nodeID] /= s
            end
        end
    end

    return edgeDensities, nodeDensities
end

function getGroupColorMap(model::AbstractLinearENAModel, plotconfig::NamedTuple)
    allGroups = sort(unique(model.metadata[!, plotconfig.groupBy]))
    groupColors = Dict(
        g => plotconfig.groupColors[i]
        for (i, g) in enumerate(allGroups)
        if i <= length(plotconfig.groupColors) # in case we don't have enough colors
    )

    return DefaultDict(colorant"black", groupColors) # in case we don't have enough colors  
end

function getSpectralColors(vals)
    # return [HSL(val*240, 1, 0.5) for val in vals]
    # return [HSL(val*240 + 120, 1, 0.5) for val in vals]
    # return [HSL(30, 1, val*.6 + .2) for val in vals]
    return [HSL(200 - val*170, .9, .25 + val*.2) for val in vals]
end

function plotConfidenceInterval(p, xs, ys, color, shape, label, ci_shape)

    # Plot the mean center
    if length(xs) > 0
        x = mean(xs)
        y = mean(ys)
        Plots.plot!(
            p, [x], [y],
            label=label,
            seriestype=:scatter,
            markersize=4,
            markershape=shape,
            markercolor=color,
            markerstrokecolor=color
        )
    end

    # Plot the confidence interval
    if length(xs) > 1
        if ci_shape == :rect
            ci_x = collect(confint(OneSampleTTest(xs)))
            ci_y = collect(confint(OneSampleTTest(ys)))
            Plots.plot!(p, [ci_x[1], ci_x[2]], [ci_y[1], ci_y[1]],
                label=nothing,
                seriestype=:line,
                linewidth=1,
                linecolor=color)

            Plots.plot!(p, [ci_x[1], ci_x[2]], [ci_y[2], ci_y[2]],
                label=nothing,
                seriestype=:line,
                linewidth=1,
                linecolor=color)

            Plots.plot!(p, [ci_x[1], ci_x[1]], [ci_y[1], ci_y[2]],
                label=nothing,
                seriestype=:line,
                linewidth=1,
                linecolor=color)

            Plots.plot!(p, [ci_x[2], ci_x[2]], [ci_y[1], ci_y[2]],
                label=nothing,
                seriestype=:line,
                linewidth=1,
                linecolor=color)
        elseif ci_shape == :density
            plot_kde_without_cbar!(p,
                kde((xs, ys)),
                levels=[1, 2, 3, 4, 5],
                color=color
            )
        # elseif ci_shape == :ellipse
        #     # remove outliers, as the ellipse is sensitive to them
        #     q1x = quantile(xs, .25)
        #     q3x = quantile(xs, .75)
        #     iqrx = q3x - q1x
        #     inx = map(xs) do x
        #         if x < q1x - 1.5 * iqrx
        #             return false
        #         elseif x > q3x + 1.5 * iqrx
        #             return false
        #         else
        #             return true
        #         end
        #     end
        #     q1y = quantile(ys, .25)
        #     q3y = quantile(ys, .75)
        #     iqry = q3y - q1y
        #     iny = map(ys) do y
        #         if y < q1y - 1.5 * iqry
        #             return false
        #         elseif y > q3y + 1.5 * iqry
        #             return false
        #         else
        #             return true
        #         end
        #     end
        #     inxy = inx .& iny
        #     xs, ys = xs[inxy], ys[inxy]

        #     # covariance matrix
        #     Σ = cov([xs ys])
        #     vals = eigvals(Σ)
        #     vecs = eigvecs(Σ)
        #     if vals[1] <= vals[2]
        #         vals[1], vals[2] = vals[2], vals[1]
        #         vecs[:, 1], vecs[:, 2] = vecs[:, 2], vecs[:, 1]
        #     end

        #     # angle
        #     θ = atan(vecs[2,1] / vecs[1,1])
        #     if θ < 0
        #         θ += 2π
        #     end

        #     # center
        #     cx = mean(xs)
        #     cy = mean(ys)

        #     # ellipse points
        #     t = range(0, 2π, length=100)
        #     confidence = 0.95
        #     quant = sqrt(quantile(Chisq(2), confidence))
        #     rx = quant * sqrt(vals[1])
        #     ry = quant * sqrt(vals[2])
        #     xe = rx .* cos.(t)
        #     ye = ry .* sin.(t)
        #     ellipse = [xe ye] * [cos(θ) sin(θ); -sin(θ) cos(θ)]
        #     pxs = cx .+ ellipse[:, 1]
        #     pys = cy .+ ellipse[:, 2]

        #     # plot it
        #     Plots.plot!(p, pxs, pys,
        #         label=nothing,
        #         color=color)
        #         # seriestype=:line,
        #         # linewidth=1,
        #         # linecolor=color)
        else
            error("Unrecognized confidenceShape $(ci_shape). Available options are :rect and :density")
        end
    end
end

function getAxisLabels(model::AbstractLinearENAModel, plotconfig::NamedTuple)
    return (
        xlabel="$(plotconfig.xlabel) ($(round(100*model.embedding[plotconfig.x, :].variance_explained, digits=1))%)",
        ylabel="$(plotconfig.ylabel) ($(round(100*model.embedding[plotconfig.y, :].variance_explained, digits=1))%)"
    )
end

function paintSortedNetwork!(
        p::Plot,
        model::AbstractLinearENAModel,
        plotconfig::NamedTuple,
        edgeWidths::Dict,
        edgeColors::Dict,
        nodeWidths::Dict
    )

    edgeOrder = first.(sort(collect(edgeWidths), by=x->x[2]))
    edgeMap = Dict(
        edge.edgeID => i
        for (i, edge) in enumerate(eachrow(model.edges))
    )

    edgePoints = Dict(
        edgeID => []
        for edgeID in edgeOrder
    )

    for edge in eachrow(model.edges)
        push!(edgePoints[edge.edgeID], fixPoint([
            model.pointsNodes[plotconfig.x, edge.ground],
            model.pointsNodes[plotconfig.y, edge.ground]
        ], plotconfig))

        if edge.kind == :count
            push!(edgePoints[edge.edgeID], fixPoint([0, 0], plotconfig))
        elseif edge.kind in [:undirected, :directed]
            if plotconfig.showWarps
                pointT = fixPoint([
                    model.embedding[plotconfig.x, edge.edgeID],
                    model.embedding[plotconfig.y, edge.edgeID]
                ], plotconfig)
                
                push!(edgePoints[edge.edgeID], pointT)
                push!(edgePoints[edge.edgeID], pointT)
                push!(edgePoints[edge.edgeID], pointT)
            end

            push!(edgePoints[edge.edgeID], fixPoint([
                model.pointsNodes[plotconfig.x, edge.response],
                model.pointsNodes[plotconfig.y, edge.response]
            ], plotconfig))
        end

        if plotconfig.fitNodesToCircle
            for point in edgePoints[edge.edgeID]
                s = sqrt(sum(point .^ 2))
                if s != 0
                    point .= point / s
                end
            end
        end
    end

    for edgeID in edgeOrder
        edge = model.edges[edgeMap[edgeID], :]
        if edgeWidths[edgeID] > 0 && edge.kind != :directed
            plot!(p,
                first.(edgePoints[edgeID]), # xs
                last.(edgePoints[edgeID]), # ys
                label=nothing,
                seriestype=:curves,
                arrows=false,
                linewidth=edgeWidths[edgeID],
                linecolor=edgeColors[edgeID]
            )
        end
    end

    for nodeID in model.nodes.nodeID
        x = fixX(model.pointsNodes[plotconfig.x, nodeID], plotconfig)
        y = fixY(model.pointsNodes[plotconfig.y, nodeID], plotconfig)
        label = text(string(nodeID), :top, default(:xtickfontsize))
        if plotconfig.fitNodesToCircle
            s = sqrt(x^2 + y^2)
            if s != 0
                x /= s
                y /= s
            end

            angle = atan(y, x) * 180 / pi
            if x < 0
                angle = atan(-y, -x) * 180 / pi
            end
            
            label = text(string(nodeID), :top, default(:xtickfontsize), rotation=angle)
        end

        if nodeWidths[nodeID] > 0
            plot!(
                p, [x], [y],
                label=nothing,
                seriestype=:scatter,
                series_annotations=[label],
                markershape=:circle,
                markersize=[nodeWidths[nodeID]],
                markercolor=:black,
                markerstrokewidth=0
            )
        end
    end

    for edgeID in reverse(edgeOrder)
        edge = model.edges[edgeMap[edgeID], :]
        if edgeWidths[edgeID] > 0 && edge.kind == :directed
            plot!(p,
                first.(edgePoints[edgeID]), # xs
                last.(edgePoints[edgeID]), # ys
                label=nothing,
                seriestype=:curves,
                arrows=true,
                linewidth=edgeWidths[edgeID],
                linecolor=edgeColors[edgeID]
            )
        end
    end
end

# Painters
abstract type AbstractLinearEdgePainter <: AbstractEdgePainter end
struct FlatLinearEdgePainter <: AbstractLinearEdgePainter
    color::Colorant
end

function paint_edges!(
        p::Plot,
        edgePainter::FlatLinearEdgePainter,
        model::AbstractLinearENAModel,
        plotconfig::NamedTuple,
        displayRows::BitVector=BitVector(repeat([true], nrow(model.accum)))
    )

    # Find widths
    allEdgeWidths, allNodeWidths = computeNetworkDensities(model)
    edgeWidths, nodeWidths = computeNetworkDensities(model, displayRows)

    # Rescale
    s = maximum(values(allEdgeWidths))
    if s != 0
        for edgeID in model.edges.edgeID
            edgeWidths[edgeID] *= GLOBAL_MAX_EDGE_SIZE / s
        end
    end

    s = maximum(values(allNodeWidths))
    if s != 0
        for nodeID in model.nodes.nodeID
            nodeWidths[nodeID] *= GLOBAL_MAX_NODE_SIZE / s
        end
    end

    edgeColors = Dict(
        edgeID => edgePainter.color
        for edgeID in model.edges.edgeID
    )

    paintSortedNetwork!(p, model, plotconfig, edgeWidths, edgeColors, nodeWidths)
end

struct TrendLinearEdgePainter <: AbstractLinearEdgePainter
    vals::Vector{<:Union{Missing,Real}}
    neg_color::Colorant
    pos_color::Colorant
end

function paint_edges!(
        p::Plot,
        edgePainter::TrendLinearEdgePainter,
        model::AbstractLinearENAModel,
        plotconfig::NamedTuple,
        displayRows::BitVector=BitVector(repeat([true], nrow(model.accum)))
    )

    # Compile data for regressions
    regressionData = DataFrame(Dict(
        :vals => edgePainter.vals[displayRows],
        [
            # why map? Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
            edgeID => map(Float64, Vector(model.accum[displayRows, edgeID]))
            for edgeID in model.edges.edgeID
        ]...
    ))

    # Color map the lines based on their correlation with the vals
    mid_color = weighted_color_mean(0.5, edgePainter.neg_color, edgePainter.pos_color)
    mid_color = weighted_color_mean(0.3, mid_color, colorant"white")
    edgeColorMap = nonlinearGradientMap(
        # weighted_color_mean(0.95, edgePainter.neg_color, colorant"black"),
        edgePainter.neg_color,
        mid_color,
        # weighted_color_mean(0.95, edgePainter.pos_color, colorant"black")
        edgePainter.pos_color
    )

    # Compute line widths as the strength (slope) between the vals and the accum network weights
    edgeWidths = Dict(edgeID => 0.0 for edgeID in model.edges.edgeID)
    edgeColors = Dict(edgeID => colorant"black" for edgeID in model.edges.edgeID)
    for edgeID in model.edges.edgeID
        f1 = FormulaTerm(term(edgeID), term(:vals))
        try
            m1 = fit(LinearModel, f1, regressionData)
            slope = coef(m1)[2]
            pearson = cor(regressionData[!, edgeID], regressionData[!, :vals])
            if isnan(pearson)
                pearson = 0
            end

            if !plotconfig.showWeakEdges && abs(pearson) < 0.3
                slope = 0
            end

            edgeWidths[edgeID] = abs(slope)
            index = 1 + round(Int, (length(edgeColorMap) - 1) * (pearson + 1) / 2)
            edgeColors[edgeID] = edgeColorMap[index]
        catch e
            # do nothing
        end
    end

    # Compute node widths
    nodeWidths = Dict(nodeID => 0.0 for nodeID in model.nodes.nodeID)
    for edge in eachrow(model.edges)
        nodeWidths[edge.ground] += edgeWidths[edge.edgeID]
        nodeWidths[edge.response] += edgeWidths[edge.edgeID]
    end

    # Rescale
    s = maximum(values(edgeWidths))
    if s != 0
        for edgeID in model.edges.edgeID
            edgeWidths[edgeID] *= GLOBAL_MAX_EDGE_SIZE / s
        end
    end

    s = maximum(values(nodeWidths))
    if s != 0
        for nodeID in model.nodes.nodeID
            nodeWidths[nodeID] *= GLOBAL_MAX_NODE_SIZE / s
        end
    end

    paintSortedNetwork!(p, model, plotconfig, edgeWidths, edgeColors, nodeWidths)

    # Optionally add color bar
    if plotconfig.colorbar
        plot!(p,
            [], [], label=false,
            marker_z=-1:1,
            color=cgrad(edgeColorMap),
            colorbar=true,
            colorbar_title="Pearson r"
        )
    end
end

struct NodesOnlyLinearEdgePainter <: AbstractLinearEdgePainter end
function paint_edges!(
        p::Plot,
        edgePainter::NodesOnlyLinearEdgePainter,
        model::AbstractLinearENAModel,
        plotconfig::NamedTuple,
        displayRows::BitVector=BitVector(repeat([true], nrow(model.accum)))
    )

    # Find widths
    allEdgeWidths, allNodeWidths = computeNetworkDensities(model)
    edgeWidths, nodeWidths = computeNetworkDensities(model, displayRows)

    # Rescale
    s = maximum(values(allEdgeWidths))
    if s != 0
        for edgeID in model.edges.edgeID
            edgeWidths[edgeID] *= GLOBAL_MAX_EDGE_SIZE / s
        end
    end

    s = maximum(values(allNodeWidths))
    if s != 0
        for nodeID in model.nodes.nodeID
            nodeWidths[nodeID] *= GLOBAL_MAX_NODE_SIZE / s
        end
    end

    # "zero-out" the edges
    edgeWidths = Dict(
        edgeID => 0.0
        for edgeID in model.edges.edgeID
    )

    edgeColors = Dict(
        edgeID => colorant"transparent"
        for edgeID in model.edges.edgeID
    )

    paintSortedNetwork!(p, model, plotconfig, edgeWidths, edgeColors, nodeWidths)
end

# Subplots
function plot_omnibus!(
        ::Type{M}, ps::Array{Plot}, model::AbstractLinearENAModel, plotconfig::NamedTuple
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    p = plot(;
        leg=plotconfig.leg,
        aspect_ratio=:equal,
        margin=plotconfig.margin,
        yrotation=90,
        size=(plotconfig.size, plotconfig.size),
        getAxisLabels(model, plotconfig)...
    )

    letter = DEFAULT_ALPHABET[length(ps)+1]
    title!(p, "($(letter)) Grand Total")
    plot_network!(M, p, model, FlatLinearEdgePainter(colorant"#aaa"), plotconfig)
    plot_units!(M, p, model, plotconfig)
    plot_means!(M, p, model, plotconfig)
    push!(ps, p)
end

function plot_trends!(
        ::Type{M}, ps::Array{Plot}, model::AbstractLinearENAModel, plotconfig::NamedTuple
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    px = plot(
        leg=plotconfig.leg,
        aspect_ratio=:equal,
        margin=plotconfig.margin,
        yrotation=90,
        size=(plotconfig.size, plotconfig.size)
    )

    py = plot(
        leg=plotconfig.leg,
        aspect_ratio=:equal,
        margin=plotconfig.margin,
        yrotation=90,
        size=(plotconfig.size, plotconfig.size)
    )

    letterx = DEFAULT_ALPHABET[length(ps)+1]
    title!(px, "($(letterx)) Rate of Change by X")

    lettery = DEFAULT_ALPHABET[length(ps)+2]
    title!(py, "($(lettery)) Rate of Change by Y")

    unitIDs = model.accum.unitID
    xs = Vector(model.points[plotconfig.x, unitIDs])
    ys = Vector(model.points[plotconfig.y, unitIDs])
    plot_network!(M, px, model, TrendLinearEdgePainter(xs, plotconfig.negColor, plotconfig.posColor), plotconfig)
    plot_network!(M, py, model, TrendLinearEdgePainter(ys, plotconfig.negColor, plotconfig.posColor), plotconfig)
    plot!(px,
        [-999], [-999],
        label="Positive Change",
        seriestype=:scatter,
        markershape=:hline,
        markersize=GLOBAL_UNIT_SIZE,
        markercolor=[plotconfig.posColor],
        markerstrokewidth=1
    )
    plot!(px,
        [-999], [-999],
        label="Negative Change",
        seriestype=:scatter,
        markershape=:hline,
        markersize=GLOBAL_UNIT_SIZE,
        markercolor=[plotconfig.negColor],
        markerstrokewidth=1
    )
    plot!(py,
        [-999], [-999],
        label="Positive Change",
        seriestype=:scatter,
        markershape=:hline,
        markersize=GLOBAL_UNIT_SIZE,
        markercolor=[plotconfig.posColor],
        markerstrokewidth=1
    )
    plot!(py,
        [-999], [-999],
        label="Negative Change",
        seriestype=:scatter,
        markershape=:hline,
        markersize=GLOBAL_UNIT_SIZE,
        markercolor=[plotconfig.negColor],
        markerstrokewidth=1
    )
    
    push!(ps, px)
    push!(ps, py)
end

function plot_groups!(
        ::Type{M}, ps::Array{Plot}, model::AbstractLinearENAModel, plotconfig::NamedTuple
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    if isnothing(plotconfig.groupBy)
        return
    end

    groupColors = getGroupColorMap(model, plotconfig)
    for group in keys(groupColors)
        p = plot(
            leg=plotconfig.leg,
            aspect_ratio=:equal,
            margin=plotconfig.margin,
            yrotation=90,
            size=(plotconfig.size, plotconfig.size)
        )

        letter = DEFAULT_ALPHABET[length(ps)+1]
        title!(p, "($(letter)) $(group)")
        groupRows = model.metadata[!, plotconfig.groupBy] .== group
        plot_network!(M, p, model, FlatLinearEdgePainter(groupColors[group]), plotconfig, groupRows)
        plot_units!(M, p, model, plotconfig, groupRows)
        plot_means!(M, p, model, plotconfig, groupRows)
        push!(ps, p)
    end
end

function plot_subtractions!(
        ::Type{M}, ps::Array{Plot}, model::AbstractLinearENAModel, plotconfig::NamedTuple
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    if isnothing(plotconfig.groupBy)
        return
    end

    groupColors = getGroupColorMap(model, plotconfig)
    for (i, group1) in enumerate(keys(groupColors))
        for (j, group2) in enumerate(keys(groupColors))
            if i < j
                p = plot(
                    leg=plotconfig.leg,
                    aspect_ratio=:equal,
                    margin=plotconfig.margin,
                    yrotation=90,
                    size=(plotconfig.size, plotconfig.size)
                )

                letter = DEFAULT_ALPHABET[length(ps)+1]
                title!(p, "($(letter)) $(group1) vs. $(group2)")
                groupRows = (model.metadata[!, plotconfig.groupBy] .== group1) .|
                            (model.metadata[!, plotconfig.groupBy] .== group2)
                vals = map(model.metadata[!, plotconfig.groupBy]) do g
                    if g == group1
                        return 1
                    elseif g == group2
                        return 2
                    else
                        return missing
                    end
                end

                # plot_network!(M, p, model, TrendLinearEdgePainter(vals, groupColors[group1], groupColors[group2]), plotconfig, groupRows)
                plot_network!(M, p, model, TrendLinearEdgePainter(vals, plotconfig.negColor, plotconfig.posColor), plotconfig, groupRows)
                plot_units!(M, p, model, plotconfig, groupRows)
                plot_means!(M, p, model, plotconfig, groupRows)
                plot!(p,
                    [-999], [-999],
                    label="$(group2)-leaning Edge",
                    seriestype=:scatter,
                    markershape=:hline,
                    markersize=GLOBAL_UNIT_SIZE,
                    markercolor=[plotconfig.posColor],
                    markerstrokewidth=1
                )
                plot!(p,
                    [-999], [-999],
                    label="$(group1)-leaning Edge",
                    seriestype=:scatter,
                    markershape=:hline,
                    markersize=GLOBAL_UNIT_SIZE,
                    markercolor=[plotconfig.negColor],
                    markerstrokewidth=1
                )
                push!(ps, p)
            end
        end
    end
end

function plot_extras!(
        ::Type{M}, ps::Array{Plot}, model::AbstractLinearENAModel, plotconfig::NamedTuple
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    if !plotconfig.showExtras
        return
    end

    if !isnothing(plotconfig.trajectoryBy)
        p = plot(
            leg=plotconfig.leg,
            aspect_ratio=:equal,
            margin=plotconfig.margin,
            yrotation=90,
            size=(plotconfig.size, plotconfig.size)
        )

        letter = DEFAULT_ALPHABET[length(ps)+1]
        title!(p, "($(letter)) Trajectory by $(plotconfig.trajectoryBy)")
        plot_network!(M, p, model, NodesOnlyLinearEdgePainter(), plotconfig)
        plot_units!(M, p, model, merge(
            plotconfig,
            (
                spectralColorBy=plotconfig.trajectoryBy,
                groupBy=nothing
            )
        ))

        plot_trajectories!(M, p, model, plotconfig)
        push!(ps, p)
    end

    if !isnothing(plotconfig.spectoryBy)
        p = plot(
            leg=plotconfig.leg,
            aspect_ratio=:equal,
            margin=plotconfig.margin,
            yrotation=90,
            size=(plotconfig.size, plotconfig.size)
        )

        letter = DEFAULT_ALPHABET[length(ps)+1]
        title!(p, "($(letter)) Spectory by $(plotconfig.spectoryBy)")
        plot_network!(M, p, model, NodesOnlyLinearEdgePainter(), plotconfig)
        plot_units!(M, p, model, merge(
            plotconfig,
            (
                spectralColorBy=plotconfig.spectoryBy,
                groupBy=nothing
            )
        ))

        plot_spectories!(M, p, model, plotconfig)
        push!(ps, p)
    end    
end

# Plot Elements
function plot_network!(
        ::Type{M},
        p::Plot,
        model::AbstractLinearENAModel,
        edgePainter::AbstractLinearEdgePainter,
        plotconfig::NamedTuple,
        displayRows::BitVector=BitVector(repeat([true], nrow(model.accum)))
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    if !plotconfig.showNetworks
        return
    end

    # Just pass the work off to the edge painter
    paint_edges!(p, edgePainter, model, plotconfig, displayRows)
end

function plot_units!(
        ::Type{M},
        p::Plot,
        model::AbstractLinearENAModel,
        plotconfig::NamedTuple,
        displayRows::BitVector=BitVector(repeat([true], nrow(model.accum)))
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    if !plotconfig.showUnits
        return
    end

    # Get unit plotted positions
    unitIDs = model.accum.unitID
    xs = [fixX(x, plotconfig) for x in model.points[plotconfig.x, unitIDs]]
    ys = [fixY(y, plotconfig) for y in model.points[plotconfig.y, unitIDs]]
    
    # Default color: black
    colors = colorant"black"
    label = plotconfig.unitLabel

    # Optionally, color spectrally
    if !isnothing(plotconfig.spectralColorBy) && isnothing(plotconfig.groupBy)
        # label = string(plotconfig.spectralColorBy) # TODO show a scale, not one random color
        label = nothing
        colData = nothing
        if plotconfig.spectralColorBy in Symbol.(names(model.accum))
            colData = model.accum
        elseif plotconfig.spectralColorBy in Symbol.(names(model.metadata))
            colData = model.metadata
        end

        if !isnothing(colData)
            colVals = Vector{Real}(colData[displayRows, plotconfig.spectralColorBy])
            allColVals = Vector{Real}(colData[!, plotconfig.spectralColorBy])
            colMin = minimum(allColVals)
            colMax = maximum(allColVals)
            colVals = colVals .- colMin
            colVals /= colMax
            colors = getSpectralColors(colVals)

            if plotconfig.colorbar
                T = 10
                plot!(p,
                    [], [], label=false,
                    marker_z=colMin:colMax,
                    color=cgrad(getSpectralColors((t - 1)/(T - 1) for t in 1:T)),
                    colorbar=true,
                    colorbar_title=string(plotconfig.spectralColorBy)
                )
            end
        end
    end
    
    # Draw the units: if grouping, have to plot each group individually
    if !isnothing(plotconfig.groupBy)
        for group in unique(model.metadata[displayRows, plotconfig.groupBy])
            groupRows = model.metadata[!, plotconfig.groupBy] .== group
            plotRows = displayRows .& groupRows
            label = string(group, " ", plotconfig.unitLabel)
            groupColors = getGroupColorMap(model, plotconfig)
            plot!(
                p, xs[plotRows], ys[plotRows],
                label=label,
                seriestype=:scatter,
                markershape=:circle,
                markersize=GLOBAL_UNIT_SIZE,
                markercolor=groupColors[group],
                markerstrokewidth=0
            )
        end
    else
        plot!(
            p, xs[displayRows], ys[displayRows],
            label=label,
            seriestype=:scatter,
            markershape=:circle,
            markersize=GLOBAL_UNIT_SIZE,
            markercolor=colors,
            markerstrokewidth=0
        )
    end
end

function plot_means!(
        ::Type{M},
        p::Plot,
        model::AbstractLinearENAModel,
        plotconfig::NamedTuple,
        displayRows::BitVector=BitVector(repeat([true], nrow(model.accum)))
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    if !plotconfig.showMeans
        return
    end

    if isnothing(plotconfig.groupBy)
        return
    end

    groupColors = getGroupColorMap(model, plotconfig)
    for group in keys(groupColors)
        groupRows = model.metadata[!, plotconfig.groupBy] .== group
        meanedRows = displayRows .& groupRows
        if any(meanedRows)
            unitIDs = model.accum[meanedRows, :unitID]
            xs = [fixX(x, plotconfig) for x in model.points[plotconfig.x, unitIDs]]
            ys = [fixY(y, plotconfig) for y in model.points[plotconfig.y, unitIDs]]
            color = groupColors[group]
            if isnothing(plotconfig.innerGroupBy)
                plotConfidenceInterval(p, xs, ys, color, :square, "$(group) Mean", plotconfig.confidenceShape)
            else
                allInnerGroups = sort(unique(model.metadata[!, plotconfig.innerGroupBy]))
                # NOTE when we run out of markers, we stop plotting inner groups
                for (igroup, imarker) in zip(allInnerGroups, DEFAULT_INNER_MEAN_MARKERS)
                    igroupRows = model.metadata[!, plotconfig.innerGroupBy] .== igroup
                    imeanedRows = meanedRows .& igroupRows
                    if any(imeanedRows)
                        iunitIDs = model.accum[imeanedRows, :unitID]
                        ixs = [fixX(x, plotconfig) for x in model.points[plotconfig.x, iunitIDs]]
                        iys = [fixY(y, plotconfig) for y in model.points[plotconfig.y, iunitIDs]]
                        plotConfidenceInterval(p, ixs, iys, color, imarker, "$(group) Mean where $(plotconfig.innerGroupBy) = $(igroup)", plotconfig.confidenceShape)
                    end
                end
            end
        end
    end
end

function plot_trajectories!(
        ::Type{M},
        p::Plot,
        model::AbstractLinearENAModel,
        plotconfig::NamedTuple,
        displayRows::BitVector=BitVector(repeat([true], nrow(model.accum)))
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    if isnothing(plotconfig.trajectoryBy)
        return
    end

    tempPositions = DataFrame(Dict(
        :unitID => model.accum.unitID,
        :pos_x => [fixX(model.points[plotconfig.x, unitID], plotconfig) for unitID in model.accum.unitID],
        :pos_y => [fixY(model.points[plotconfig.y, unitID], plotconfig) for unitID in model.accum.unitID]
    ))

    smoothingData = innerjoin(model.accum[displayRows, :], model.metadata[displayRows, :], on=:unitID)
    smoothingData = innerjoin(smoothingData, tempPositions, on=:unitID)
    if plotconfig.trajectoryBy in Symbol.(names(smoothingData))
        tb = plotconfig.trajectoryBy
        if length(unique(smoothingData[!, plotconfig.trajectoryBy])) > plotconfig.trajectoryBins
            tb = :TIEDRANK
            ranks = tiedrank(smoothingData[!, plotconfig.trajectoryBy])
            ranks /= maximum(ranks) + 1
            ranks *= plotconfig.trajectoryBins
            smoothingData[!, tb] = floor.(ranks)
        end

        smoothingData = combine(
            groupby(smoothingData, tb),
            tb => median => tb,
            :pos_x => median => :pos_x,
            :pos_y => median => :pos_y
        )
            
        if nrow(smoothingData) > 3
            smoothingData = sort(smoothingData, tb)
            ts = smoothingData[!, tb]
            ps = transpose(Matrix{Float64}(smoothingData[!, [:pos_x, :pos_y]]))
            bspline = ParametricSpline(ts, ps, k=3, bc="nearest")
            smooth_ts = range(ts[1], stop=ts[end], length=500)
            smooth_ps = transpose(bspline(smooth_ts))
            plot!(
                p,
                smooth_ps[:, 1],
                smooth_ps[:, 2],
                linecolor=:black,
                arrow=:closed,
                label=nothing
            )
        end
    end 
end

function plot_spectories!(
        ::Type{M},
        p::Plot,
        model::AbstractLinearENAModel,
        plotconfig::NamedTuple,
        displayRows::BitVector=BitVector(repeat([true], nrow(model.accum)))
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    if isnothing(plotconfig.spectoryBy)
        return
    end

    tempPositions = DataFrame(Dict(
        :unitID => model.accum.unitID,
        :pos_x => [fixX(model.points[plotconfig.x, unitID], plotconfig) for unitID in model.accum.unitID],
        :pos_y => [fixY(model.points[plotconfig.y, unitID], plotconfig) for unitID in model.accum.unitID]
    ))

    smoothingData = innerjoin(model.accum[displayRows, :], model.metadata[displayRows, :], on=:unitID)
    smoothingData = innerjoin(smoothingData, tempPositions, on=:unitID)
    if plotconfig.spectoryBy in Symbol.(names(smoothingData))
        sb = plotconfig.spectoryBy
        sort!(smoothingData, sb)
        sb_min = first(smoothingData[!, sb])
        sb_max = last(smoothingData[!, sb])
        T = floor(Int, 2 * (1 / plotconfig.spectoryPercent) - 1)
        for t in 1:T
            left = floor(Int, (t - 1) / T * nrow(smoothingData) + 1)
            right = min(floor(Int, t / T * nrow(smoothingData) + 1), nrow(smoothingData))
            mid = round(Int, (left+right) / 2)
            sb_mid = smoothingData[mid, sb]
            xs = Vector(smoothingData[left:right, :pos_x])
            ys = Vector(smoothingData[left:right, :pos_y])
            color = first(getSpectralColors([(sb_mid - sb_min)/(sb_max - sb_min)]))
            plot_kde_without_cbar!(p,
                kde((xs, ys)),
                levels=[1],
                color=color
            )
        end
    end 
end

# BUGFIX existing contour plots have a cbar by default.
# That existing cbar conflicts with other cbars we might want to do, eg., one
# based on the unit points themselves, subtraction network edges, etc
function plot_kde_without_cbar!(p::Plot, U::BivariateKDE; levels::Vector{<:Real}=[1], color::Colorant=:black)
    for level in Contour.levels(contours(U.x, U.y, U.density, levels))
        for line in lines(level)

            plot!(p, vertices(line), linecolor=color, label=nothing)
        end
    end
end