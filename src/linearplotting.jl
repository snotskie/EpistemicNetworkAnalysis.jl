# Defaults
const DEFAULT_NEG_COLOR = colorant"#cc423a"
const DEFAULT_POS_COLOR = colorant"#218ebf"
const DEFAULT_EXTRA_COLORS = [
    colorant"#56BD7C", colorant"#EF691B", colorant"#9d5dbb",
    colorant"#FBC848", colorant"#D0386C",
    colorant"#f18e9f", colorant"#9A9EAB", colorant"#ff8c39",
    colorant"#346B88"
]

const DEFAULT_ALPHABET = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZθωερτψυπασδφγηςκλζχξωβνμΘΩΨΥΠΣΔΦΓΛΞΩ"
const GLOBAL_MAX_EDGE_SIZE = 8
const GLOBAL_MAX_NODE_SIZE = 8
const GLOBAL_UNIT_SIZE = 3

function defaultplotkwargs(
        ::Type{M},
        model::AbstractLinearENAModel;
        x::Int=1,
        y::Int=2,
        margin::PlotMeasures.AbsoluteLength=10mm,
        size::Real=600,
        meanCenter::Bool=model.config.sphereNormalize,
        origin::Array{<:Real}=(meanCenter ?  [mean(model.points[x, :]), mean(model.points[y, :])] : [0,0]),
        lims::Real=1,
        xticks::Array{<:Real}=round.([origin[1]-lims, origin[1], origin[1]+lims], digits=4),
        yticks::Array{<:Real}=round.([origin[2]-lims, origin[2], origin[2]+lims], digits=4),
        xlims::Array{<:Real}=xticks[[1, end]],
        ylims::Array{<:Real}=yticks[[1, end]],
        titles::Array{<:AbstractString}=String[],
        xlabel::AbstractString=model.embedding[x, :label],
        ylabel::AbstractString=model.embedding[y, :label],
        unitLabel::AbstractString="Unit",
        leg::Union{Symbol,Bool}=:topleft,
        negColor::Colorant=DEFAULT_NEG_COLOR,
        posColor::Colorant=DEFAULT_POS_COLOR,
        extraColors::Array{<:Colorant,1}=DEFAULT_EXTRA_COLORS,
        alphabet::String=DEFAULT_ALPHABET,
        flipX::Bool=false,
        flipY::Bool=false,
        groupBy::Union{Symbol,Nothing}=nothing,
        spectralColorBy::Union{Symbol,Nothing}=nothing,
        showExtras::Bool=false,
        showNetworks::Bool=true,
        showUnits::Bool=true,
        showMeans::Bool=true,
        showWarps::Bool=false,
        fitNodesToCircle::Bool=false,
        weakLinks::Bool=true,
        # TODO the options I cut from the model type
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
        extraColors=extraColors,
        alphabet=alphabet,
        flipX=flipX,
        flipY=flipY,
        groupBy=groupBy,
        spectralColorBy=spectralColorBy,
        showExtras=showExtras,
        showNetworks=showNetworks,
        showUnits=showUnits,
        showMeans=showMeans,
        showWarps=showWarps,
        fitNodesToCircle=fitNodesToCircle,
        weakLinks=weakLinks,
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
    # plot_extras!(M, ps, model, plotconfig)
    return layoutSubplots(ps, plotconfig)
end

# Helpers
function fixX(x::Real, model::AbstractLinearENAModel, plotconfig::NamedTuple)
    # NOTE mean centering happens in the plotkwargs via xticks etc.
    # if plotconfig.meanCenter
    #     x -= mean(model.points[plotconfig.x, :])
    # end

    return x * (plotconfig.flipX ? -1 : 1)
end

function fixY(y::Real, model::AbstractLinearENAModel, plotconfig::NamedTuple)
    # NOTE mean centering happens in the plotkwargs via xticks etc.
    # if plotconfig.meanCenter
    #     y -= mean(model.points[plotconfig.y, :])
    # end

    return y * (plotconfig.flipY ? -1 : 1)
end

function fixPoint(point::Array{<:Real}, model::AbstractLinearENAModel, plotconfig::NamedTuple)
    return [
        fixX(point[1], model, plotconfig),
        fixY(point[2], model, plotconfig)
    ]
end

function getGroupColorMap(model::AbstractLinearENAModel, plotconfig::NamedTuple)
    allGroups = sort(unique(model.metadata[!, plotconfig.groupBy]))
    groupColors = Dict(
        g => plotconfig.extraColors[i]
        for (i, g) in enumerate(allGroups)
        if i <= length(plotconfig.extraColors) # in case we don't have enough colors
    )

    return DefaultDict(colorant"black", groupColors) # in case we don't have enough colors  
end

function plotConfidenceInterval(p, xs, ys, color, shape, label)

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
        ], model, plotconfig))

        if edge.kind == :count
            push!(edgePoints[edge.edgeID], fixPoint([0, 0], model, plotconfig))
        elseif edge.kind in [:undirected, :undirected]
            if plotconfig.showWarps
                pointT = fixPoint([
                    model.embedding[plotconfig.x, edge.edgeID],
                    model.embedding[plotconfig.y, edge.edgeID]
                ], model, plotconfig)
                
                push!(edgePoints[edge.edgeID], pointT)
                push!(edgePoints[edge.edgeID], pointT)
                push!(edgePoints[edge.edgeID], pointT)
            end

            push!(edgePoints[edge.edgeID], fixPoint([
                model.pointsNodes[plotconfig.x, edge.response],
                model.pointsNodes[plotconfig.y, edge.response]
            ], model, plotconfig))
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
        if edgeWidths[edgeID] > 0
            plot!(p,
                first.(edgePoints[edgeID]), # xs
                last.(edgePoints[edgeID]), # ys
                label=nothing,
                seriestype=:curves,
                arrows=(edge.kind == :directed),
                linewidth=edgeWidths[edgeID],
                linecolor=edgeColors[edgeID]
            )
        end
    end

    for nodeID in model.nodes.nodeID
        x = fixX(model.pointsNodes[plotconfig.x, nodeID], model, plotconfig)
        y = fixY(model.pointsNodes[plotconfig.y, nodeID], model, plotconfig)
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
    vals::Vector{<:Real}
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
    mid_color = weighted_color_mean(0.5, RGB(edgePainter.neg_color), RGB(edgePainter.pos_color))
    mid_color = weighted_color_mean(0.3, RGB(mid_color), colorant"white")
    edgeColorMap = nonlinearGradientMap(
        weighted_color_mean(0.95, edgePainter.neg_color, colorant"black"),
        mid_color,
        weighted_color_mean(0.95, edgePainter.pos_color, colorant"black"),
        curve=2.5
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

            if !plotconfig.weakLinks && abs(pearson) < 0.3
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
end

# Subplots
function plot_omnibus!(
        ::Type{M}, ps::Array{Plot}, model::AbstractLinearENAModel, plotconfig::NamedTuple
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    p = plot(;
        leg=plotconfig.leg,
        margin=plotconfig.margin,
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
        margin=plotconfig.margin,
        size=(plotconfig.size, plotconfig.size)
    )

    py = plot(
        leg=plotconfig.leg,
        margin=plotconfig.margin,
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
            margin=plotconfig.margin,
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
                    margin=plotconfig.margin,
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

                plot_network!(M, p, model, TrendLinearEdgePainter(vals, groupColors[group1], groupColors[group2]), plotconfig, groupRows)
                plot_units!(M, p, model, plotconfig, groupRows)
                plot_means!(M, p, model, plotconfig, groupRows)
                push!(ps, p)
            end
        end
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

    # TODO move to plot_extras or _trends or _groups or idk?
    # #### Optional: illustrate a trajectory by a continuous, non-repeating value
    # if !isnothing(showTrajectoryBy)
    #     smoothingData = innerjoin(ena.accumModel[displayRows, :], ena.metadata[displayRows, :], on=:unitID)
    #     if showTrajectoryBy in Symbol.(names(smoothingData))
    #         smoothingData = combine(
    #             groupby(smoothingData, showTrajectoryBy),
    #             showTrajectoryBy => mean => showTrajectoryBy,
    #             :pos_x => mean => :pos_x,
    #             :pos_y => mean => :pos_y
    #         )
                
    #         if nrow(smoothingData) > 3
    #             smoothingData = sort(smoothingData, showTrajectoryBy)
    #             ts = smoothingData[!, showTrajectoryBy]
    #             ps = transpose(Matrix{Float64}(smoothingData[!, [:pos_x, :pos_y]]))
    #             bspline = ParametricSpline(ts, ps, k=3, bc="nearest")
    #             smooth_ts = range(ts[1], stop=ts[end], length=500)
    #             smooth_ps = transpose(bspline(smooth_ts))
    #             plot!(p, smooth_ps[:, 1] * (flipX ? -1 : 1), smooth_ps[:, 2] * (flipY ? -1 : 1), linecolor=:black, arrow=:closed)
    #         end
    #     end 
    # end
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
    xs = [fixX(x, model, plotconfig) for x in model.points[plotconfig.x, unitIDs]]
    ys = [fixY(y, model, plotconfig) for y in model.points[plotconfig.y, unitIDs]]
    
    # Default color: black
    colors = colorant"black"
    label = plotconfig.unitLabel

    # Optionally, color spectrally
    if !isnothing(plotconfig.spectralColorBy)
        label = string(plotconfig.spectralColorBy)
        if plotconfig.spectralColorBy in Symbol.(names(model.accum))
            colVals = Vector{Real}(model.accum[displayRows, plotconfig.spectralColorBy])
            allColVals = colVals = Vector{Real}(model.accum[!, plotconfig.spectralColorBy])
            colVals = colVals .- minimum(allColVals)
            colVals /= maximum(allColVals)
            colors = [HSL(colVal*240, 1, 0.5) for colVal in colVals]
        elseif plotconfig.spectralColorBy in Symbol.(names(model.metadata))
            colVals = Vector{Real}(model.metadata[displayRows, plotconfig.spectralColorBy])
            allColVals = Vector{Real}(model.metadata[!, plotconfig.spectralColorBy])
            colVals = colVals .- minimum(allColVals)
            colVals /= maximum(allColVals)
            colors = [HSL(colVal*240, 1, 0.5) for colVal in colVals]
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
        unitIDs = model.accum[meanedRows, :unitID]
        xs = [fixX(x, model, plotconfig) for x in model.points[plotconfig.x, unitIDs]]
        ys = [fixX(y, model, plotconfig) for y in model.points[plotconfig.y, unitIDs]]
        color = groupColors[group]
        plotConfidenceInterval(p, xs, ys, color, :square, "$(group) Mean")
    end
end

# TODO below here

# ## Text display
# function Base.display(ena::AbstractENAModel) # TODO use show, not display

#     ### Show plotted points
#     println("Units (plotted points):")
#     show(ena.accumModel, allrows=true)
#     println()

#     ### Show centroids
#     println("Units (centroids):")
#     show(ena.centroidModel, allrows=true)
#     println()

#     ### Show codes
#     println("Codes:")
#     show(ena.codeModel, allrows=true)
#     println()

#     ### Show network
#     println("Network:")
#     show(ena.networkModel, allrows=true)
#     println()

#     ### Show every test result we have
#     results = test(ena)
#     for key in keys(results)
#         println("$key:")
#         println(results[key])
#         println()
#     end
# end


## Plotting
### Common defaults and globals


# ### Top-level wrapper
# function plot(ena::AbstractENAModel;
# margin=10mm, size=600, lims=1, ticks=[-lims, 0, lims],
# titles=[], xlabel="X", ylabel="Y", leg=:topleft,
# negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
# extraColors::Array{<:Colorant,1}=DEFAULT_EXTRA_COLORS,
# flipX=false, flipY=false,
# singleUnit=nothing, groupBy=nothing,
# showExtras::Bool=false, showNetworks::Bool=true, showUnits::Bool=true, showMeans::Bool=true,
# kwargs...)

# #### Combine the kwargs to make them easier to pass without needing
# #### knowledge of what my child functions need
# kwargs = (
#     margin=margin, size=size, lims=lims, ticks=ticks,
#     titles=titles, xlabel=xlabel, ylabel=ylabel, leg=leg,
#     negColor=negColor, posColor=posColor, extraColors=extraColors,
#     flipX=flipX, flipY=flipY, singleUnit=singleUnit, groupBy=groupBy,
#     kwargs...
# )

# #### Initialize usual subplots
# ps = [
#     plot(leg=leg, margin=margin, size=(size, size)), # omnibus
#     plot(leg=false, margin=margin, size=(size, size)),  # predictive x
#     plot(leg=false, margin=margin, size=(size, size))  # predictive y
# ]

# #### Start usual subplots: Distribution
# title!(ps[1], "(a) " * get(titles, 1, "Distribution"))

# #### Draw usual subplots: Dynamics
# title!(ps[2], "(b) " * get(titles, 2, "Rate of Change by X"))
# plot_predictive!(ps[2], ena, :pos_x; kwargs...)
# title!(ps[3], "(c) " * get(titles, 3, "Rate of Change by Y"))
# plot_predictive!(ps[3], ena, :pos_y; kwargs...)

# #### If we need group-wise subplots...
# if !isnothing(groupBy)

#     #### ...continue usual subplots: Distribution
#     allRows = [true for row in eachrow(ena.metadata)]
#     if showNetworks
#         plot_network!(ps[1], ena, allRows; kwargs...)
#     end

#     #### ...then for each group...
#     groups = sort(unique(ena.metadata[!, groupBy]))
#     letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZθωερτψυπασδφγηςκλζχξωβνμΘΩΨΥΠΣΔΦΓΛΞΩ"[(1+length(ps)):end]
#     letterIndices = collect(eachindex(letters)) # because we're possibly dealing with unicode
#     for (g, group) in enumerate(groups)

#         #### Initialize and draw them, until we run out of letters
#         if g <= length(extraColors) && g <= length(letterIndices)
#             p = plot(leg=false, margin=margin, size=(size, size))
#             groupRows = [row[groupBy] == group for row in eachrow(ena.metadata)]
#             if showNetworks
#                 plot_network!(p, ena, groupRows; color=extraColors[g], kwargs...)
#             end

#             if showUnits
#                 plot_units!(p, ena, groupRows; unitLabel="$(group) Units", color=extraColors[g], kwargs...)
#                 plot_units!(ps[1], ena, groupRows; unitLabel="$(group) Units", color=extraColors[g], kwargs...)
#             end

#             if showExtras
#                 plot_extras!(p, ena, groupRows; color=extraColors[g], kwargs...)
#             end
            
#             if showMeans
#                 plot_cis!(p, ena, groupRows, group; color=extraColors[g], kwargs...)
#                 plot_cis!(ps[1], ena, groupRows, group; color=extraColors[g], kwargs...)
#             end

#             title!(p, "($(letters[letterIndices[g]])) " * string(get(titles, 2+g, group)))
#             push!(ps, p)
#         end
#     end

#     #### ...then for each pair of groups...
#     n = length(groups)+1
#     for i in 1:length(groups)
#         for j in (i+1):length(groups)
#             if n <= length(letterIndices) && j <= length(extraColors)
#                 posGroup = groups[j]
#                 negGroup = groups[i]
#                 p = plot(leg=false, margin=margin, size=(size, size))
#                 plot_subtraction!(p, ena, groupBy, groups[i], groups[j];
#                     kwargs..., negColor=extraColors[i], posColor=extraColors[j])

#                 posGroupRows = [row[groupBy] == posGroup for row in eachrow(ena.metadata)]
#                 negGroupRows = [row[groupBy] == negGroup for row in eachrow(ena.metadata)]
#                 if showUnits
#                     plot_units!(p, ena, posGroupRows; color=extraColors[j], kwargs...)
#                     plot_units!(p, ena, negGroupRows; color=extraColors[i], kwargs...)
#                 end

#                 if showMeans
#                     plot_cis!(p, ena, posGroupRows, posGroup; color=extraColors[j], kwargs...)
#                     plot_cis!(p, ena, negGroupRows, negGroup; color=extraColors[i], kwargs...)
#                 end

#                 temp = "$(groups[i]) vs. $(groups[j])"
#                 title!(p, "($(letters[letterIndices[n]])) " * string(get(titles, 2+n, temp)))
#                 push!(ps, p)
#                 n += 1
#             end
#         end
#     end

#     #### then finish usual subplots: Distribution
#     if showExtras
#         plot_extras!(ps[1], ena, allRows; kwargs...)
#     end
# else
#     #### Else just draw usual subplots: Distribution
#     allRows = [true for row in eachrow(ena.metadata)]
#     if showNetworks
#         plot_network!(ps[1], ena, allRows; kwargs...)
#     end

#     if showUnits
#         plot_units!(ps[1], ena, allRows; kwargs...)
#     end

#     if showExtras
#         plot_extras!(ps[1], ena, allRows; kwargs...)
#     end
# end

# #### Layout the subplots
# results = test(ena)
# for p in ps

    
#     if !isnan(results[:variance_x])
#         xlabel!(p, "$xlabel ($(round(Int, results[:variance_x]*100))%)")
#     else
#         xlabel!(p, xlabel)
#     end

#     if !isnan(results[:variance_y])
#         ylabel!(p, "$ylabel ($(round(Int, results[:variance_y]*100))%)")
#     else
#         ylabel!(p, ylabel)
#     end
# end


# end

# ### Helper - Draw the dots
# function plot_units!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
# color::Colorant=colorant"black", spectralColorBy::Union{Symbol,Nothing}=nothing,
# flipX::Bool=false, flipY::Bool=false,
# unitLabel::String="Units",
# kwargs...)

# #### Get the x/y positions
# xs, ys = help_xs_and_ys(ena, displayRows, flipX, flipY)

# #### Optional: color code by a continuous value
# if !isnothing(spectralColorBy)
#     if spectralColorBy in Symbol.(names(ena.accumModel))
#         colVals = Vector{Float64}(ena.accumModel[displayRows, spectralColorBy])
#         allColVals = colVals = Vector{Float64}(ena.accumModel[!, spectralColorBy])
#         colVals = colVals .- minimum(allColVals)
#         colVals /= maximum(allColVals)
#         color = [HSL(colVal*240, 1, 0.5) for colVal in colVals]
#     elseif spectralColorBy in Symbol.(names(ena.metadata))
#         colVals = Vector{Float64}(ena.metadata[displayRows, spectralColorBy])
#         allColVals = Vector{Float64}(ena.metadata[!, spectralColorBy])
#         colVals = colVals .- minimum(allColVals)
#         colVals /= maximum(allColVals)
#         color = [HSL(colVal*240, 1, 0.5) for colVal in colVals]
#     end
# end

# #### Draw the units, in black by default
# plot!(p, xs, ys,
#     label=unitLabel,
#     seriestype=:scatter,
#     markershape=:circle,
#     markersize=GLOBAL_UNIT_SIZE,
#     markercolor=color,
#     markerstrokewidth=0)
# end

# ### Helper - Draw the lines
# function plot_network!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
# color::Colorant=colorant"#aaa",
# flipX::Bool=false, flipY::Bool=false, showWarps::Bool=false, showCodeLabels::Bool=true,
# showArrows::Bool=false, showTrajectoryBy::Union{Symbol,Nothing}=nothing, showNetworkLines::Bool=true,
# rotateCodeLabels::Bool=false,
# kwargs...)

# #### Find the true weight on each line
# allLineWidths = map(eachrow(ena.networkModel)) do networkRow
#     return sum(ena.accumModel[!, networkRow[:relationship]])
# end

# displayAccums = ena.accumModel[displayRows, :]
# lineWidths = map(eachrow(ena.networkModel)) do networkRow
#     return sum(displayAccums[!, networkRow[:relationship]])
# end

# #### Rescale the lines
# lineWidths *= GLOBAL_MAX_EDGE_SIZE / maximum(allLineWidths)

# #### Initialize code widths, compute while we visit each line
# codeWidths = zeros(nrow(ena.codeModel))

# #### For each line...
# for (i, networkRow) in enumerate(eachrow(ena.networkModel))
#     j, k = ena.relationshipMap[networkRow[:relationship]]

#     #### ...add to its code weights
#     codeWidths[j] += lineWidths[i]
#     codeWidths[k] += lineWidths[i]

#     #### and plot that line
#     if showNetworkLines
#         pointA = [ena.codeModel[j, :pos_x] * (flipX ? -1 : 1), ena.codeModel[j, :pos_y] * (flipY ? -1 : 1)]
#         pointB = [ena.codeModel[k, :pos_x] * (flipX ? -1 : 1), ena.codeModel[k, :pos_y] * (flipY ? -1 : 1)]
#         if j == k
#             pointB = [0.0, 0.0]
#         end
        
#         pointT = (pointA+pointB)/2
#         if showWarps
#             pointT = [networkRow[:weight_x] * (flipX ? -1 : 1), networkRow[:weight_y] * (flipY ? -1 : 1)]
#         end

#         points = hcat(pointA, pointT, pointT, pointT, pointB)
#         # arrows = nothing
#         # if showArrows
#         #     arrows = true #lineWidths[i] #arrow(:closed, :head, lineWidths[i], lineWidths[i])
#         # end

#         plot!(p,
#             points[1, :],
#             points[2, :],
#             label=nothing,
#             seriestype=:curves,
#             arrows=showArrows,
#             linewidth=lineWidths[i],
#             linecolor=color)
#     end
# end

# #### Rescale and draw the codes
# codeWidths *= GLOBAL_MAX_NODE_SIZE / maximum(codeWidths)
# x = ena.codeModel[!, :pos_x] * (flipX ? -1 : 1)
# y = ena.codeModel[!, :pos_y] * (flipY ? -1 : 1)
# if showCodeLabels
#     labels = map(zip(ena.codeModel[!, :code], x, y)) do (label, xi, yi)
#         if rotateCodeLabels
#             return rotatedLabel(label, xi, yi)
#         else
#             return text(label, :top, default(:xtickfontsize))
#         end
#     end

#     plot!(p, x, y,
#         label=nothing,
#         seriestype=:scatter,
#         series_annotations=labels,
#         markershape=:circle,
#         markersize=codeWidths,
#         markercolor=:black,
#         markerstrokewidth=0)
# else
#     plot!(p, x, y,
#         label=nothing,
#         seriestype=:scatter,
#         # series_annotations=labels,
#         markershape=:circle,
#         markersize=codeWidths,
#         markercolor=:black,
#         markerstrokewidth=0)
# end

# #### Optional: illustrate a trajectory by a continuous, non-repeating value
# if !isnothing(showTrajectoryBy)
#     smoothingData = innerjoin(ena.accumModel[displayRows, :], ena.metadata[displayRows, :], on=:unitID)
#     if showTrajectoryBy in Symbol.(names(smoothingData))
#         smoothingData = combine(
#             groupby(smoothingData, showTrajectoryBy),
#             showTrajectoryBy => mean => showTrajectoryBy,
#             :pos_x => mean => :pos_x,
#             :pos_y => mean => :pos_y
#         )
            
#         if nrow(smoothingData) > 3
#             smoothingData = sort(smoothingData, showTrajectoryBy)
#             ts = smoothingData[!, showTrajectoryBy]
#             ps = transpose(Matrix{Float64}(smoothingData[!, [:pos_x, :pos_y]]))
#             bspline = ParametricSpline(ts, ps, k=3, bc="nearest")
#             smooth_ts = range(ts[1], stop=ts[end], length=500)
#             smooth_ps = transpose(bspline(smooth_ts))
#             plot!(p, smooth_ps[:, 1] * (flipX ? -1 : 1), smooth_ps[:, 2] * (flipY ? -1 : 1), linecolor=:black, arrow=:closed)
#         end
#     end 
# end
# end

# ### Helper - Draw the predictive lines
# function plot_predictive!(p::Plot, ena::AbstractENAModel, targetCol::Symbol;
# negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
# flipX::Bool=false, flipY::Bool=false, weakLinks::Bool=true, showWarps::Bool=false,
# showCodeLabels::Bool=true, showArrows::Bool=false, reverseLineSort::Bool=false,
# rotateCodeLabels::Bool=false,
# kwargs...)

# ### Grab the data we need as one data frame
# regressionData = hcat(ena.accumModel, ena.metadata, makeunique=true)
# xs, ys = help_xs_and_ys(ena, !, flipX, flipY)
# regressionData[!, :pos_x] = xs
# regressionData[!, :pos_y] = ys

# ### Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
# for networkRow in eachrow(ena.networkModel)
#     regressionData[!, networkRow[:relationship]] = map(Float64, regressionData[!, networkRow[:relationship]])
# end

# ### Compute line widths as the strength (slope) between the xpos and the accum network weights
# f1 = FormulaTerm(term(:y), term(targetCol))
# lineData = map(eachrow(ena.networkModel)) do networkRow
#     r = networkRow[:relationship]
#     f1 = FormulaTerm(term(r), f1.rhs)
#     try
#         m1 = fit(LinearModel, f1, regressionData)
#         slope = coef(m1)[2]
#         pearson = cor(regressionData[!, targetCol], regressionData[!, r])
#         return (slope, pearson)
#     catch e
#         return (0, 0)
#     end
# end

# ### Color the lines based on their correlation with the x position
# midColor = weighted_color_mean(0.5, RGB(negColor), RGB(posColor))
# midColor = weighted_color_mean(0.3, RGB(midColor), colorant"white")
# lineColorMap = help_nonlinear_gradient(weighted_color_mean(0.95, negColor, colorant"black"),
#                                        midColor,
#                                        weighted_color_mean(0.95, posColor, colorant"black"),
#                                        curve=2.5)
# lineColors = map(lineData) do (slope, pearson)
#     if isnan(pearson)
#         pearson = 0
#     end

#     if flipX
#         pearson *= -1
#     end        

#     index = 1 + round(Int, (length(lineColorMap) - 1) * (pearson + 1) / 2)
#     return lineColorMap[index]
# end

# ### Size the lines based on their slope with the x position
# lineWidths = map(lineData) do (slope, pearson)
#     return abs(slope)
# end

# ### Normalize
# lineWidths *= GLOBAL_MAX_EDGE_SIZE / maximum(lineWidths)

# ### Placeholder, let's compute code weights as we visit each line
# codeWidths = zeros(nrow(ena.codeModel))
# codeVisible = Bool.(zeros(nrow(ena.codeModel)))

# ### For each line...
# networkData = hcat(ena.networkModel, DataFrame(:width => lineWidths, :color => lineColors, :pearson => last.(lineData)))
# for networkRow in sort(eachrow(networkData), by=row->row[:width], rev=reverseLineSort)

#     ### ...contribute to the code weights...
#     j, k = ena.relationshipMap[networkRow[:relationship]]
#     codeWidths[j] += networkRow[:width]
#     codeWidths[k] += networkRow[:width]

#     ### ...and if that line should be shown...
#     if weakLinks || abs(networkRow[:pearson]) >= 0.3
#         ### ...plot it in the right width and color
#         codeVisible[j] = true
#         codeVisible[k] = true
#         pointA = [ena.codeModel[j, :pos_x] * (flipX ? -1 : 1), ena.codeModel[j, :pos_y] * (flipY ? -1 : 1)]
#         pointB = [ena.codeModel[k, :pos_x] * (flipX ? -1 : 1), ena.codeModel[k, :pos_y] * (flipY ? -1 : 1)]
#         if j == k
#             pointB = [0.0, 0.0]
#         end
        
#         pointT = (pointA+pointB)/2
#         if showWarps
#             pointT = [networkRow[:weight_x] * (flipX ? -1 : 1), networkRow[:weight_y] * (flipY ? -1 : 1)]
#         end

#         points = hcat(pointA, pointT, pointT, pointT, pointB)
#         plot!(p,
#             points[1, :],
#             points[2, :],
#             label=nothing,
#             seriestype=:curves,
#             arrows=showArrows,
#             linewidth=networkRow[:width],
#             linecolor=networkRow[:color])
#     end
# end

# ### Rescale the code widths
# codeWidths *= GLOBAL_MAX_NODE_SIZE / maximum(codeWidths)

# ### And plot the codes and we're done
# x = ena.codeModel[codeVisible, :pos_x] * (flipX ? -1 : 1)
# y = ena.codeModel[codeVisible, :pos_y] * (flipY ? -1 : 1)
# if showCodeLabels
#     labels = map(zip(ena.codeModel[codeVisible, :code], x, y)) do (label, xi, yi)
#         if rotateCodeLabels
#             return rotatedLabel(label, xi, yi)
#         else
#             return text(label, :top, default(:xtickfontsize))
#         end
#     end

#     plot!(p, x, y,
#         label=nothing,
#         seriestype=:scatter,
#         series_annotations=labels,
#         markershape=:circle,
#         markersize=codeWidths,
#         markercolor=:black,
#         markerstrokewidth=0)
# else
#     plot!(p, x, y,
#         label=nothing,
#         seriestype=:scatter,
#         # series_annotations=labels,
#         markershape=:circle,
#         markersize=codeWidths,
#         markercolor=:black,
#         markerstrokewidth=0)
# end
# end

# ### Helper - Draw the subtraction lines (nearly identical to plot_predictive)
# function plot_subtraction!(p::Plot, ena::AbstractENAModel, groupVar::Symbol, negGroup::Any, posGroup::Any;
# negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
# flipX::Bool=false, flipY::Bool=false, weakLinks::Bool=true, showWarps::Bool=false,
# showCodeLabels::Bool=true, showArrows::Bool=false, reverseLineSort::Bool=false,
# rotateCodeLabels::Bool=false,
# kwargs...)

# ### Grab the data we need as one data frame
# regressionData = hcat(ena.accumModel, ena.metadata, makeunique=true)
# xs, ys = help_xs_and_ys(ena, !, flipX, flipY)
# regressionData[!, :pos_x] = xs
# regressionData[!, :pos_y] = ys

# ### Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
# for networkRow in eachrow(ena.networkModel)
#     regressionData[!, networkRow[:relationship]] = map(Float64, regressionData[!, networkRow[:relationship]])
# end

# regressionData[!, :SubtractionVar] = map(eachrow(regressionData)) do row
#     if row[groupVar] == posGroup
#         return 1
#     elseif row[groupVar] == negGroup
#         return 0
#     else
#         return missing
#     end
# end

# rowsForCor = map(eachrow(regressionData)) do row
#     return !ismissing(row[:SubtractionVar])
# end

# ### Compute line widths as the strength (slope) between the xpos and the accum network weights
# f1 = @formula(y ~ SubtractionVar)
# lineData = map(eachrow(ena.networkModel)) do networkRow
#     r = networkRow[:relationship]
#     f1 = FormulaTerm(term(r), f1.rhs)
#     try
#         m1 = fit(LinearModel, f1, regressionData)
#         slope = coef(m1)[2]
#         pearson = cor(regressionData[rowsForCor, :SubtractionVar], regressionData[rowsForCor, r])
#         return (slope, pearson)
#     catch e
#         return (0, 0)
#     end
# end

# ### Color the lines based on their correlation with the x position
# midColor = weighted_color_mean(0.5, RGB(negColor), RGB(posColor))
# midColor = weighted_color_mean(0.3, RGB(midColor), colorant"white")
# lineColorMap = help_nonlinear_gradient(weighted_color_mean(0.95, negColor, colorant"black"),
#                                        midColor,
#                                        weighted_color_mean(0.95, posColor, colorant"black"),
#                                        curve=1)

# lineColors = map(lineData) do (slope, pearson)
#     if isnan(pearson)
#         pearson = 0
#     end

#     if flipX
#         pearson *= -1
#     end        

#     index = 1 + round(Int, (length(lineColorMap) - 1) * (pearson + 1) / 2)
#     return lineColorMap[index]
# end

# ### Size the lines based on their slope with the x position
# lineWidths = map(lineData) do (slope, pearson)
#     return abs(slope)
# end

# ### Normalize
# lineWidths *= GLOBAL_MAX_EDGE_SIZE / maximum(lineWidths)

# ### Placeholder, let's compute code weights as we visit each line
# codeWidths = zeros(nrow(ena.codeModel))
# codeVisible = Bool.(zeros(nrow(ena.codeModel)))

# ### For each line...
# networkData = hcat(ena.networkModel, DataFrame(:width => lineWidths, :color => lineColors, :pearson => last.(lineData)))
# for networkRow in sort(eachrow(networkData), by=row->row[:width], rev=reverseLineSort)

#     ### ...contribute to the code weights...
#     j, k = ena.relationshipMap[networkRow[:relationship]]
#     codeWidths[j] += networkRow[:width]
#     codeWidths[k] += networkRow[:width]

#     ### ...and if that line should be shown...
#     if weakLinks || abs(networkRow[:pearson]) >= 0.3
#         ### ...plot it in the right width and color
#         codeVisible[j] = true
#         codeVisible[k] = true
#         pointA = [ena.codeModel[j, :pos_x] * (flipX ? -1 : 1), ena.codeModel[j, :pos_y] * (flipY ? -1 : 1)]
#         pointB = [ena.codeModel[k, :pos_x] * (flipX ? -1 : 1), ena.codeModel[k, :pos_y] * (flipY ? -1 : 1)]
#         if j == k
#             pointB = [0.0, 0.0]
#         end
        
#         pointT = (pointA+pointB)/2
#         if showWarps
#             pointT = [networkRow[:weight_x] * (flipX ? -1 : 1), networkRow[:weight_y] * (flipY ? -1 : 1)]
#         end

#         points = hcat(pointA, pointT, pointT, pointT, pointB)
#         plot!(p,
#             points[1, :],
#             points[2, :],
#             label=nothing,
#             seriestype=:curves,
#             arrows=showArrows,
#             linewidth=networkRow[:width],
#             linecolor=networkRow[:color])
#     end
# end

# ### Rescale the code widths
# codeWidths *= GLOBAL_MAX_NODE_SIZE / maximum(codeWidths)

# ### And plot the codes and we're done
# x = ena.codeModel[codeVisible, :pos_x] * (flipX ? -1 : 1)
# y = ena.codeModel[codeVisible, :pos_y] * (flipY ? -1 : 1)
# if showCodeLabels
#     labels = map(zip(ena.codeModel[codeVisible, :code], x, y)) do (label, xi, yi)
#         if rotateCodeLabels
#             return rotatedLabel(label, xi, yi)
#         else
#             return text(label, :top, default(:xtickfontsize))
#         end
#     end

#     plot!(p, x, y,
#         label=nothing,
#         seriestype=:scatter,
#         series_annotations=labels,
#         markershape=:circle,
#         markersize=codeWidths,
#         markercolor=:black,
#         markerstrokewidth=0)
# else
#     plot!(p, x, y,
#         label=nothing,
#         seriestype=:scatter,
#         # series_annotations=labels,
#         markershape=:circle,
#         markersize=codeWidths,
#         markercolor=:black,
#         markerstrokewidth=0)
# end


# # ### Placeholder, let's compute code weights as we visit each line
# # codeWidths = zeros(nrow(ena.codeModel))

# # ### For each line...
# # networkData = hcat(ena.networkModel, DataFrame(:width => lineWidths, :color => lineColors))
# # for networkRow in sort(eachrow(networkData), by=row->row[:width], rev=reverseLineSort)

# #     ### ...contribute to the code weights...
# #     j, k = ena.relationshipMap[networkRow[:relationship]]
# #     codeWidths[j] += networkRow[:width]
# #     codeWidths[k] += networkRow[:width]

# #     ### ...and plot that line, in the right width and color
# #     x = ena.codeModel[[j, k], :pos_x] * (flipX ? -1 : 1)
# #     y = ena.codeModel[[j, k], :pos_y] * (flipY ? -1 : 1)
# #     plot!(p, x, y,
# #         label=nothing,
# #         seriestype=:line,
# #         linewidth=networkRow[:width],
# #         linecolor=networkRow[:color])
# # end

# # ### Rescale the code widths
# # codeWidths *= GLOBAL_MAX_NODE_SIZE / maximum(codeWidths)

# # ### And plot the codes and we're done
# # x = ena.codeModel[!, :pos_x] * (flipX ? -1 : 1)
# # y = ena.codeModel[!, :pos_y] * (flipY ? -1 : 1)
# # labels = map(label->text(label, :top, 8), ena.codeModel[!, :code])
# # plot!(p, x, y,
# #     label=nothing,
# #     seriestype=:scatter,
# #     series_annotations=labels,
# #     markershape=:circle,
# #     markersize=codeWidths,
# #     markercolor=:black,
# #     markerstrokewidth=0)
    
# end

# ### Helper Placeholder - extras to add to the distribution subplot
# function plot_extras!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
# kwargs...)

# # do nothing
# end

# ### Helper - draw the confidence intervals
# function plot_cis!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1}, groupName::Any;
# color::Colorant=colorant"black",
# flipX::Bool=false, flipY::Bool=false,
# showCIs::Bool=true,
# kwargs...)

# xs, ys = help_xs_and_ys(ena, displayRows, flipX, flipY)
# help_plot_ci(p, xs, ys, color, :square, "$(groupName) Mean", showCIs)
# end
