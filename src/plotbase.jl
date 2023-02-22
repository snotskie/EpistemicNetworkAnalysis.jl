## Plotting
### Common defaults and globals
const DEFAULT_NEG_COLOR = colorant"#cc423a"
const DEFAULT_POS_COLOR = colorant"#218ebf"
const DEFAULT_EXTRA_COLORS = [
    colorant"#56BD7C", colorant"#EF691B", colorant"#9d5dbb",
    colorant"#FBC848", colorant"#D0386C",
    colorant"#f18e9f", colorant"#9A9EAB", colorant"#ff8c39",
    colorant"#346B88"
]

const GLOBAL_MAX_LINE_SIZE = 8
const GLOBAL_MAX_CODE_SIZE = 8
const GLOBAL_UNIT_SIZE = 3

### Top-level wrapper
function plot(ena::AbstractENAModel;
    margin=10mm, size=600, lims=1, ticks=[-lims, 0, lims],
    titles=[], xlabel="X", ylabel="Y", leg=:topleft,
    negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
    extraColors::Array{<:Colorant,1}=DEFAULT_EXTRA_COLORS,
    flipX=false, flipY=false,
    singleUnit=nothing, groupBy=nothing,
    showExtras::Bool=false, showNetworks::Bool=true, showUnits::Bool=true, showMeans::Bool=true,
    kwargs...)

    #### Combine the kwargs to make them easier to pass without needing
    #### knowledge of what my child functions need
    kwargs = (
        margin=margin, size=size, lims=lims, ticks=ticks,
        titles=titles, xlabel=xlabel, ylabel=ylabel, leg=leg,
        negColor=negColor, posColor=posColor, extraColors=extraColors,
        flipX=flipX, flipY=flipY, singleUnit=singleUnit, groupBy=groupBy,
        kwargs...
    )

    #### Initialize usual subplots
    ps = [
        plot(leg=leg, margin=margin, size=(size, size)), # omnibus
        plot(leg=false, margin=margin, size=(size, size)),  # predictive x
        plot(leg=false, margin=margin, size=(size, size))  # predictive y
    ]
    
    #### Start usual subplots: Distribution
    title!(ps[1], "(a) " * get(titles, 1, "Distribution"))

    #### Draw usual subplots: Dynamics
    title!(ps[2], "(b) " * get(titles, 2, "Rate of Change by X"))
    plot_predictive!(ps[2], ena, :pos_x; kwargs...)
    title!(ps[3], "(c) " * get(titles, 3, "Rate of Change by Y"))
    plot_predictive!(ps[3], ena, :pos_y; kwargs...)

    #### If we need group-wise subplots...
    if !isnothing(groupBy)

        #### ...continue usual subplots: Distribution
        allRows = [true for row in eachrow(ena.metadata)]
        if showNetworks
            plot_network!(ps[1], ena, allRows; kwargs...)
        end

        #### ...then for each group...
        groups = sort(unique(ena.metadata[!, groupBy]))
        letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZθωερτψυπασδφγηςκλζχξωβνμΘΩΨΥΠΣΔΦΓΛΞΩ"[(1+length(ps)):end]
        letterIndices = collect(eachindex(letters)) # because we're possibly dealing with unicode
        for (g, group) in enumerate(groups)

            #### Initialize and draw them, until we run out of letters
            if g <= length(extraColors) && g <= length(letterIndices)
                p = plot(leg=false, margin=margin, size=(size, size))
                groupRows = [row[groupBy] == group for row in eachrow(ena.metadata)]
                if showNetworks
                    plot_network!(p, ena, groupRows; color=extraColors[g], kwargs...)
                end

                if showUnits
                    plot_units!(p, ena, groupRows; unitLabel="$(group) Units", color=extraColors[g], kwargs...)
                    plot_units!(ps[1], ena, groupRows; unitLabel="$(group) Units", color=extraColors[g], kwargs...)
                end

                if showExtras
                    plot_extras!(p, ena, groupRows; color=extraColors[g], kwargs...)
                end
                
                if showMeans
                    plot_cis!(p, ena, groupRows, group; color=extraColors[g], kwargs...)
                    plot_cis!(ps[1], ena, groupRows, group; color=extraColors[g], kwargs...)
                end

                title!(p, "($(letters[letterIndices[g]])) " * string(get(titles, 2+g, group)))
                push!(ps, p)
            end
        end

        #### ...then for each pair of groups...
        n = length(groups)+1
        for i in 1:length(groups)
            for j in (i+1):length(groups)
                if n <= length(letterIndices) && j <= length(extraColors)
                    posGroup = groups[j]
                    negGroup = groups[i]
                    p = plot(leg=false, margin=margin, size=(size, size))
                    plot_subtraction!(p, ena, groupBy, groups[i], groups[j];
                        kwargs..., negColor=extraColors[i], posColor=extraColors[j])

                    posGroupRows = [row[groupBy] == posGroup for row in eachrow(ena.metadata)]
                    negGroupRows = [row[groupBy] == negGroup for row in eachrow(ena.metadata)]
                    if showUnits
                        plot_units!(p, ena, posGroupRows; color=extraColors[j], kwargs...)
                        plot_units!(p, ena, negGroupRows; color=extraColors[i], kwargs...)
                    end

                    if showMeans
                        plot_cis!(p, ena, posGroupRows, posGroup; color=extraColors[j], kwargs...)
                        plot_cis!(p, ena, negGroupRows, negGroup; color=extraColors[i], kwargs...)
                    end

                    temp = "$(groups[i]) vs. $(groups[j])"
                    title!(p, "($(letters[letterIndices[n]])) " * string(get(titles, 2+n, temp)))
                    push!(ps, p)
                    n += 1
                end
            end
        end

        #### then finish usual subplots: Distribution
        if showExtras
            plot_extras!(ps[1], ena, allRows; kwargs...)
        end
    else
        #### Else just draw usual subplots: Distribution
        allRows = [true for row in eachrow(ena.metadata)]
        if showNetworks
            plot_network!(ps[1], ena, allRows; kwargs...)
        end

        if showUnits
            plot_units!(ps[1], ena, allRows; kwargs...)
        end

        if showExtras
            plot_extras!(ps[1], ena, allRows; kwargs...)
        end
    end

    #### Layout the subplots
    results = test(ena)
    for p in ps
        xticks!(p, ticks)
        yticks!(p, ticks)
        if lims > 0
            xlims!(p, -lims, lims)
            ylims!(p, -lims, lims)
        end
        
        if !isnan(results[:variance_x])
            xlabel!(p, "$xlabel ($(round(Int, results[:variance_x]*100))%)")
        else
            xlabel!(p, xlabel)
        end

        if !isnan(results[:variance_y])
            ylabel!(p, "$ylabel ($(round(Int, results[:variance_y]*100))%)")
        else
            ylabel!(p, ylabel)
        end
    end

    N = ceil(Int, sqrt(length(ps)))
    M = ceil(Int, length(ps)/N)
    layout = grid(N, M)
    while length(ps) < N*M
        push!(ps, plot(legend=false,grid=false,foreground_color_subplot=:white))
    end

    p = plot(ps..., size=(size*M, size*N), layout=layout)

    #### And done
    return p
end

### Helper - Draw the dots
function plot_units!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
    color::Colorant=colorant"black", spectralColorBy::Union{Symbol,Nothing}=nothing,
    flipX::Bool=false, flipY::Bool=false,
    unitLabel::String="Units",
    kwargs...)

    #### Get the x/y positions
    xs, ys = help_xs_and_ys(ena, displayRows, flipX, flipY)
    
    #### Optional: color code by a continuous value
    if !isnothing(spectralColorBy)
        if spectralColorBy in Symbol.(names(ena.accumModel))
            colVals = Vector{Float64}(ena.accumModel[displayRows, spectralColorBy])
            allColVals = colVals = Vector{Float64}(ena.accumModel[!, spectralColorBy])
            colVals = colVals .- minimum(allColVals)
            colVals /= maximum(allColVals)
            color = [HSL(colVal*240, 1, 0.5) for colVal in colVals]
        elseif spectralColorBy in Symbol.(names(ena.metadata))
            colVals = Vector{Float64}(ena.metadata[displayRows, spectralColorBy])
            allColVals = Vector{Float64}(ena.metadata[!, spectralColorBy])
            colVals = colVals .- minimum(allColVals)
            colVals /= maximum(allColVals)
            color = [HSL(colVal*240, 1, 0.5) for colVal in colVals]
        end
    end

    #### Draw the units, in black by default
    plot!(p, xs, ys,
        label=unitLabel,
        seriestype=:scatter,
        markershape=:circle,
        markersize=GLOBAL_UNIT_SIZE,
        markercolor=color,
        markerstrokewidth=0)
end

### Helper - Draw the lines
function plot_network!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
    color::Colorant=colorant"#aaa",
    flipX::Bool=false, flipY::Bool=false, showWarps::Bool=false, showCodeLabels::Bool=true,
    showArrows::Bool=false, showTrajectoryBy::Union{Symbol,Nothing}=nothing, showNetworkLines::Bool=true,
    rotateCodeLabels::Bool=false,
    kwargs...)

    #### Find the true weight on each line
    allLineWidths = map(eachrow(ena.networkModel)) do networkRow
        return sum(ena.accumModel[!, networkRow[:relationship]])
    end

    displayAccums = ena.accumModel[displayRows, :]
    lineWidths = map(eachrow(ena.networkModel)) do networkRow
        return sum(displayAccums[!, networkRow[:relationship]])
    end

    #### Rescale the lines
    lineWidths *= GLOBAL_MAX_LINE_SIZE / maximum(allLineWidths)

    #### Initialize code widths, compute while we visit each line
    codeWidths = zeros(nrow(ena.codeModel))

    #### For each line...
    for (i, networkRow) in enumerate(eachrow(ena.networkModel))
        j, k = ena.relationshipMap[networkRow[:relationship]]

        #### ...add to its code weights
        codeWidths[j] += lineWidths[i]
        codeWidths[k] += lineWidths[i]

        #### and plot that line
        if showNetworkLines
            pointA = [ena.codeModel[j, :pos_x] * (flipX ? -1 : 1), ena.codeModel[j, :pos_y] * (flipY ? -1 : 1)]
            pointB = [ena.codeModel[k, :pos_x] * (flipX ? -1 : 1), ena.codeModel[k, :pos_y] * (flipY ? -1 : 1)]
            if j == k
                pointB = [0.0, 0.0]
            end
            
            pointT = (pointA+pointB)/2
            if showWarps
                pointT = [networkRow[:weight_x] * (flipX ? -1 : 1), networkRow[:weight_y] * (flipY ? -1 : 1)]
            end

            points = hcat(pointA, pointT, pointT, pointT, pointB)
            # arrows = nothing
            # if showArrows
            #     arrows = true #lineWidths[i] #arrow(:closed, :head, lineWidths[i], lineWidths[i])
            # end

            plot!(p,
                points[1, :],
                points[2, :],
                label=nothing,
                seriestype=:curves,
                arrows=showArrows,
                linewidth=lineWidths[i],
                linecolor=color)
        end
    end

    #### Rescale and draw the codes
    codeWidths *= GLOBAL_MAX_CODE_SIZE / maximum(codeWidths)
    x = ena.codeModel[!, :pos_x] * (flipX ? -1 : 1)
    y = ena.codeModel[!, :pos_y] * (flipY ? -1 : 1)
    if showCodeLabels
        labels = map(ena.codeModel[!, :code], x, y) do (label, xi, yi)
            if rotateCodeLabels
                return text(label, :top, default(:xtickfontsize), rotation=help_font_angle(xi, yi))
            else
                return text(label, :top, default(:xtickfontsize))
            end
        end

        plot!(p, x, y,
            label=nothing,
            seriestype=:scatter,
            series_annotations=labels,
            markershape=:circle,
            markersize=codeWidths,
            markercolor=:black,
            markerstrokewidth=0)
    else
        plot!(p, x, y,
            label=nothing,
            seriestype=:scatter,
            # series_annotations=labels,
            markershape=:circle,
            markersize=codeWidths,
            markercolor=:black,
            markerstrokewidth=0)
    end

    #### Optional: illustrate a trajectory by a continuous, non-repeating value
    if !isnothing(showTrajectoryBy)
        smoothingData = innerjoin(ena.accumModel[displayRows, :], ena.metadata[displayRows, :], on=:ENA_UNIT)
        if showTrajectoryBy in Symbol.(names(smoothingData))
            smoothingData = combine(
                groupby(smoothingData, showTrajectoryBy),
                showTrajectoryBy => mean => showTrajectoryBy,
                :pos_x => mean => :pos_x,
                :pos_y => mean => :pos_y
            )
                
            if nrow(smoothingData) > 3
                smoothingData = sort(smoothingData, showTrajectoryBy)
                ts = smoothingData[!, showTrajectoryBy]
                ps = transpose(Matrix{Float64}(smoothingData[!, [:pos_x, :pos_y]]))
                bspline = ParametricSpline(ts, ps, k=3, bc="nearest")
                smooth_ts = range(ts[1], stop=ts[end], length=500)
                smooth_ps = transpose(bspline(smooth_ts))
                plot!(p, smooth_ps[:, 1] * (flipX ? -1 : 1), smooth_ps[:, 2] * (flipY ? -1 : 1), linecolor=:black, arrow=:closed)
            end
        end 
    end
end

### Helper - Draw the predictive lines
function plot_predictive!(p::Plot, ena::AbstractENAModel, targetCol::Symbol;
    negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
    flipX::Bool=false, flipY::Bool=false, weakLinks::Bool=true, showWarps::Bool=false,
    showCodeLabels::Bool=true, showArrows::Bool=false, reverseLineSort::Bool=false,
    rotateCodeLabels::Bool=false,
    kwargs...)

    ### Grab the data we need as one data frame
    regressionData = hcat(ena.accumModel, ena.metadata, makeunique=true)
    xs, ys = help_xs_and_ys(ena, !, flipX, flipY)
    regressionData[!, :pos_x] = xs
    regressionData[!, :pos_y] = ys

    ### Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
    for networkRow in eachrow(ena.networkModel)
        regressionData[!, networkRow[:relationship]] = map(Float64, regressionData[!, networkRow[:relationship]])
    end

    ### Compute line widths as the strength (slope) between the xpos and the accum network weights
    f1 = FormulaTerm(term(:y), term(targetCol))
    lineData = map(eachrow(ena.networkModel)) do networkRow
        r = networkRow[:relationship]
        f1 = FormulaTerm(term(r), f1.rhs)
        try
            m1 = fit(LinearModel, f1, regressionData)
            slope = coef(m1)[2]
            pearson = cor(regressionData[!, targetCol], regressionData[!, r])
            return (slope, pearson)
        catch e
            return (0, 0)
        end
    end

    ### Color the lines based on their correlation with the x position
    midColor = weighted_color_mean(0.5, RGB(negColor), RGB(posColor))
    midColor = weighted_color_mean(0.3, RGB(midColor), colorant"white")
    lineColorMap = help_nonlinear_gradient(weighted_color_mean(0.95, negColor, colorant"black"),
                                           midColor,
                                           weighted_color_mean(0.95, posColor, colorant"black"),
                                           curve=2.5)
    lineColors = map(lineData) do (slope, pearson)
        if isnan(pearson)
            pearson = 0
        end

        if flipX
            pearson *= -1
        end        

        index = 1 + round(Int, (length(lineColorMap) - 1) * (pearson + 1) / 2)
        return lineColorMap[index]
    end

    ### Size the lines based on their slope with the x position
    lineWidths = map(lineData) do (slope, pearson)
        return abs(slope)
    end

    ### Normalize
    lineWidths *= GLOBAL_MAX_LINE_SIZE / maximum(lineWidths)

    ### Placeholder, let's compute code weights as we visit each line
    codeWidths = zeros(nrow(ena.codeModel))
    codeVisible = Bool.(zeros(nrow(ena.codeModel)))

    ### For each line...
    networkData = hcat(ena.networkModel, DataFrame(:width => lineWidths, :color => lineColors, :pearson => last.(lineData)))
    for networkRow in sort(eachrow(networkData), by=row->row[:width], rev=reverseLineSort)

        ### ...contribute to the code weights...
        j, k = ena.relationshipMap[networkRow[:relationship]]
        codeWidths[j] += networkRow[:width]
        codeWidths[k] += networkRow[:width]

        ### ...and if that line should be shown...
        if weakLinks || abs(networkRow[:pearson]) >= 0.3
            ### ...plot it in the right width and color
            codeVisible[j] = true
            codeVisible[k] = true
            pointA = [ena.codeModel[j, :pos_x] * (flipX ? -1 : 1), ena.codeModel[j, :pos_y] * (flipY ? -1 : 1)]
            pointB = [ena.codeModel[k, :pos_x] * (flipX ? -1 : 1), ena.codeModel[k, :pos_y] * (flipY ? -1 : 1)]
            if j == k
                pointB = [0.0, 0.0]
            end
            
            pointT = (pointA+pointB)/2
            if showWarps
                pointT = [networkRow[:weight_x] * (flipX ? -1 : 1), networkRow[:weight_y] * (flipY ? -1 : 1)]
            end

            points = hcat(pointA, pointT, pointT, pointT, pointB)
            plot!(p,
                points[1, :],
                points[2, :],
                label=nothing,
                seriestype=:curves,
                arrows=showArrows,
                linewidth=networkRow[:width],
                linecolor=networkRow[:color])
        end
    end

    ### Rescale the code widths
    codeWidths *= GLOBAL_MAX_CODE_SIZE / maximum(codeWidths)

    ### And plot the codes and we're done
    x = ena.codeModel[codeVisible, :pos_x] * (flipX ? -1 : 1)
    y = ena.codeModel[codeVisible, :pos_y] * (flipY ? -1 : 1)
    if showCodeLabels
        labels = map(ena.codeModel[codeVisible, :code], x, y) do (label, xi, yi)
            if rotateCodeLabels
                return text(label, :top, default(:xtickfontsize), rotation=help_font_angle(xi, yi))
            else
                return text(label, :top, default(:xtickfontsize))
            end
        end

        plot!(p, x, y,
            label=nothing,
            seriestype=:scatter,
            series_annotations=labels,
            markershape=:circle,
            markersize=codeWidths,
            markercolor=:black,
            markerstrokewidth=0)
    else
        plot!(p, x, y,
            label=nothing,
            seriestype=:scatter,
            # series_annotations=labels,
            markershape=:circle,
            markersize=codeWidths,
            markercolor=:black,
            markerstrokewidth=0)
    end
end

### Helper - Draw the subtraction lines (nearly identical to plot_predictive)
function plot_subtraction!(p::Plot, ena::AbstractENAModel, groupVar::Symbol, negGroup::Any, posGroup::Any;
    negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
    flipX::Bool=false, flipY::Bool=false, weakLinks::Bool=true, showWarps::Bool=false,
    showCodeLabels::Bool=true, showArrows::Bool=false, reverseLineSort::Bool=false,
    rotateCodeLabels::Bool=false,
    kwargs...)

    ### Grab the data we need as one data frame
    regressionData = hcat(ena.accumModel, ena.metadata, makeunique=true)
    xs, ys = help_xs_and_ys(ena, !, flipX, flipY)
    regressionData[!, :pos_x] = xs
    regressionData[!, :pos_y] = ys

    ### Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
    for networkRow in eachrow(ena.networkModel)
        regressionData[!, networkRow[:relationship]] = map(Float64, regressionData[!, networkRow[:relationship]])
    end

    regressionData[!, :SubtractionVar] = map(eachrow(regressionData)) do row
        if row[groupVar] == posGroup
            return 1
        elseif row[groupVar] == negGroup
            return 0
        else
            return missing
        end
    end

    rowsForCor = map(eachrow(regressionData)) do row
        return !ismissing(row[:SubtractionVar])
    end

    ### Compute line widths as the strength (slope) between the xpos and the accum network weights
    f1 = @formula(y ~ SubtractionVar)
    lineData = map(eachrow(ena.networkModel)) do networkRow
        r = networkRow[:relationship]
        f1 = FormulaTerm(term(r), f1.rhs)
        try
            m1 = fit(LinearModel, f1, regressionData)
            slope = coef(m1)[2]
            pearson = cor(regressionData[rowsForCor, :SubtractionVar], regressionData[rowsForCor, r])
            return (slope, pearson)
        catch e
            return (0, 0)
        end
    end

    ### Color the lines based on their correlation with the x position
    midColor = weighted_color_mean(0.5, RGB(negColor), RGB(posColor))
    midColor = weighted_color_mean(0.3, RGB(midColor), colorant"white")
    lineColorMap = help_nonlinear_gradient(weighted_color_mean(0.95, negColor, colorant"black"),
                                           midColor,
                                           weighted_color_mean(0.95, posColor, colorant"black"),
                                           curve=1)
    
    lineColors = map(lineData) do (slope, pearson)
        if isnan(pearson)
            pearson = 0
        end

        if flipX
            pearson *= -1
        end        

        index = 1 + round(Int, (length(lineColorMap) - 1) * (pearson + 1) / 2)
        return lineColorMap[index]
    end

    ### Size the lines based on their slope with the x position
    lineWidths = map(lineData) do (slope, pearson)
        return abs(slope)
    end

    ### Normalize
    lineWidths *= GLOBAL_MAX_LINE_SIZE / maximum(lineWidths)

    ### Placeholder, let's compute code weights as we visit each line
    codeWidths = zeros(nrow(ena.codeModel))
    codeVisible = Bool.(zeros(nrow(ena.codeModel)))

    ### For each line...
    networkData = hcat(ena.networkModel, DataFrame(:width => lineWidths, :color => lineColors, :pearson => last.(lineData)))
    for networkRow in sort(eachrow(networkData), by=row->row[:width], rev=reverseLineSort)

        ### ...contribute to the code weights...
        j, k = ena.relationshipMap[networkRow[:relationship]]
        codeWidths[j] += networkRow[:width]
        codeWidths[k] += networkRow[:width]

        ### ...and if that line should be shown...
        if weakLinks || abs(networkRow[:pearson]) >= 0.3
            ### ...plot it in the right width and color
            codeVisible[j] = true
            codeVisible[k] = true
            pointA = [ena.codeModel[j, :pos_x] * (flipX ? -1 : 1), ena.codeModel[j, :pos_y] * (flipY ? -1 : 1)]
            pointB = [ena.codeModel[k, :pos_x] * (flipX ? -1 : 1), ena.codeModel[k, :pos_y] * (flipY ? -1 : 1)]
            if j == k
                pointB = [0.0, 0.0]
            end
            
            pointT = (pointA+pointB)/2
            if showWarps
                pointT = [networkRow[:weight_x] * (flipX ? -1 : 1), networkRow[:weight_y] * (flipY ? -1 : 1)]
            end

            points = hcat(pointA, pointT, pointT, pointT, pointB)
            plot!(p,
                points[1, :],
                points[2, :],
                label=nothing,
                seriestype=:curves,
                arrows=showArrows,
                linewidth=networkRow[:width],
                linecolor=networkRow[:color])
        end
    end

    ### Rescale the code widths
    codeWidths *= GLOBAL_MAX_CODE_SIZE / maximum(codeWidths)

    ### And plot the codes and we're done
    x = ena.codeModel[codeVisible, :pos_x] * (flipX ? -1 : 1)
    y = ena.codeModel[codeVisible, :pos_y] * (flipY ? -1 : 1)
    if showCodeLabels
        labels = map(ena.codeModel[codeVisible, :code], x, y) do (label, xi, yi)
            if rotateCodeLabels
                return text(label, :top, default(:xtickfontsize), rotation=help_font_angle(xi, yi))
            else
                return text(label, :top, default(:xtickfontsize))
            end
        end

        plot!(p, x, y,
            label=nothing,
            seriestype=:scatter,
            series_annotations=labels,
            markershape=:circle,
            markersize=codeWidths,
            markercolor=:black,
            markerstrokewidth=0)
    else
        plot!(p, x, y,
            label=nothing,
            seriestype=:scatter,
            # series_annotations=labels,
            markershape=:circle,
            markersize=codeWidths,
            markercolor=:black,
            markerstrokewidth=0)
    end


    # ### Placeholder, let's compute code weights as we visit each line
    # codeWidths = zeros(nrow(ena.codeModel))

    # ### For each line...
    # networkData = hcat(ena.networkModel, DataFrame(:width => lineWidths, :color => lineColors))
    # for networkRow in sort(eachrow(networkData), by=row->row[:width], rev=reverseLineSort)

    #     ### ...contribute to the code weights...
    #     j, k = ena.relationshipMap[networkRow[:relationship]]
    #     codeWidths[j] += networkRow[:width]
    #     codeWidths[k] += networkRow[:width]

    #     ### ...and plot that line, in the right width and color
    #     x = ena.codeModel[[j, k], :pos_x] * (flipX ? -1 : 1)
    #     y = ena.codeModel[[j, k], :pos_y] * (flipY ? -1 : 1)
    #     plot!(p, x, y,
    #         label=nothing,
    #         seriestype=:line,
    #         linewidth=networkRow[:width],
    #         linecolor=networkRow[:color])
    # end

    # ### Rescale the code widths
    # codeWidths *= GLOBAL_MAX_CODE_SIZE / maximum(codeWidths)

    # ### And plot the codes and we're done
    # x = ena.codeModel[!, :pos_x] * (flipX ? -1 : 1)
    # y = ena.codeModel[!, :pos_y] * (flipY ? -1 : 1)
    # labels = map(label->text(label, :top, 8), ena.codeModel[!, :code])
    # plot!(p, x, y,
    #     label=nothing,
    #     seriestype=:scatter,
    #     series_annotations=labels,
    #     markershape=:circle,
    #     markersize=codeWidths,
    #     markercolor=:black,
    #     markerstrokewidth=0)
        
end

### Helper Placeholder - extras to add to the distribution subplot
function plot_extras!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
    kwargs...)

    # do nothing
end

### Helper - draw the confidence intervals
function plot_cis!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1}, groupName::Any;
    color::Colorant=colorant"black",
    flipX::Bool=false, flipY::Bool=false,
    showCIs::Bool=true,
    kwargs...)

    xs, ys = help_xs_and_ys(ena, displayRows, flipX, flipY)
    help_plot_ci(p, xs, ys, color, :square, "$(groupName) Mean", showCIs)
end
