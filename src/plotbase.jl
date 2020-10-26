## Plotting
### Top-level wrapper
function plot(ena::AbstractENAModel;
    margin=10mm, size=500, lims=1, ticks=[-1, 0, 1],
    titles=[], xlabel="", ylabel="", leg=:bottomleft,
    negColor::Colorant=colorant"red", posColor::Colorant=colorant"blue",
    extraColors=[colorant"purple", colorant"orange", colorant"green",
                 colorant"pink", colorant"cyan", colorant"tan",
                 colorant"brown", colorant"orchid", colorant"navy",
                 colorant"olive", colorant"magenta"],
    flipX=false, flipY=false,
    singleUnit=nothing, groupBy=nothing,
    kwargs...)

    #### Combine the kwargs
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
        plot(leg=leg, margin=margin, size=(size, size))  # predictive
        # plot(leg=leg, margin=margin, size=(size, size))  # errors
    ]

    #### Draw usual subplots
    allRows = [true for row in eachrow(ena.metadata)]
    plot_units!(ps[1], ena, allRows; color=colorant"black", kwargs...)
    plot_network!(ps[1], ena, allRows; color=colorant"black", kwargs...)
    plot_extras!(ps[1], ena, allRows; kwargs...)
    title!(ps[1], "(a) " * get(titles, 1, "Omnibus"))
    results = test(ena)
    xlabel!(ps[1], "$xlabel ($(round(Int, results[:variance_x]*100))%)")
    ylabel!(ps[1], "$ylabel ($(round(Int, results[:variance_y]*100))%)")

    title!(ps[2], "(b) " * get(titles, 2, "Predictive"))
    plot_predictive!(ps[2], ena; kwargs...)

    # title!(ps[3], "(c) " * get(titles, 3, "Errors"))
    # plot_errors!(ps[3], ena, allRows; kwargs...)

    #### Initialize group-wise subplots
    if !isnothing(groupBy)
        groups = sort(unique(ena.metadata[!, groupBy]))
        for (g, group) in enumerate(groups)

            #### Draw group-wise subplots
            letters = "abcdefghijklmnopqrstuvwxyz"[(1+length(ps)):end]
            if g <= length(extraColors) && g <= length(letters)
                p = plot(leg=leg, margin=margin, size=(size, size))
                groupRows = [row[groupBy] == group for row in eachrow(ena.metadata)]
                plot_units!(p, ena, groupRows; color=extraColors[g], kwargs...)
                plot_network!(p, ena, groupRows; color=extraColors[g], kwargs...)
                xs = ena.centroidModel[groupRows, :pos_x] * (flipX ? -1 : 1)
                ys = ena.centroidModel[groupRows, :pos_y] * (flipY ? -1 : 1)
                help_plot_ci(p, xs, ys, extraColors[g], :square, "$(group) Mean")
                help_plot_ci(ps[1], xs, ys, extraColors[g], :square, "$(group) Mean")
                title!(p, "($(letters[g])) " * string(get(titles, 2+g, group)))
                push!(ps, p)
            end
        end
    end

    #### Layout the subplots
    for p in ps
        xticks!(p, ticks)
        yticks!(p, ticks)
        xlims!(p, -lims, lims)
        ylims!(p, -lims, lims)
    end

    N = ceil(Int, sqrt(length(ps)))
    M = ceil(Int, length(ps)/N)
    layout = grid(N, M)
    while length(ps) < N*M
        push!(ps, plot(legend=false,grid=false,foreground_color_subplot=:white))
    end

    p = plot(ps..., size=(size*M, size*N), layout=layout)
    # title!(p, title)
    return p
end

### Helper - Draw the dots
function plot_units!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
    color::Colorant=colorant"black",
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    #### Get the x/y positions
    displayCentroids = ena.centroidModel[displayRows, :]
    x = displayCentroids[!, :pos_x] * (flipX ? -1 : 1)
    y = displayCentroids[!, :pos_y] * (flipY ? -1 : 1)

    #### Draw them in black
    plot!(p, x, y,
        label="Units",
        seriestype=:scatter,
        markershape=:circle,
        markersize=2,
        markercolor=color,
        markerstrokecolor=color)
end

### Helper - Draw the lines
function plot_network!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
    color::Colorant=colorant"black",
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    #### Find the true weight on each line
    displayAccums = ena.accumModel[displayRows, :]
    lineWidths = map(eachrow(ena.networkModel)) do networkRow
        return sum(displayAccums[!, networkRow[:relationship]])
    end

    #### Rescale the lines
    lineWidths *= 2 / maximum(lineWidths)

    #### Initialize code widths, compute while we visit each line
    codeWidths = zeros(nrow(ena.codeModel))

    #### For each line...
    for (i, networkRow) in enumerate(eachrow(ena.networkModel))
        j, k = ena.relationshipMap[networkRow[:relationship]]

        #### ...add to its code weights
        codeWidths[j] += lineWidths[i]
        codeWidths[k] += lineWidths[i]

        #### and plot that line
        x = ena.codeModel[[j, k], :pos_x] * (flipX ? -1 : 1)
        y = ena.codeModel[[j, k], :pos_y] * (flipY ? -1 : 1)
        plot!(p, x, y,
            label=nothing,
            seriestype=:line,
            linewidth=lineWidths[i],
            linecolor=color)
    end

    #### Rescale and draw the codes
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

### Helper - Draw the predictive lines
function plot_predictive!(p::Plot, ena::AbstractENAModel;
    negColor::Colorant=colorant"red", posColor::Colorant=colorant"blue",
    flipX::Bool=false, flipY::Bool=false, VOIMode::Bool=false,
    kwargs...)

    ### Grab the data we need as one data frame
    regressionData = hcat(ena.accumModel, ena.metadata, makeunique=true)
    regressionData[!, :pos_x] = ena.centroidModel[!, :pos_x] * (flipX ? -1 : 1)
    regressionData[!, :pos_y] = ena.centroidModel[!, :pos_y] * (flipY ? -1 : 1)

    ### Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
    for networkRow in eachrow(ena.networkModel)
        regressionData[!, networkRow[:relationship]] = map(Float64, regressionData[!, networkRow[:relationship]])
    end

    ### TEMP: Testing two different options
    altNetworkModel = ena.networkModel
    if VOIMode
        altNetworkModel = copy(ena.networkModel)
        rotate!(ena.rotation, altNetworkModel, ena.accumModel, ena.metadata)
    end

    ### Compute line widths as the strength (slope) between the xpos and the accum network weights
    f1 = @formula(y ~ pos_x)
    lineWidths = map(eachrow(altNetworkModel)) do networkRow
        if VOIMode
            return networkRow[:weight_x]
        else
            f1 = FormulaTerm(term(networkRow[:relationship]), f1.rhs)
            try
                ## The function call is different when we have contrasts
                m1 = fit(LinearModel, f1, regressionData)
                slope = coef(m1)[2]
                return slope
            catch e
                println(e)
                error("""
                An error occured running a regression during the rotation step of this ENA model.
                Usually, this occurs because the data, the regression model, and regression formula are not in agreement.
                If you are using a MeansRotation, then this usually means that your accidentally grouped your
                units on a different variable than the variable you passed to your MeansRotation.
                """)
            end
        end
    end

    ### Use purple for negative widths and orange for positive
    lineColors = map(lineWidths) do width
        if width < 0
            return negColor
        else
            return posColor
        end
    end

    ### Convert negatives back to positive, then rescale
    lineWidths = abs.(lineWidths)
    lineWidths *= 2 / maximum(lineWidths)

    ### Placeholder, let's compute code weights as we visit each line
    codeWidths = zeros(nrow(ena.codeModel))

    ### For each line...
    for (i, networkRow) in enumerate(eachrow(ena.networkModel))

        ### ...contribute to the code weights...
        j, k = ena.relationshipMap[networkRow[:relationship]]
        codeWidths[j] += lineWidths[i]
        codeWidths[k] += lineWidths[i]

        ### ...and plot that line, in the right width and color
        x = ena.codeModel[[j, k], :pos_x] * (flipX ? -1 : 1)
        y = ena.codeModel[[j, k], :pos_y] * (flipY ? -1 : 1)
        plot!(p, x, y,
            label=nothing,
            seriestype=:line,
            linestyle=:dash,
            linewidth=lineWidths[i],
            linecolor=lineColors[i])
    end

    ### Rescale then plot the codes
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

# ### Helper - plot the errors between centroid and accum model
# function plot_errors!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
#     flipX::Bool=false, flipY::Bool=false,
#     kwargs...)
    
#     for i in 1:nrow(ena.accumModel)
#       accumPos = ena.accumModel[i, :pos_x] > 0
#       centroidPos = ena.centroidModel[i, :pos_x] > 0
#       if accumPos != centroidPos
#         x = [ena.accumModel[i, :pos_x], ena.centroidModel[i, :pos_x]] * (flipX ? -1 : 1)
#         y = [ena.accumModel[i, :pos_y], ena.centroidModel[i, :pos_y]] * (flipY ? -1 : 1)
#         plot!(p, x, y,
#             label=nothing,
#             seriestype=:line,
#             linewidth=1,
#             linecolor=:black)
        
#         x = [ena.centroidModel[i, :pos_x]] * (flipX ? -1 : 1)
#         y = [ena.centroidModel[i, :pos_y]] * (flipY ? -1 : 1)
#         labels = [text(ena.centroidModel[i, :ENA_UNIT], :top, 4)]
#         plot!(p, x, y,
#             label=nothing,
#             seriestype=:scatter,
#             series_annotations=labels,
#             markershape=:circle,
#             markersize=2,
#             markercolor=:black,
#             markerstrokecolor=:black)
#       end
#     end
# end

### Helper Placeholder - extras to add to the omnibus subplot
function plot_extras!(p::Plot, ena::AbstractENAModel, displayRows::Array{Bool,1};
    kwargs...)

    # do nothing
end