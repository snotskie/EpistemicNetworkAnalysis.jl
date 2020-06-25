# TODO show the raw data for just one relationship

function Base.display(ena::ENAModel) # TODO should this be print, display, or show?
    println("Units (refit positions):")
    show(ena.refitUnitModel[!, [:ENA_UNIT, :pos_x, :pos_y]], allrows=true)
    println()
    println("Codes:")
    show(ena.codeModel, allrows=true)
    println()
    println("Network:")
    show(ena.networkModel, allrows=true)
    println()
    # println("Model fit (p value):")
    # println(round(ena.pvalue, digits=4))
    # println()
    println("Model fit (Pearson):")
    println(round(ena.pearson, digits=4))
    println()
    println()
    println("Model fit (Variance explained on x-axis):")
    println(round(ena.variance_x, digits=4))
    println()
    println()
    println("Model fit (Variance explained on y-axis):")
    println(round(ena.variance_y, digits=4))
    println()
end

function Plots.plot(ena::ENAModel;
    xaxisname::String="X", yaxisname::String="Y",
    artist::ENAArtist=DefaultArtist(),
    showoriginalspace::Bool=false,
    showunits::Bool=true,
    showlines::Bool=true,
    showcodes::Bool=true,
    showconfidence::Bool=false)

    # Plot
    p = plot(leg=false, margin=10mm, size=(500, 500))

    # Artist (to choose the colors etc.)
    artist(ena, p) do networkColors,
                   unitColors,
                   codeColors,
                   unitShapes,
                   codeShapes,
                   networkLineWidths,
                   unitMarkerSizes,
                   codeMarkerSizes,
                   confidenceIntervals
    
        # Draw Units
        if showoriginalspace
            plot!(p, ena.unitModel[!, :pos_x], ena.unitModel[!, :pos_y],
                seriestype=:scatter,
                markershape=unitShapes,
                markersize=unitMarkerSizes,
                markercolor=:grey,
                markerstrokecolor=:grey)
        end

        if showunits
            plot!(p, ena.refitUnitModel[!, :pos_x], ena.refitUnitModel[!, :pos_y],
                seriestype=:scatter,
                markershape=unitShapes,
                markersize=unitMarkerSizes,
                markercolor=unitColors,
                markerstrokecolor=unitColors)
        end
        
        # Draw Lines
        if showlines
            for (i, networkRow) in enumerate(eachrow(ena.networkModel))
                j, k = ena.relationshipMap[networkRow[:relationship]]
                plot!(p, ena.codeModel[[j, k], :pos_x], ena.codeModel[[j, k], :pos_y],
                    seriestype=:line,
                    linewidth=networkLineWidths[i],
                    linecolor=networkColors[i])
            end
        end

        # Draw Codes (dots)
        if showcodes
            plot!(p, ena.codeModel[!, :pos_x], ena.codeModel[!, :pos_y],
                seriestype=:scatter,
                series_annotations=map(label->text(label, :top, 5), ena.codeModel[!, :code]),
                markershape=codeShapes,
                markersize=codeMarkerSizes,
                markercolor=codeColors,
                markerstrokecolor=codeColors)
        end

        # Draw CI Boxes
        if showconfidence
            for CI in confidenceIntervals
                (mu_x, mu_y, ci_x, ci_y, color, shape, size) = CI
                plot!(p, [mu_x], [mu_y], 
                    seriestype=:scatter,
                    markersize=size,
                    markershape=shape,
                    markercolor=color,
                    markerstrokecolor=color)

                plot!(p, [ci_x[1], ci_x[2]], [ci_y[1], ci_y[1]], 
                    seriestype=:line,
                    linewidth=1,
                    linestyle=:dash,
                    linecolor=color)

                plot!(p, [ci_x[1], ci_x[2]], [ci_y[2], ci_y[2]], 
                    seriestype=:line,
                    linewidth=1,
                    linestyle=:dash,
                    linecolor=color)

                plot!(p, [ci_x[1], ci_x[1]], [ci_y[1], ci_y[2]], 
                    seriestype=:line,
                    linewidth=1,
                    linestyle=:dash,
                    linecolor=color)

                plot!(p, [ci_x[2], ci_x[2]], [ci_y[1], ci_y[2]], 
                    seriestype=:line,
                    linewidth=1,
                    linestyle=:dash,
                    linecolor=color)
            end
        end

        ## Tidy the plot
        xticks!(p, [-1, 0, 1])
        yticks!(p, [-1, 0, 1])
        xlims!(p, -1.5, 1.5)
        ylims!(p, -1.5, 1.5)
        xlabel!(p, "$xaxisname ($(round(Int, ena.variance_x*100))%)")
        ylabel!(p, "$yaxisname ($(round(Int, ena.variance_y*100))%)")
    end

    return p
end