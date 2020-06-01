# TODO show the raw data for just one relationship

function Base.display(ena::ENAModel) # TODO should this be print, display, or show?
    println("Units (full):")
    display(ena.unitModel)
    println()
    println("Units (position):")
    show(ena.unitModel[!, [:ENA_UNIT, :dim_x, :fit_x, :dim_y, :fit_y]], allrows=true)
    println()
    println("Codes:")
    show(ena.codeModel, allrows=true)
    println()
    println("Network:")
    show(ena.networkModel, allrows=true)
    println()
    println("Model fit (p value):")
    println(ena.pvalue)
    println()
end

function Makie.plot(ena::ENAModel;
    xaxisname::String="", yaxisname::String="",
    artist::ENAArtist=DefaultArtist(),
    showprojection::Bool=false,
    unitscale::Real=1, linescale::Real=1, codescale::Real=1,
    meanscale::Real=1, textscale::Real=1)

    # Scene
    scene = Scene()

    # Artist (to choose the colors etc.)
    artist(ena, scene) do networkColors,
                   unitColors,
                   codeColors,
                   unitShapes,
                   codeShapes,
                   networkLineWidths,
                   unitMarkerSizes,
                   codeMarkerSizes,
                   confidenceIntervals
    
        # Draw Units
        if showprojection
            scatter!(scene, ena.unitModel[!, :dim_x], ena.unitModel[!, :dim_y],
                     markersize=unitMarkerSizes.*unitscale, color=:grey)
        end

        scatter!(scene, ena.unitModel[!, :fit_x], ena.unitModel[!, :fit_y],
                 markersize=unitMarkerSizes.*unitscale, color=unitColors)
        
        # Draw Lines
        for (i, networkRow) in enumerate(eachrow(ena.networkModel))
            j, k = ena.relationshipMap[networkRow[:relationship]]
            lines!(scene, ena.codeModel[[j, k], :fit_x], ena.codeModel[[j, k], :fit_y],
                   linewidth=networkLineWidths[i]*linescale, color=networkColors[i])
        end

        # Draw Codes (dots)
        scatter!(scene, ena.codeModel[!, :fit_x], ena.codeModel[!, :fit_y],
                 markersize=codeMarkerSizes.*codescale, color=codeColors)
        
        # Draw Codes (labels)
        for codeRow in eachrow(ena.codeModel)
            text!(scene, string(" ", codeRow[:code]), position=(codeRow[:fit_x], codeRow[:fit_y]),
                  textsize=0.2*textscale, align=(:left, :center))
        end

        # BUGFIX: In Cairo, the first scatter after a text goes with the text, not where we tell it too
        # This flushes that bug with a 0-sized scatter
        scatter!(scene, [0], [0], markersize=0)

        # Draw CI Boxes
        for CI in confidenceIntervals
            (mu_x, mu_y, ci_x, ci_y, color, shape, size) = CI
            scatter!(scene, [mu_x], [mu_y], markersize=size*meanscale, marker=shape, color=color)
            lines!(scene, [ci_x[1], ci_x[2]], [ci_y[1], ci_y[1]], linewidth=1, linestyle=:dash, color=color)
            lines!(scene, [ci_x[1], ci_x[2]], [ci_y[2], ci_y[2]], linewidth=1, linestyle=:dash, color=color)
            lines!(scene, [ci_x[1], ci_x[1]], [ci_y[1], ci_y[2]], linewidth=1, linestyle=:dash, color=color)
            lines!(scene, [ci_x[2], ci_x[2]], [ci_y[1], ci_y[2]], linewidth=1, linestyle=:dash, color=color)
        end

        # Label the Axes
        axis = scene[Axis]
        axis[:names, :axisnames] = (xaxisname, yaxisname)
    end

    return scene
end