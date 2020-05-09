# TODO show the data for just one relationship

function Base.display(ena::ENAModel)
    println("Units (full):")
    display(ena.unitModel)
    println()
    println("Units (position):")
    display(ena.unitModel[!, [:ENA_UNIT, :dim_x, :fit_x, :dim_y, :fit_y]])
    println()
    println("Codes:")
    display(ena.codeModel)
    println()
    println("Network:")
    display(ena.networkModel)
    println()
end

function Makie.plot(ena::ENAModel;
    xaxisname::String="", yaxisname::String="",
    groups::Bool=false, showprojection::Bool=false,
    markerscale::Real=1, linescale::Real=1, textscale::Real=1)

    # Tidying up
    if groups && isnothing(ena.groupVar)
        error("ENA groups plots require the ENAModel to have a groupVar.")
    end

    # Pre-splitting the groups when we're in groups mode
    controlUnits = ena.unitModel
    treatmentUnits = ena.unitModel
    if groups
        controlUnits = filter(ena.unitModel) do unitRow
            if unitRow[ena.groupVar] == ena.controlGroup
                return true
            else
                return false
            end
        end

        treatmentUnits = filter(ena.unitModel) do unitRow
            if unitRow[ena.groupVar] == ena.treatmentGroup
                return true
            else
                return false
            end
        end
    end

    # Scene
    scene = Scene()

    # Units
    if showprojection
        scatter!(scene, ena.unitModel[!, :dim_x], ena.unitModel[!, :dim_y], markersize=0.05*markerscale, color=:grey)
    end

    color = :black
    if groups
        color = map(eachrow(ena.unitModel)) do unitRow
            if unitRow[ena.groupVar] == ena.controlGroup
                return :purple
            elseif unitRow[ena.groupVar] == ena.treatmentGroup
                return :orange
            else
                return :black
            end
        end
    end
    
    scatter!(scene, ena.unitModel[!, :fit_x], ena.unitModel[!, :fit_y], markersize=0.05*markerscale, color=color)

    # Connections
    ## Default size and color
    linewidth = ena.networkModel[!, :thickness]
    color = map(linewidth) do x
        return :black
    end

    ## Groups mode size and color
    if groups

        # ## Count the size of the groups
        # numControl = sum(map(eachrow(ena.unitModel)) do unitRow
        #     if unitRow[ena.groupVar] == ena.controlGroup
        #         return 1
        #     else
        #         return 0
        #     end
        # end)

        # numTreatment = sum(map(eachrow(ena.unitModel)) do unitRow
        #     if unitRow[ena.groupVar] == ena.treatmentGroup
        #         return 1
        #     else
        #         return 0
        #     end
        # end)

        # ## Find the "direction" of the size
        # linewidth = map(eachrow(ena.networkModel)) do networkRow
        #     r = networkRow[:relationship]
        #     return sum(map(eachrow(ena.unitModel)) do unitRow
        #         if unitRow[ena.groupVar] == ena.controlGroup
        #             return -unitRow[r] / numControl
        #         elseif unitRow[ena.groupVar] == ena.treatmentGroup
        #             return +unitRow[r] / numTreatment
        #         else
        #             return 0
        #         end
        #     end)
        # end

        ## Find the "direction", "strength", and "angle" for the size
        ## Thicker lines are those whose rotation weights are more towards
        ## one side of the difference of the means
        mu_x_control = mean(controlUnits[!, :fit_x])
        mu_y_control = mean(controlUnits[!, :fit_y])
        mu_x_treatment = mean(treatmentUnits[!, :fit_x])
        mu_y_treatment = mean(treatmentUnits[!, :fit_y])
        mu_x_all = (mu_x_treatment + mu_x_control) / 2
        mu_y_all = (mu_y_treatment + mu_y_control) / 2
        vt = Vector{Float64}([
            mu_x_treatment - mu_x_all,
            mu_y_treatment - mu_y_all
        ])

        norm_vt = sqrt(sum(vt .* vt))
        for (i, r) in enumerate(keys(ena.relationships))
            j, k = ena.relationships[r]
            vl = Vector{Float64}([
                ena.networkModel[i, :weight_x],
                ena.networkModel[i, :weight_y]
            ])

            linewidth[i] = sum(vl .* vt) / norm_vt
        end

        ## Map the "direction" of the size to the colors
        color = map(linewidth) do x
            if x < 0
                return :purple
            else
                return :orange
            end
        end

        ## Renormalize the sizes
        linewidth = linewidth .^ 2
        s = maximum(linewidth)
        linewidth /= s
    end

    ## Draw the lines with their size and color
    for (i, r) in enumerate(keys(ena.relationships))
        j, k = ena.relationships[r]
        lines!(scene, ena.codeModel[[j, k], :fit_x], ena.codeModel[[j, k], :fit_y], linewidth=2*linewidth[i]*linescale, color=color[i])
    end

    # Codes with labels
    scatter!(scene, ena.codeModel[!, :fit_x], ena.codeModel[!, :fit_y], markersize=0.1*ena.codeModel[!, :thickness]*markerscale)
    for codeRow in eachrow(ena.codeModel)
        text!(scene, string(" ", codeRow[:code]), position=(codeRow[:fit_x], codeRow[:fit_y]), textsize=0.2*textscale, align=(:left, :center))
    end

    # Confidence Interval
    if !groups
        ## Whole group
        mu_x = mean(ena.unitModel[!, :fit_x])
        mu_y = mean(ena.unitModel[!, :fit_y])
        ci_x = collect(confint(OneSampleTTest(ena.unitModel[!, :fit_x])))
        ci_y = collect(confint(OneSampleTTest(ena.unitModel[!, :fit_y])))
        scatter!(scene, [mu_x], [mu_y], markersize=0.05*markerscale, marker="■")
        lines!(scene, [ci_x[1], ci_x[2]], [ci_y[1], ci_y[1]], linewidth=1, linestyle=:dash)
        lines!(scene, [ci_x[1], ci_x[2]], [ci_y[2], ci_y[2]], linewidth=1, linestyle=:dash)
        lines!(scene, [ci_x[1], ci_x[1]], [ci_y[1], ci_y[2]], linewidth=1, linestyle=:dash)
        lines!(scene, [ci_x[2], ci_x[2]], [ci_y[1], ci_y[2]], linewidth=1, linestyle=:dash)
    else
        ## Control group
        mu_x = mean(controlUnits[!, :fit_x])
        mu_y = mean(controlUnits[!, :fit_y])
        ci_x = collect(confint(OneSampleTTest(controlUnits[!, :fit_x])))
        ci_y = collect(confint(OneSampleTTest(controlUnits[!, :fit_y])))
        scatter!(scene, [mu_x], [mu_y], markersize=0.05*markerscale, marker="■", color=:purple)
        lines!(scene, [ci_x[1], ci_x[2]], [ci_y[1], ci_y[1]], linewidth=1, linestyle=:dash, color=:purple)
        lines!(scene, [ci_x[1], ci_x[2]], [ci_y[2], ci_y[2]], linewidth=1, linestyle=:dash, color=:purple)
        lines!(scene, [ci_x[1], ci_x[1]], [ci_y[1], ci_y[2]], linewidth=1, linestyle=:dash, color=:purple)
        lines!(scene, [ci_x[2], ci_x[2]], [ci_y[1], ci_y[2]], linewidth=1, linestyle=:dash, color=:purple)

        ## Treatment group
        mu_x = mean(treatmentUnits[!, :fit_x])
        mu_y = mean(treatmentUnits[!, :fit_y])
        ci_x = collect(confint(OneSampleTTest(treatmentUnits[!, :fit_x])))
        ci_y = collect(confint(OneSampleTTest(treatmentUnits[!, :fit_y])))
        scatter!(scene, [mu_x], [mu_y], markersize=0.05*markerscale, marker="■", color=:orange)
        lines!(scene, [ci_x[1], ci_x[2]], [ci_y[1], ci_y[1]], linewidth=1, linestyle=:dash, color=:orange)
        lines!(scene, [ci_x[1], ci_x[2]], [ci_y[2], ci_y[2]], linewidth=1, linestyle=:dash, color=:orange)
        lines!(scene, [ci_x[1], ci_x[1]], [ci_y[1], ci_y[2]], linewidth=1, linestyle=:dash, color=:orange)
        lines!(scene, [ci_x[2], ci_x[2]], [ci_y[1], ci_y[2]], linewidth=1, linestyle=:dash, color=:orange)
    end

    # Axes
    axis = scene[Axis]
    axis[:names, :axisnames] = (xaxisname, yaxisname)
    # TODO change the grid to just a crosshair at the origin

    return scene
end