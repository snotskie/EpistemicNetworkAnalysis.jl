# TODO graphical representation, default
# TODO graphical, with options
# TODO option to select what line thickness to plot, the :thickness, or the :weight_x or :weight_y, color coded for +/-

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

function Makie.plot(ena::ENAModel)
    scene = Scene()

    # TODO axes

    # Units
    scatter!(scene, ena.unitModel[!, :dim_x], ena.unitModel[!, :dim_y], markersize=0.05, color=:red)
    scatter!(scene, ena.unitModel[!, :fit_x], ena.unitModel[!, :fit_y], markersize=0.05, color=:grey)

    # Connections
    for (i, r) in enumerate(keys(ena.relationships))
        j, k = ena.relationships[r]
        linewidth = ena.networkModel[i, :thickness]
        color = :black
        # linewidth = abs.(ena.networkModel[i, :weight_x])
        # color = map(ena.networkModel[i, :weight_x]) do x
        #     if x < 0
        #         return :blue
        #     else
        #         return :red
        #     end
        # end
        
        lines!(scene, ena.codeModel[[j, k], :fit_x], ena.codeModel[[j, k], :fit_y], linewidth=2*linewidth, color=color)
    end

    # Codes with labels
    scatter!(scene, ena.codeModel[!, :fit_x], ena.codeModel[!, :fit_y], markersize=0.1*ena.codeModel[!, :thickness])
    for codeRow in eachrow(ena.codeModel)
        text!(scene, string(" ", codeRow[:code]), position=(codeRow[:fit_x], codeRow[:fit_y]), textsize=0.2, align=(:left, :center))
    end

    # Confidence Interval
    mu_x = mean(ena.unitModel[!, :fit_x])
    mu_y = mean(ena.unitModel[!, :fit_y])
    ci_x = collect(confint(OneSampleTTest(ena.unitModel[!, :fit_x])))
    ci_y = collect(confint(OneSampleTTest(ena.unitModel[!, :fit_y])))
    scatter!(scene, [mu_x], [mu_y], markersize=0.05, marker="■")
    lines!(scene, [ci_x[1], ci_x[2]], [ci_y[1], ci_y[1]], linewidth=1, linestyle=:dash)
    lines!(scene, [ci_x[1], ci_x[2]], [ci_y[2], ci_y[2]], linewidth=1, linestyle=:dash)
    lines!(scene, [ci_x[1], ci_x[1]], [ci_y[1], ci_y[2]], linewidth=1, linestyle=:dash)
    lines!(scene, [ci_x[2], ci_x[2]], [ci_y[1], ci_y[2]], linewidth=1, linestyle=:dash)

    return scene
end