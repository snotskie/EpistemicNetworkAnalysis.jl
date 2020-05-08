# TODO text representation
# TODO graphical representation, default
# TODO graphical, with options
# TODO option to select what line thickness to plot, the :thickness, or the :weight_x or :weight_y, color coded for +/-

function Base.display(ena::ENAModel)
    scene = Scene()
    scatter!(scene, ena.unitModel[!, :dim_x], ena.unitModel[!, :dim_y], markersize=0.05, color=:red)
    scatter!(scene, ena.unitModel[!, :fit_x], ena.unitModel[!, :fit_y], markersize=0.05)
    for (i, r) in enumerate(keys(ena.relationships))
        j, k = ena.relationships[r]
        linewidth = ena.networkModel[i, :thickness]
        lines!(scene, ena.codeModel[[j, k], :fit_x], ena.codeModel[[j, k], :fit_y], linewidth=linewidth)
    end

    scatter!(scene, ena.codeModel[!, :fit_x], ena.codeModel[!, :fit_y], markersize=0.1*ena.codeModel[!, :thickness])
    for codeRow in eachrow(ena.codeModel)
        text!(scene, string(" ", codeRow[:code]), position=(codeRow[:fit_x], codeRow[:fit_y]), textsize=0.2, align=(:left, :center))
    end

    # TODO confidence intervals

    display(scene)
end