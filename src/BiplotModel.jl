struct BiplotModel{T} <: AbstractBiplotModel{T}
    # Inherits:
    codes::Array{Symbol,1}
    conversations::Array{Symbol,1}
    units::Array{Symbol,1}
    rotation::T
    accumModel::DataFrame # all the unit-level data we compute
    centroidModel::DataFrame # accumModel with re-approximated relationship columns
    metadata::DataFrame
    codeModel::DataFrame # all the code-level data we compute
    networkModel::DataFrame # all the connections-level data we compute
    relationshipMap::Any

    # Adds:
    windowSize::Int
end

function BiplotModel(data::DataFrame, codes::Array{Symbol,1}, conversations::Array{Symbol,1}, units::Array{Symbol,1};
                     kwargs...)

    # Let ENA do all the work
    ena = ENAModel(
        data, codes, conversations, units;
        kwargs..., 
        relationshipFilter=(i, j, ci, cj)->(i == j),
        windowSize=1
    )

    return BiplotModel(
        ena.codes, ena.conversations, ena.units, ena.rotation,
        ena.accumModel, ena.centroidModel, ena.metadata, ena.codeModel, ena.networkModel,
        ena.relationshipMap,
        ena.windowSize
    )
end

####### Plot base overrides : nearly identical line methods to the base, except one end of each line is at the origin ######

### Helper - Draw the lines
function plot_network!(p::Plot, ena::AbstractBiplotModel, displayRows::Array{Bool,1};
    color::Colorant=colorant"#aaa",
    flipX::Bool=false, flipY::Bool=false,
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

        #### and plot that line
        x = [ena.codeModel[j, :pos_x] * (flipX ? -1 : 1), 0]
        y = [ena.codeModel[j, :pos_y] * (flipY ? -1 : 1), 0]
        plot!(p, x, y,
            label=nothing,
            seriestype=:line,
            linewidth=lineWidths[i],
            linecolor=color)
    end

    #### Rescale and draw the codes
    codeWidths *= GLOBAL_MAX_CODE_SIZE / maximum(codeWidths)
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
        markerstrokewidth=0)
end

### Helper - Draw the predictive lines
function plot_predictive!(p::Plot, ena::AbstractBiplotModel;
    negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    ### Grab the data we need as one data frame
    regressionData = hcat(ena.accumModel, ena.metadata, makeunique=true)
    regressionData[!, :pos_x] = ena.centroidModel[!, :pos_x] * (flipX ? -1 : 1)
    regressionData[!, :pos_y] = ena.centroidModel[!, :pos_y] * (flipY ? -1 : 1)

    ### Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
    for networkRow in eachrow(ena.networkModel)
        regressionData[!, networkRow[:relationship]] = map(Float64, regressionData[!, networkRow[:relationship]])
    end

    ### Compute line widths as the strength (slope) between the xpos and the accum network weights
    f1 = @formula(y ~ pos_x)
    lineData = map(eachrow(ena.networkModel)) do networkRow
        r = networkRow[:relationship]
        f1 = FormulaTerm(term(r), f1.rhs)
        try
            m1 = fit(LinearModel, f1, regressionData)
            slope = coef(m1)[2]
            pearson = cor(regressionData[!, :pos_x], regressionData[!, r])
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

    ### For each line...
    networkData = hcat(ena.networkModel, DataFrame(:width => lineWidths, :color => lineColors))
    for networkRow in sort(eachrow(networkData), by=row->row[:width])

        ### ...contribute to the code weights...
        j, k = ena.relationshipMap[networkRow[:relationship]]
        codeWidths[j] += networkRow[:width]

        ### ...and plot that line, in the right width and color
        x = [ena.codeModel[j, :pos_x] * (flipX ? -1 : 1), 0]
        y = [ena.codeModel[j, :pos_y] * (flipY ? -1 : 1), 0]
        plot!(p, x, y,
            label=nothing,
            seriestype=:line,
            linewidth=networkRow[:width],
            linecolor=networkRow[:color])
    end

    ### Rescale the code widths
    codeWidths *= GLOBAL_MAX_CODE_SIZE / maximum(codeWidths)

    ### And plot the codes and we're done
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
        markerstrokewidth=0)
end

### Helper - Draw the subtraction lines (nearly identical to plot_predictive)
function plot_subtraction!(p::Plot, ena::AbstractBiplotModel, groupVar::Symbol, negGroup::Any, posGroup::Any;
    negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    ### Grab the data we need as one data frame
    regressionData = hcat(ena.accumModel, ena.metadata, makeunique=true)
    regressionData[!, :pos_x] = ena.centroidModel[!, :pos_x] * (flipX ? -1 : 1)
    regressionData[!, :pos_y] = ena.centroidModel[!, :pos_y] * (flipY ? -1 : 1)

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

    ### For each line...
    networkData = hcat(ena.networkModel, DataFrame(:width => lineWidths, :color => lineColors))
    for networkRow in sort(eachrow(networkData), by=row->row[:width])

        ### ...contribute to the code weights...
        j, k = ena.relationshipMap[networkRow[:relationship]]
        codeWidths[j] += networkRow[:width]

        ### ...and plot that line, in the right width and color
        x = [ena.codeModel[j, :pos_x] * (flipX ? -1 : 1), 0]
        y = [ena.codeModel[j, :pos_y] * (flipY ? -1 : 1), 0]
        plot!(p, x, y,
            label=nothing,
            seriestype=:line,
            linewidth=networkRow[:width],
            linecolor=networkRow[:color])
    end

    ### Rescale the code widths
    codeWidths *= GLOBAL_MAX_CODE_SIZE / maximum(codeWidths)

    ### And plot the codes and we're done
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
        markerstrokewidth=0)
end

### Helper Placeholder - extras to add to the distribution subplot
function plot_extras!(p::Plot, ena::AbstractBiplotModel, displayRows::Array{Bool,1};
    kwargs...)

    # do nothing
end