struct FormulaRotation{T <: RegressionModel} <: AbstractFormulaRotation
    regression_model::Type{T}
    coefindex::Int
    f1::FormulaTerm
    contrasts::Union{Nothing,Dict}
end

# Implement rotation
function rotate!(rotation::AbstractFormulaRotation, networkModel::DataFrame, codeModel::DataFrame, metadata::DataFrame, subspaceModel::DataFrame)

    # Check assumptions
    # if nrow(subspaceModel) != nrow(metadata)
    #     error("Cannot perform a Formula-based rotation when rotateOn=:codeModel")
    # end

    ## TODO check assumptions about f1

    ## Grab the data we need as one data frame
    regressionData = innerjoin(subspaceModel, metadata, on=:ENA_UNIT)

    ## Bugfix
    rhs = rotation.f1.rhs
    if rhs isa Term
        rhs = [rhs]
    end

    ## Filter out rows with missing data
    for t in rhs
        if isa(t, Term)
            col = Symbol(t)
            goodRows = completecases(regressionData[!, [col]])
            regressionData = regressionData[goodRows, :]
        end
    end

    ## Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        regressionData[!, r] = map(Float64, regressionData[!, r])
    end

    ## For each relationship, find the effect of the requested predictor, use those as the axis weights
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        f1 = FormulaTerm(term(r), rotation.f1.rhs)
        try
            ## The function call is different depending on type and optional contrasts dict
            if rotation.regression_model <: GeneralizedLinearModel
                if isnothing(rotation.contrasts)
                    m1 = fit(rotation.regression_model, f1, regressionData, Normal())
                    slope = coef(m1)[rotation.coefindex]
                    networkRow[:weight_x] = slope
                else
                    m1 = fit(rotation.regression_model, f1, regressionData, Normal(), contrasts=rotation.contrasts)
                    slope = coef(m1)[rotation.coefindex]
                    networkRow[:weight_x] = slope
                end
            else
                if isnothing(rotation.contrasts)
                    m1 = fit(rotation.regression_model, f1, regressionData)
                    slope = coef(m1)[rotation.coefindex]
                    networkRow[:weight_x] = slope
                else
                    m1 = fit(rotation.regression_model, f1, regressionData, contrasts=rotation.contrasts)
                    slope = coef(m1)[rotation.coefindex]
                    networkRow[:weight_x] = slope
                end
            end
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

    help_one_vector(networkModel, subspaceModel)
end

## Units - we can color them by the coef variable, if a simple term
function plot_units!(p::Plot, ena::AbstractENAModel{<:AbstractFormulaRotation}, displayRows::Array{Bool,1};
    flipX::Bool=false, flipY::Bool=false, minLabel::Union{Nothing,String}=nothing, maxLabel::Union{Nothing,String}=nothing,
    negColor::Colorant=colorant"red", posColor::Colorant=colorant"blue",
    groupBy=nothing,
    kwargs...)

    ### Punt to the parent when we have a groupBy
    if !isnothing(groupBy)
        return invoke(plot_units!, Tuple{Plot, AbstractENAModel{<:AbstractENARotation}, Array{Bool,1}},
            p, ena, displayRows; flipX=flipX, flipY=flipY, minLabel=minLabel, maxLabel=maxLabel, negColor=negColor,
            posColor=posColor, groupBy=groupBy, kwargs...)
    end

    ### Grab filtered values
    displayMetadata = ena.metadata[displayRows, :]
    displayUnits = ena.centroidModel[displayRows, :]
    if ena.rotateOn == :accumModel
        displayUnits = ena.accumModel[displayRows, :]
    end

    ### Default, when we don't have a good column to use for a numeric variable, use all black
    unitColors = [colorant"black" for unitRow in eachrow(displayUnits)]

    ### Grab the name of the potential column as a Symbol
    col = Symbol(ena.rotation.f1.rhs[ena.rotation.coefindex])

    ### Default, don't do anything fancy with the legend
    legend_col = nothing

    ### If the column exists in the metadata...
    if col in Symbol.(names(ena.metadata))

        ### ...and the first non-missing value overall is a number...
        vals = filter(x->!ismissing(x), ena.metadata[!, col])
        if first(vals) isa Number

            ### ...and there is a non-zero range of values
            lo = minimum(vals)
            hi = maximum(vals)
            if hi != lo

                ### ...then let's use that column in the legend...
                legend_col = col

                ### ...and color-code the units based on a gradient, using black for those with missing values
                midColor = weighted_color_mean(0.5, RGB(negColor), RGB(posColor))
                midColor = weighted_color_mean(0.1, RGB(midColor), colorant"white")
                colorMap = help_nonlinear_gradient(negColor, midColor, posColor)
                unitColors = map(eachrow(displayMetadata)) do unitRow
                    if !ismissing(unitRow[col])
                        index = 1 + round(Int, (length(colorMap) - 1) * (unitRow[col] - lo) / (hi - lo))
                        return colorMap[index]
                    else
                        return colorant"#111"
                    end
                end
            end
        end
    end

    ### Find the x/y positions
    x = displayUnits[!, :pos_x] * (flipX ? -1 : 1)
    y = displayUnits[!, :pos_y] * (flipY ? -1 : 1)

    ### When we used gradient color-coding...
    if !isnothing(legend_col)

        ### ...use meaningful legend labels
        if isnothing(minLabel)
            minLabel = "Min $legend_col"
        end

        if isnothing(maxLabel)
            maxLabel = "Max $legend_col"
        end

        ### ...and add those to the legend
        plot!(p, [-999], [-999],
            label=minLabel,
            seriestype=:scatter,
            markershape=:circle,
            markersize=GLOBAL_UNIT_SIZE,
            markercolor=negColor,
            markerstrokewidth=0)
        
        plot!(p, [-999], [-999],
            label=maxLabel,
            seriestype=:scatter,
            markershape=:circle,
            markersize=GLOBAL_UNIT_SIZE,
            markercolor=posColor,
            markerstrokewidth=0)
        
        ### ...and plot the points with the correct gradient colors
        plot!(p, x, y,
            label=nothing,
            seriestype=:scatter,
            markershape=:circle,
            markersize=GLOBAL_UNIT_SIZE,
            markercolor=unitColors,
            markerstrokewidth=0)
    else
        ### Otherwise, just show units in black
        plot!(p, x, y,
            label="Units",
            seriestype=:scatter,
            markershape=:circle,
            markersize=GLOBAL_UNIT_SIZE,
            markercolor=unitColors,
            markerstrokewidth=0)
    end
end

# ## Extras - we can add a "litmus" strip to illustrate the strength of the continuous effect
# function plot_extras!(p::Plot, ena::AbstractENAModel{<:AbstractFormulaRotation}, displayRows::Array{Bool,1};
#     flipX::Bool=false, flipY::Bool=false, lims::Real=1,
#     negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
#     kwargs...)
    
#     ### Grab filtered values
#     displayUnits = ena.centroidModel[displayRows, :]
#     displayMetadata = ena.metadata[displayRows, :]

#     ### Constant: how many bins for the histogram, on each side of the axis
#     bins = 8

#     ### Grab the name of the potential column as a Symbol
#     col = Symbol(ena.rotation.f1.rhs[ena.rotation.coefindex])

#     ### If the column exists in the metadata...
#     if col in Symbol.(names(ena.metadata))

#         ### ...and the first non-missing value overall is a number...
#         vals = filter(x->!ismissing(x), ena.metadata[!, col])
#         if first(vals) isa Number

#             ### ...and there is a non-zero range of values
#             lo = minimum(vals)
#             hi = maximum(vals)
#             if hi != lo

#                 ### ...then let's prepare our gradient (same gradient used on the rings in plot_units)
#                 midColor = weighted_color_mean(0.5, RGB(negColor), RGB(posColor))
#                 midColor = weighted_color_mean(0.1, RGB(midColor), colorant"white")
#                 binMap = help_nonlinear_gradient(negColor, midColor, posColor)

#                 ### ...and prepare our placeholders
#                 negBinColors = [colorant"black" for i in 1:bins] # average color in that bin
#                 posBinColors = [colorant"black" for i in 1:bins]
#                 negBinSizes = [0.0 for i in 1:bins] # number of units in that bin
#                 posBinSizes = [0.0 for i in 1:bins]
#                 negBinScalars = [0.0 for i in 1:bins] # 0 means all were min, 1 means all were max
#                 posBinScalars = [0.0 for i in 1:bins]

#                 ### ...then for each bin...
#                 for i in 1:bins

#                     ### ...find the rows in that bin, on both sides of the axis
#                     left = -1 + (i-1)/bins
#                     right = -1 + i/bins
#                     rowsInNegRange = map(eachrow(displayUnits)) do unitRow
#                         if left - .5/bins <= unitRow[:pos_x] && unitRow[:pos_x] < right + .5/bins
#                             return true
#                         else
#                             return false
#                         end
#                     end

#                     rowsInPosRange = map(eachrow(displayUnits)) do unitRow
#                         if -right - .5/bins <= unitRow[:pos_x] && unitRow[:pos_x] < -left + .5/bins
#                             return true
#                         else
#                             return false
#                         end
#                     end

#                     ### ...then draw a piece of the litmus using the average color for that bin
#                     if any(rowsInNegRange)
#                         unitsInNegRange = displayMetadata[rowsInNegRange, :]
#                         avgInNegRange = mean(unitRow[col] for unitRow in eachrow(unitsInNegRange) if !ismissing(unitRow[col]))
#                         index = round(Int, (length(binMap) - 1) * (avgInNegRange - lo) / (hi - lo) + 1)
#                         negBinColors[i] = binMap[index]
#                         negBinSizes[i] = nrow(unitsInNegRange)
#                         negBinScalars[i] = (avgInNegRange - lo) / (hi - lo)
#                     end

#                     if any(rowsInPosRange)
#                         unitsInPosRange = displayMetadata[rowsInPosRange, :]
#                         avgInPosRange = mean(unitRow[col] for unitRow in eachrow(unitsInPosRange) if !ismissing(unitRow[col]))
#                         index = round(Int, (length(binMap) - 1) * (avgInPosRange - lo) / (hi - lo) + 1)
#                         posBinColors[i] = binMap[index]
#                         posBinSizes[i] = nrow(unitsInPosRange)
#                         posBinScalars[i] = (avgInPosRange - lo) / (hi - lo)
#                     end
#                 end

#                 ### ...then rescale the bins
#                 s = maximum(vcat(negBinSizes, posBinSizes)) / 100
#                 negBinSizes /= s
#                 posBinSizes /= s

#                 ### ...then draw those bins
#                 for i in 1:bins
#                     left = -1 + (i-1)/bins
#                     right = -1 + i/bins
#                     if true # TEMP
#                         plot!(p, [left, right], [-lims, -lims],
#                             label=nothing,
#                             seriestype=:line,
#                             linewidth=negBinSizes[i],
#                             linecolor=negBinColors[i])
                        
#                         plot!(p, [-right, -left], [-lims, -lims],
#                             label=nothing,
#                             seriestype=:line,
#                             linewidth=posBinSizes[i],
#                             linecolor=posBinColors[i])
#                     else
#                         plot!(p, [left, right], [-lims, -lims],
#                             label=nothing,
#                             seriestype=:line,
#                             linewidth=negBinSizes[i],
#                             linecolor=negColor)
                        
#                         plot!(p, [-right, -left], [-lims, -lims],
#                             label=nothing,
#                             seriestype=:line,
#                             linewidth=posBinSizes[i],
#                             linecolor=negColor)

#                         plot!(p, [left, right], [-lims, -lims],
#                             label=nothing,
#                             seriestype=:line,
#                             linewidth=negBinSizes[i] * negBinScalars[i],
#                             linecolor=posColor)
                        
#                         plot!(p, [-right, -left], [-lims, -lims],
#                             label=nothing,
#                             seriestype=:line,
#                             linewidth=posBinSizes[i] * posBinScalars[i],
#                             linecolor=posColor)
#                     end
#                 end
#             end
#         end
#     end
# end