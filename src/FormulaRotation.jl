"""
TODO document
"""
struct FormulaRotation{T <: RegressionModel} <: AbstractFormulaRotation
    regression_model::Type{T}
    coefindex::Int
    f1::FormulaTerm
    contrasts::Union{Nothing,Dict}
end

# Implement rotation
function rotate!(rotation::AbstractFormulaRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)
    ## TODO check assumptions about f1

    ## Grab the data we need as one data frame
    regressionData = hcat(unitModel, metadata, makeunique=true)

    ## Filter out rows with missing data
    for t in rotation.f1.rhs
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
            ## The function call is different when we have contrasts
            if isnothing(rotation.contrasts)
                m1 = fit(rotation.regression_model, f1, regressionData)
                slope = coef(m1)[rotation.coefindex]
                networkRow[:weight_x] = slope
            else
                m1 = fit(rotation.regression_model, f1, regressionData, contrasts=rotation.contrasts)
                slope = coef(m1)[rotation.coefindex]
                networkRow[:weight_x] = slope
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

    ## Normalize the axis weights
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    ## Find the first svd dim of the data orthogonal to the x weights, use these as the y weights
    xAxis = Matrix{Float64}(unitModel[!, networkModel[!, :relationship]]) *
            Matrix{Float64}(networkModel[!, [:weight_x]])
    xAxis = xAxis .- mean(xAxis)
    controlModel = DataFrame(xAxis)
    pcaModel = projection(help_deflating_svd(networkModel, unitModel, controlModel))
    networkModel[!, :weight_y] = pcaModel[:, 1]
end

# Override tests
function test(ena::AbstractENAModel{<:AbstractFormulaRotation})

    ## Get results from parent
    results = invoke(test, Tuple{AbstractENAModel{<:AbstractENARotation}}, ena)

    ## Grab data we need as a single data frame
    regressionData = hcat(ena.centroidModel, ena.metadata, makeunique=true)

    ## Construct regression formulas
    fxab = FormulaTerm(term(:pos_x), ena.rotation.f1.rhs)
    fxa = FormulaTerm(term(:pos_x), ena.rotation.f1.rhs[1:end .!= ena.rotation.coefindex])

    ## Placeholders
    variance_xab = 0
    variance_xa = 0
    pvalue_x = 1

    ## Run regression models; the function call is different when we have constrasts
    if isnothing(ena.rotation.contrasts)
        mxab = fit(ena.rotation.regression_model, fxab, regressionData)
        mxa = fit(ena.rotation.regression_model, fxa, regressionData)
        variance_xab = var(predict(mxab)) / var(regressionData[!, :pos_x])
        variance_xa = var(predict(mxa)) / var(regressionData[!, :pos_x])
        pvalue_x = coeftable(mxab).cols[4][ena.rotation.coefindex]
    else
        mxab = fit(ena.rotation.regression_model, fxab, regressionData, contrasts=ena.rotation.contrasts)
        mxa = fit(ena.rotation.regression_model, fxa, regressionData, contrasts=ena.rotation.contrasts)
        variance_xab = var(predict(mxab)) / var(regressionData[!, :pos_x])
        variance_xa = var(predict(mxa)) / var(regressionData[!, :pos_x])
        pvalue_x = coeftable(mxab).cols[4][ena.rotation.coefindex]
    end

    ## Compute f^2 and add our values to the results and return
    f2_x = (variance_xab - variance_xa) / (1 - variance_xab)
    results[:f2_x] = f2_x
    results[:pvalue_x] = pvalue_x
    return results
end

# Override plotting pieces
## Labels - we should also report the p-value and effect size
function plot_labels!(p::Plot, ena::AbstractENAModel{<:AbstractFormulaRotation};
    xlabel="X", ylabel="Y",
    kwargs...)

    ### Run tests, then put the values into the axis labels
    results = test(ena)
    xlabel!(p, "$xlabel ($(round(Int, results[:variance_x]*100))%, p<$(ceil(results[:pvalue_x], digits=4)), f²=$(round(results[:f2_x], digits=4)))")
    ylabel!(p, "$ylabel ($(round(Int, results[:variance_y]*100))%)")
end

## Units - we can color them by the coef variable, if a simple term
function plot_units!(p::Plot, ena::AbstractENAModel{<:AbstractFormulaRotation}, displayRows::Array{Bool,1};
    flipX::Bool=false, flipY::Bool=false, minLabel::Union{Nothing,String}=nothing, maxLabel::Union{Nothing,String}=nothing,
    kwargs...)

    ### Grab filtered values
    displayCentroids = ena.centroidModel[displayRows, :]
    displayMetadata = ena.metadata[displayRows, :]

    ### Default, when we don't have a good column to use for a numeric variable, use all black
    unitColors = [:black for unitRow in eachrow(displayCentroids)]
    unitRings = [:black for unitRow in eachrow(displayCentroids)]

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
                colorMap = help_nonlinear_gradient(colorant"purple", colorant"white", colorant"orange")
                ringMap = help_nonlinear_gradient(colorant"purple", colorant"grey", colorant"orange")
                unitColors = map(eachrow(displayMetadata)) do unitRow
                    if !ismissing(unitRow[col])
                        index = round(Int, (length(colorMap) - 1) * (unitRow[col] - lo) / (hi - lo) + 1)
                        return colorMap[index]
                    else
                        return :black
                    end
                end

                unitRings = map(eachrow(displayMetadata)) do unitRow
                    if !ismissing(unitRow[col])
                        index = round(Int, (length(ringMap) - 1) * (unitRow[col] - lo) / (hi - lo) + 1)
                        return ringMap[index]
                    else
                        return :black
                    end
                end
            end
        end
    end

    ### Find the x/y positions
    x = displayCentroids[!, :pos_x] * (flipX ? -1 : 1)
    y = displayCentroids[!, :pos_y] * (flipY ? -1 : 1)

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
            markersize=4,
            markerstrokewidth=1,
            markercolor=:purple,
            markerstrokecolor=:purple)
        
        plot!(p, [-999], [-999],
            label=maxLabel,
            seriestype=:scatter,
            markershape=:circle,
            markersize=4,
            markerstrokewidth=1,
            markercolor=:orange,
            markerstrokecolor=:orange)
        
        ### ...and plot the points with the correct gradient colors
        plot!(p, x, y,
            label=nothing,
            seriestype=:scatter,
            markershape=:circle,
            markersize=4,
            markerstrokewidth=1,
            markercolor=unitColors,
            markerstrokecolor=unitRings)
    else
        ### Otherwise, just show units in black
        plot!(p, x, y,
            label="Units",
            seriestype=:scatter,
            markershape=:circle,
            markersize=2,
            markercolor=unitColors,
            markerstrokecolor=unitColors)
    end
end

# ## Extras - we can add a "litmus" strip to illustrate the strength of the continuous effect
# function plot_extras!(p::Plot, ena::AbstractENAModel{<:AbstractFormulaRotation}, displayRows::Array{Bool,1};
#     flipX::Bool=false, flipY::Bool=false,
#     kwargs...)
    
#     ### Grab filtered values
#     displayCentroids = ena.centroidModel[displayRows, :]
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

#                 ### ...then let's prepare our gradient
#                 # colorMap = range(colorant"purple", colorant"orange", length=101)
#                 med = mean(vals) #lo + (hi-lo)/2 #median(vals)
#                 medIndex = round(Int, 100 * (med - lo) / (hi - lo) + 1)
#                 # colorMap = vcat(
#                 #     range(colorant"purple", colorant"white", length=medIndex-1),
#                 #     [colorant"white"],
#                 #     range(colorant"white", colorant"orange", length=100-medIndex+1)
#                 # )
#                 colorMap = FORMULA_ROTATION_COLOR_MAP

#                 ### ...and for each bin...
#                 for i in 1:bins

#                     ### ...find the rows in that bin, on both sides of the axis
#                     left = -1 + (i-1)/bins
#                     right = -1 + i/bins
#                     rowsInNegRange = map(eachrow(displayCentroids)) do unitRow
#                         if left - .5/bins <= unitRow[:pos_x] && unitRow[:pos_x] < right + .5/bins
#                             return true
#                         else
#                             return false
#                         end
#                     end

#                     rowsInPosRange = map(eachrow(displayCentroids)) do unitRow
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
#                         # index = round(Int, 100 * (avgInNegRange - lo) / (hi - lo) + 1)
#                         index = round(Int, (length(colorMap) - 1) * (avgInNegRange - lo) / (hi - lo) + 1)
#                         color = colorMap[index]
#                         plot!(p, [left, right], [-.9, -.9],
#                             label=nothing,
#                             seriestype=:line,
#                             linewidth=8,
#                             linecolor=color)
#                     end

#                     if any(rowsInPosRange)
#                         unitsInPosRange = displayMetadata[rowsInPosRange, :]
#                         avgInPosRange = mean(unitRow[col] for unitRow in eachrow(unitsInPosRange) if !ismissing(unitRow[col]))
#                         # index = round(Int, 100 * (avgInPosRange - lo) / (hi - lo) + 1)
#                         index = round(Int, (length(colorMap) - 1) * (avgInPosRange - lo) / (hi - lo) + 1)
#                         color = colorMap[index]
#                         plot!(p, [-right, -left], [-.9, -.9],
#                             label=nothing,
#                             seriestype=:line,
#                             linewidth=8,
#                             linecolor=color)
#                     end
#                 end
#             end
#         end
#     end
# end