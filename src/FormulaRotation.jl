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

    ## Filter missing data
    regressionData = hcat(unitModel, metadata, makeunique=true)
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

    ## For each relationship, find the effect of the first predictor after the intercept
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        f1 = FormulaTerm(term(r), rotation.f1.rhs)
        try
            if rotation.contrasts isa Nothing
                m1 = fit(rotation.regression_model, f1, regressionData)
                slope = coef(m1)[rotation.coefindex]
                networkRow[:weight_x] = slope
            else
                m1 = fit(rotation.regression_model, f1, regressionData, contrasts=rotation.contrasts)
                slope = coef(m1)[rotation.coefindex]
                networkRow[:weight_x] = slope
            end
        catch e
            error("""
            An error occured running a regression during the rotation step of this ENA model.
            Usually, this occurs because the data, the regression model, and regression formula are not in agreement.
            If you are using a MeansRotation, then this usually means that your accidentally grouped your
            units on a different variable than the variable you passed to your MeansRotation.
            """)
        end
    end

    ## Normalize the weights
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
    results = invoke(test, Tuple{AbstractENAModel{<:AbstractENARotation}}, ena)
    regressionData = hcat(ena.centroidModel, ena.metadata, makeunique=true)
    fxab = FormulaTerm(term(:pos_x), ena.rotation.f1.rhs)
    fxa = FormulaTerm(term(:pos_x), ena.rotation.f1.rhs[1:end .!= ena.rotation.coefindex])
    variance_xab = 0
    variance_xa = 0
    pvalue_x = 1
    if ena.rotation.contrasts isa Nothing
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

    results = test(ena)
    xlabel!(p, "$xlabel ($(round(Int, results[:variance_x]*100))%, p<$(ceil(results[:pvalue_x], digits=4)), f²=$(round(results[:f2_x], digits=4)))")
    ylabel!(p, "$ylabel ($(round(Int, results[:variance_y]*100))%)")
end

## Units - we can color them by the coef variable, if a simple term
function plot_units!(p::Plot, ena::AbstractENAModel{<:AbstractFormulaRotation}, displayRows::Array{Bool,1};
    flipX::Bool=false, flipY::Bool=false, minLabel::Union{Nothing,String}=nothing, maxLabel::Union{Nothing,String}=nothing,
    kwargs...)

    displayCentroids = ena.centroidModel[displayRows, :]
    displayMetadata = ena.metadata[displayRows, :]
    unitColors = [:black for unitRow in eachrow(displayCentroids)]
    col = Symbol(ena.rotation.f1.rhs[ena.rotation.coefindex])
    legend_col = nothing
    if col in Symbol.(names(ena.metadata))
        vals = filter(x->!ismissing(x), ena.metadata[!, col])
        if first(vals) isa Number
            legend_col = col
            colorMap = range(colorant"purple", colorant"orange", length=101)
            lo = minimum(vals)
            hi = maximum(vals)
            if hi != lo
                unitColors = map(eachrow(displayMetadata)) do unitRow
                    if !ismissing(unitRow[col])
                        index = round(Int, 100 * (unitRow[col] - lo) / (hi - lo) + 1)
                        return colorMap[index]
                    else
                        return :black
                    end
                end
            end
        end
    end

    x = displayCentroids[!, :pos_x] * (flipX ? -1 : 1)
    y = displayCentroids[!, :pos_y] * (flipY ? -1 : 1)
    if !isnothing(legend_col)
        if isnothing(minLabel)
            minLabel = "Min $legend_col"
        end

        if isnothing(maxLabel)
            maxLabel = "Max $legend_col"
        end

        plot!(p, x, y,
            label=minLabel,
            seriestype=:scatter,
            markershape=:circle,
            markersize=2,
            markercolor=:purple,
            markerstrokecolor=:purple)
        
        plot!(p, x, y,
            label=maxLabel,
            seriestype=:scatter,
            markershape=:circle,
            markersize=2,
            markercolor=:orange,
            markerstrokecolor=:orange)
        
        plot!(p, x, y,
            label=nothing,
            seriestype=:scatter,
            markershape=:circle,
            markersize=2,
            markercolor=unitColors,
            markerstrokecolor=unitColors)
    else
        plot!(p, x, y,
            label="Units",
            seriestype=:scatter,
            markershape=:circle,
            markersize=2,
            markercolor=unitColors,
            markerstrokecolor=unitColors)
    end
end