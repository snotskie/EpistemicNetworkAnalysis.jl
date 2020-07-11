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
function rotate!(rotation::FormulaRotation, networkModel::DataFrame, centroidModel::DataFrame)
    ## TODO check assumptions about f1

    ## Filter missing data
    filteredUnitModel = centroidModel
    for t in rotation.f1.rhs
        if isa(t, Term)
            col = Symbol(t)
            goodRows = completecases(filteredUnitModel[!, [col]])
            filteredUnitModel = filteredUnitModel[goodRows, :]
        end
    end

    ## Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        filteredUnitModel[!, r] = map(Float64, filteredUnitModel[!, r])
    end

    ## For each relationship, find the effect of the first predictor after the intercept
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        f1 = FormulaTerm(term(r), rotation.f1.rhs)
        try
            if rotation.contrasts isa Nothing
                m1 = fit(rotation.regression_model, f1, filteredUnitModel)
                slope = coef(m1)[rotation.coefindex]
                networkRow[:weight_x] = slope
            else
                m1 = fit(rotation.regression_model, f1, filteredUnitModel, contrasts=rotation.contrasts)
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
    xAxis = Matrix{Float64}(centroidModel[!, networkModel[!, :relationship]]) *
            Matrix{Float64}(networkModel[!, [:weight_x]])
    xAxis = xAxis .- mean(xAxis)
    controlModel = DataFrame(xAxis)
    pcaModel = projection(help_deflating_svd(networkModel, centroidModel, controlModel))
    networkModel[!, :weight_y] = pcaModel[:, 1]
end

# Override plotting pieces
## Labels - we should also report the p-value and effect size
function plot_labels!(p::Plot, ena::AbstractENAModel{<:AbstractFormulaRotation};
    xlabel="X", ylabel="Y",
    kwargs...)

    total_variance = sum(var.(eachcol(centroids(ena)[!, network(ena)[!, :relationship]])))
    variance_x = var(centroids(ena)[!, :pos_x]) / total_variance
    variance_y = var(centroids(ena)[!, :pos_y]) / total_variance

    fxab = FormulaTerm(term(:pos_x), ena.rotation.f1.rhs)
    fxa = FormulaTerm(term(:pos_x), ena.rotation.f1.rhs[1:end .!= ena.rotation.coefindex])
    variance_xab = 0
    variance_xa = 0
    pvalue_x = 1
    if ena.rotation.contrasts isa Nothing
        mxab = fit(ena.rotation.regression_model, fxab, centroids(ena))
        mxa = fit(ena.rotation.regression_model, fxa, centroids(ena))
        variance_xab = var(predict(mxab)) / var(centroids(ena)[!, :pos_x])
        variance_xa = var(predict(mxa)) / var(centroids(ena)[!, :pos_x])
        pvalue_x = coeftable(mxab).cols[4][ena.rotation.coefindex]
    else
        mxab = fit(ena.rotation.regression_model, fxab, centroids(ena), contrasts=ena.rotation.contrasts)
        mxa = fit(ena.rotation.regression_model, fxa, centroids(ena), contrasts=ena.rotation.contrasts)
        variance_xab = var(predict(mxab)) / var(centroids(ena)[!, :pos_x])
        variance_xa = var(predict(mxa)) / var(centroids(ena)[!, :pos_x])
        pvalue_x = coeftable(mxab).cols[4][ena.rotation.coefindex]
    end

    f2_x = (variance_xab - variance_xa) / (1 - variance_xab)
    xlabel!(p, "$xlabel ($(round(Int, variance_x*100))%, p<$(ceil(pvalue_x, digits=4)), f²=$(round(f2_x, digits=4)))")
    ylabel!(p, "$ylabel ($(round(Int, variance_y*100))%)")
end