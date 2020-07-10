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
## CIs - we can mark how well the formula predicts the axes
function plot_intervals!(p::Plot, ena::AbstractENAModel{<:AbstractFormulaRotation}, displayCentroids::DataFrame, displayCounts::DataFrame;
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    xslope, xlow, xhigh = 0, -1, 1
    fx = FormulaTerm(term(:pos_x), ena.rotation.f1.rhs)
    if ena.rotation.contrasts isa Nothing
        mx = fit(ena.rotation.regression_model, fx, centroids(ena))
        xslope = coef(mx)[ena.rotation.coefindex] / 2
        xlow = confint(mx)[ena.rotation.coefindex,1] / 2
        xhigh = confint(mx)[ena.rotation.coefindex,2] / 2
    else
        mx = fit(ena.rotation.regression_model, fx, centroids(ena), contrasts=ena.rotation.contrasts)
        xslope = coef(mx)[ena.rotation.coefindex] / 2
        xlow = confint(mx)[ena.rotation.coefindex,1] / 2
        xhigh = confint(mx)[ena.rotation.coefindex,2] / 2
    end
    
    if flipX
        xlow, xhigh = -xlow, -xhigh
    end

    x1 = -0.25 * xslope / xlow
    x2 = -0.25 * xslope / xhigh
    x3 = +0.25 * xslope / xlow
    x4 = +0.25 * xslope / xhigh
    plot!(p, [x1, x2], [0, 0],
        seriestype=:line,
        linewidth=4,
        linecolor=:purple)
    
    plot!(p, [x3, x4], [0, 0],
        seriestype=:line,
        linewidth=4,
        linecolor=:orange)
end