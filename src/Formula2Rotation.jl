"""
TODO document
"""
struct Formula2Rotation{T <: RegressionModel, U <: RegressionModel} <: AbstractFormula2Rotation
    regression_model::Type{T}
    coefindex::Int
    f1::FormulaTerm
    contrasts::Union{Nothing,Dict}
    regression_model2::Type{U}
    coefindex2::Int
    f2::FormulaTerm
    contrasts2::Union{Nothing,Dict}
end

# Implement rotation
function rotate!(rotation::Formula2Rotation, networkModel::DataFrame, centroidModel::DataFrame)
    ## TODO check assumptions about f1 and f2

    ## Filter missing data
    filteredUnitModel = centroidModel
    for t in rotation.f1.rhs
        if isa(t, Term)
            col = Symbol(t)
            goodRows = completecases(filteredUnitModel[!, [col]])
            filteredUnitModel = filteredUnitModel[goodRows, :]
        end
    end

    for t in rotation.f2.rhs
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
            """)
        end
    end

    ## Again, but for the y-axis and using the second formula
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        f2 = FormulaTerm(term(r), rotation.f2.rhs)
        try
            if rotation.contrasts2 isa Nothing
                m2 = fit(rotation.regression_model2, f2, filteredUnitModel)
                slope = coef(m2)[rotation.coefindex2]
                networkRow[:weight_y] = slope
            else
                m2 = fit(rotation.regression_model2, f2, filteredUnitModel, contrasts=rotation.contrasts2)
                slope = coef(m2)[rotation.coefindex2]
                networkRow[:weight_y] = slope
            end
        catch e
            error("""
            An error occured running a regression during the rotation step of this ENA model.
            Usually, this occurs because the data, the regression model, and regression formula are not in agreement.
            """)
        end
    end

    ## Orthogonalization: replace y weights with their rejection from the x weights
    before = copy(networkModel[!, :weight_y])
    scalar = dot(networkModel[!, :weight_y], networkModel[!, :weight_x]) / dot(networkModel[!, :weight_x], networkModel[!, :weight_x])
    networkModel[!, :weight_y] -= scalar * networkModel[!, :weight_x]    
    after = copy(networkModel[!, :weight_y])

    ## Raise a warning about interpreting the y-axis when before and after have a large angle between them
    theta = dot(before, after)
    theta /= sqrt(dot(before, before))
    theta /= sqrt(dot(after, after))
    angle = acos(theta) * 180 / pi
    if abs(angle) > 5
        @warn """The angle between the y-axis and the direction of the requested effect is larger than 5 degrees ($angle degrees).
This can undermine interpreting the y-axis in terms of the requested effect."""
    end

    ## Normalize the weights
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
    if s != 0
        networkModel[!, :weight_y] /= s
    end
end

# Override plotting pieces
## Labels - we should also report the p-value and effect size
function plot_labels!(p::Plot, ena::AbstractENAModel{<:AbstractFormula2Rotation};
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

    fyab = FormulaTerm(term(:pos_y), ena.rotation.f2.rhs)
    fya = FormulaTerm(term(:pos_y), ena.rotation.f2.rhs[1:end .!= ena.rotation.coefindex2])
    variance_yab = 0
    variance_ya = 0
    pvalue_y = 1
    if ena.rotation.contrasts isa Nothing
        myab = fit(ena.rotation.regression_model, fyab, centroids(ena))
        mya = fit(ena.rotation.regression_model, fya, centroids(ena))
        variance_yab = var(predict(myab)) / var(centroids(ena)[!, :pos_y])
        variance_ya = var(predict(mya)) / var(centroids(ena)[!, :pos_y])
        pvalue_y = coeftable(myab).cols[4][ena.rotation.coefindex2]
    else
        myab = fit(ena.rotation.regression_model, fyab, centroids(ena), contrasts=ena.rotation.contrasts)
        mya = fit(ena.rotation.regression_model, fya, centroids(ena), contrasts=ena.rotation.contrasts)
        variance_yab = var(predict(myab)) / var(centroids(ena)[!, :pos_y])
        variance_ya = var(predict(mya)) / var(centroids(ena)[!, :pos_y])
        pvalue_y = coeftable(myab).cols[4][ena.rotation.coefindex2]
    end

    f2_x = (variance_xab - variance_xa) / (1 - variance_xab)
    f2_y = (variance_yab - variance_ya) / (1 - variance_yab)
    xlabel!(p, "$xlabel ($(round(Int, variance_x*100))%, p<$(ceil(pvalue_x, digits=4)), f²=$(round(f2_x, digits=4)))")
    ylabel!(p, "$ylabel ($(round(Int, variance_y*100))%, p<$(ceil(pvalue_y, digits=4)), f²=$(round(f2_y, digits=4)))")
end