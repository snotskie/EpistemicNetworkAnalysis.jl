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
## CIs - we can mark how well the formula predicts the axes
function plot_intervals!(p::Plot, ena::AbstractENAModel{<:AbstractFormula2Rotation}, displayCentroids::DataFrame, displayCounts::DataFrame;
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

    yslope, ylow, yhigh = 0, -1, 1
    fy = FormulaTerm(term(:pos_y), ena.rotation.f2.rhs)
    if ena.rotation.contrasts isa Nothing
        my = fit(ena.rotation.regression_model2, fy, centroids(ena))
        yslope = coef(my)[ena.rotation.coefindex2] / 2
        ylow = confint(my)[ena.rotation.coefindex2,1] / 2
        yhigh = confint(my)[ena.rotation.coefindex2,2] / 2
    else
        my = fit(ena.rotation.regression_model2, fy, centroids(ena), contrasts=ena.rotation.contrasts2)
        yslope = coef(my)[ena.rotation.coefindex2] / 2
        ylow = confint(my)[ena.rotation.coefindex2,1] / 2
        yhigh = confint(my)[ena.rotation.coefindex2,2] / 2
    end

    if flipY
        ylow, yhigh = -ylow, -yhigh
    end

    y1 = -0.25 * yslope / ylow
    y2 = -0.25 * yslope / yhigh
    y3 = +0.25 * yslope / ylow
    y4 = +0.25 * yslope / yhigh
    plot!(p, [0, 0], [y1, y2],
        seriestype=:line,
        linewidth=4,
        linecolor=:green)
    
    plot!(p, [0, 0], [y3, y4],
        seriestype=:line,
        linewidth=4,
        linecolor=:blue)
end