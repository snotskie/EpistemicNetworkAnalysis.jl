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
function rotate!(rotation::AbstractFormula2Rotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame, codeModel::DataFrame)

    # Check assumptions
    if nrow(unitModel) != nrow(metadata)
        error("Cannot perform a Formula-based rotation when rotateOn=:codeModel")
    end

    ## TODO check assumptions about f1 and f2

    ## Collect data we need into one data frame
    regressionData = hcat(unitModel, metadata, makeunique=true)

    ## Bugfix
    rhs1 = rotation.f1.rhs
    if rhs1 isa Term
        rhs1 = [rhs1]
    end

    rhs2 = rotation.f2.rhs
    if rhs2 isa Term
        rhs2 = [rhs2]
    end

    ## Filter out rows with missing data
    for t in rhs1
        if isa(t, Term)
            col = Symbol(t)
            goodRows = completecases(regressionData[!, [col]])
            regressionData = regressionData[goodRows, :]
        end
    end

    for t in rhs2
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

    ## For each relationship, find the effect of the chosen predictor, use that as the axis weights
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
            println(e)
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
                m2 = fit(rotation.regression_model2, f2, regressionData)
                slope = coef(m2)[rotation.coefindex2]
                networkRow[:weight_y] = slope
            else
                m2 = fit(rotation.regression_model2, f2, regressionData, contrasts=rotation.contrasts2)
                slope = coef(m2)[rotation.coefindex2]
                networkRow[:weight_y] = slope
            end
        catch e
            println(e)
            error("""
            An error occured running a regression during the rotation step of this ENA model.
            Usually, this occurs because the data, the regression model, and regression formula are not in agreement.
            """)
        end
    end

    help_two_vectors(networkModel)
end