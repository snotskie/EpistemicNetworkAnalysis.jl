abstract type AbstractCarlsFormula2Rotation <: AbstractFormulaRotation
    # fields: (inherit), regression_model2, coefindex2, f2, contrasts
    # plot accepts: (inherit)
    # test reports: (inherit), pvalue2, effect_size2
end

"""
TODO document
"""
struct CarlsFormula2Rotation{T <: RegressionModel, U <: RegressionModel} <: AbstractCarlsFormula2Rotation
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
function rotate!(rotation::AbstractCarlsFormula2Rotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)
    ## TODO check assumptions about f1 and f2

    ## Collect data we need into one data frame
    regressionData = hcat(unitModel, metadata, makeunique=true)

    ## Filter out rows with missing data
    for t in rotation.f1.rhs
        if isa(t, Term)
            col = Symbol(t)
            goodRows = completecases(regressionData[!, [col]])
            regressionData = regressionData[goodRows, :]
        end
    end

    for t in rotation.f2.rhs
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

    ## Normalize the weights for the x-axes
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    ## Compute the x-axis so we can deflate it out of the data
    xAxis = map(eachrow(unitModel)) do unitRow
        return sum(
            networkRow[:weight_x] * unitRow[networkRow[:relationship]]
            for networkRow in eachrow(networkModel)
        )
    end

    ## Create a deflated copy of the data
    A = Matrix{Float64}(regressionData[!, networkModel[!, :relationship]])
    x = Vector{Float64}(networkModel[!, :weight_x])
    deflatedData = copy(regressionData)
    deflatedData[!, networkModel[!, :relationship]] = A - A*x*transpose(x)
        
    ## Factorize the deflated data
    S = Matrix{Float64}(deflatedData[!, networkModel[!, :relationship]])
    F = svd(S)
    V = F.V

    ## Keep factors with >1 eigenvalue, short circuit if there are none, else rotate that part of the regression data and continue
    nDims = count(s > 1 for s in F.S)
    V = V[:, 1:nDims]
    regressionData[!, networkModel[1:nDims, :relationship]] = A * V

    ## Run regressions again, but for the y-axis and using the second formula and using the factored data
    visitedDimensions = 0
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

        visitedDimensions += 1
        if visitedDimensions >= nDims
            break
        end
    end

    ## Un-rotate the rotation matrix to use for the y-axis
    networkModel[!, :weight_y] = Matrix{Float64}(V) * Vector{Float64}(networkModel[1:nDims, :weight_y])

    ## Normalize the weights for the y-axis
    s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
    if s != 0
        networkModel[!, :weight_y] /= s
    end
end

# Override tests
function test(ena::AbstractENAModel{<:AbstractCarlsFormula2Rotation})

    ## Get results from parent
    results = invoke(test, Tuple{AbstractENAModel{<:AbstractFormulaRotation}}, ena)

    ## Grab data we need as a single data frame
    regressionData = hcat(ena.centroidModel, ena.metadata, makeunique=true)

    ## Construct formulas
    fyab = FormulaTerm(term(:pos_y), ena.rotation.f2.rhs)
    fya = FormulaTerm(term(:pos_y), ena.rotation.f2.rhs[1:end .!= ena.rotation.coefindex2])

    ## Placeholders
    variance_yab = 0
    variance_ya = 0
    pvalue_y = 1

    ## Run the regression models; the function call is different when we have constrasts
    if isnothing(ena.rotation.contrasts)
        myab = fit(ena.rotation.regression_model, fyab, regressionData)
        mya = fit(ena.rotation.regression_model, fya, regressionData)
        variance_yab = var(predict(myab)) / var(regressionData[!, :pos_y])
        variance_ya = var(predict(mya)) / var(regressionData[!, :pos_y])
        pvalue_y = coeftable(myab).cols[4][ena.rotation.coefindex2]
    else
        myab = fit(ena.rotation.regression_model, fyab, regressionData, contrasts=ena.rotation.contrasts)
        mya = fit(ena.rotation.regression_model, fya, regressionData, contrasts=ena.rotation.contrasts)
        variance_yab = var(predict(myab)) / var(regressionData[!, :pos_y])
        variance_ya = var(predict(mya)) / var(regressionData[!, :pos_y])
        pvalue_y = coeftable(myab).cols[4][ena.rotation.coefindex2]
    end

    ## Compute f^2 and report our results and return
    f2_y = (variance_yab - variance_ya) / (1 - variance_yab)
    results[:f2_y] = f2_y
    results[:pvalue_y] = pvalue_y
    return results
end

# Override plotting pieces
## Labels - we should also report the p-value and effect size
function plot_labels!(p::Plot, ena::AbstractENAModel{<:AbstractCarlsFormula2Rotation};
    xlabel="X", ylabel="Y",
    kwargs...)
    
    ### Run tests, then put the values into the axis labels
    results = test(ena)
    xlabel!(p, "$xlabel ($(round(Int, results[:variance_x]*100))%, p<$(ceil(results[:pvalue_x], digits=4)), f²=$(round(results[:f2_x], digits=4)))")
    ylabel!(p, "$ylabel ($(round(Int, results[:variance_y]*100))%, p<$(ceil(results[:pvalue_y], digits=4)), f²=$(round(results[:f2_y], digits=4)))")
end