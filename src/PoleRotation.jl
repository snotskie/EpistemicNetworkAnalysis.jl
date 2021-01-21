"""
TODO document
"""
struct PoleRotation <: AbstractPoleRotation
    regression_model::Type{LinearModel}
    codes::Array{Symbol,1}
    coefindex1::Int
    coefindex2::Int
end

# Simplified constructor
function PoleRotation(codes::Array{Symbol,1}, left::Int, right::Int)
    return PoleRotation(LinearModel, codes, left, right)
end

# Implement rotate
function rotate!(rotation::AbstractPoleRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)
    ## TODO check assumptions about f1

    ## Grab the data we need as one data frame
    regressionData = hcat(unitModel, metadata, makeunique=true)

    ## Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        regressionData[!, r] = map(Float64, regressionData[!, r])
    end

    ## For each relationship, find TODO
    f1 = FormulaTerm(term(:col), term(1) + term(rotation.codes[rotation.coefindex1]) + term(rotation.codes[rotation.coefindex2]))
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        f1 = FormulaTerm(term(r), f1.rhs)
        try
            m1 = fit(rotation.regression_model, f1, regressionData)
            slope = coef(m1)[2] - coef(m1)[3]
            networkRow[:weight_x] = slope
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