# TODO verify
function help_deflating_svd(networkModel::DataFrame, unitModel::DataFrame, controlModel::Union{Nothing,DataFrame}=nothing)
    # mean center X
    rawCounts = unitModel[!, networkModel[!, :relationship]]
    countMeans = collect(mean(rawCounts[!, i]) for i in 1:size(rawCounts)[2])
    X = Matrix{Float64}(DataFrame(rawCounts .- transpose(countMeans)))

    # for each column in the controls, and for each column in the unit data,
    # deflate the unit data by subtracting out of it the projection of itself
    # onto that control column
    if !isnothing(controlModel)

        # mean center Y
        controlMeans = collect(mean(controlModel[!, i]) for i in 1:size(controlModel)[2])
        Y = Matrix{Float64}(DataFrame(controlModel .- transpose(controlMeans)))
        for i in 1:size(X)[2]
            xcol = X[:, i]
            for j in 1:size(Y)[2]
                v = Y[:, j]
                scalar = dot(xcol, v) / dot(v, v)
                xcol -= scalar * v
            end

            X[:, i] = xcol
        end
    end

    # then, once we've deflated or not, we run an SVD on the data
    Z = transpose(X) * X
    pcaModel = projection(fit(PCA, Matrix{Float64}(Z), pratio=1.0))
    return pcaModel
end