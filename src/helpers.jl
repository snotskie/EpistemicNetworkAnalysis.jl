function help_deflating_svd(networkModel::DataFrame, unitModel::DataFrame, controlModel::Union{Nothing,DataFrame}=nothing)
    X = Matrix{Float64}(unitModel[!, networkModel[!, :relationship]])
    for i in 1:size(X)[2]
        xcol = X[:, i]
        xcol = xcol .- mean(xcol) # mean center
        X[:, i] = xcol
    end

    if !isnothing(controlModel)
        C = Matrix{Float64}(controlModel)
        for i in 1:size(C)[2]
            ccol = C[:, i]
            ccol = ccol .- mean(ccol) # mean center
            C[:, i] = ccol
            for j in 1:size(X)[2] # deflate
                xcol = X[:, j]
                scalar = dot(xcol, ccol) / dot(ccol, ccol)
                xcol -= scalar * ccol
                X[:, j] = xcol
            end
        end
    end

    # then, once we've deflated or not, we run an SVD on the data
    pcaModel = projection(fit(PCA, X', pratio=1.0))
    return pcaModel
end