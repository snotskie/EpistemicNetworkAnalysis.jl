# see https://www.pnas.org/content/113/51/14662
function help_ac_svd(networkModel::DataFrame, unitModel::DataFrame, controlModel::Union{Nothing,DataFrame}=nothing, lambda::Real=5)
    rawCounts = unitModel[!, [networkRow[:relationship] for networkRow in eachrow(networkModel)]]
    countMeans = collect(mean(rawCounts[!, i]) for i in 1:size(rawCounts)[2])
    X = Matrix{Float64}(rawCounts .- transpose(countMeans))
    Z = transpose(X) * X
    if !isnothing(controlModel)
        controlMeans = collect(mean(controlModel[!, i]) for i in 1:size(controlModel)[2])
        Y = Matrix{Float64}(controlModel .- transpose(controlMeans))
        K = Y * transpose(Y)
        Z = transpose(X) * X - lambda * transpose(X) * K * X
    end

    pcaModel = fit(PCA, Matrix{Float64}(Z), mean=0, pratio=1.0)
    return projection(pcaModel)
end