struct CodeNetworkRotation <: AbstractCodeNetworkRotation
    alg::Module # any module that defines a layout(adj_matrix[, n_dim]) function that returns an array of n_dim vectors
    f::Function # a unary function for how to transform the density of a line, defaults to inv
end

# Simplified constructor
function CodeNetworkRotation(alg::Module)
    return CodeNetworkRotation(alg, inv)
end

# Implement rotation
function rotate!(rotation::AbstractCodeNetworkRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)

    N = max(maximum(networkModel[!, :response]), maximum(networkModel[!, :referent]))
    # M = median(networkModel[!, :density])
    # N = nrow(unitModel)
    adj = zeros(N, N)
    # for (i, A) in enumerate(eachrow(unitModel))
    #     for (j, B) in enumerate(eachrow(unitModel))
    #         if i == j
    #             adj[i, j] = adj[j, i] = 0
    #         elseif i < j
    #             vecA = Vector{Float64}(A[networkModel[!, :relationship]])
    #             vecB = Vector{Float64}(B[networkModel[!, :relationship]])
    #             adj[i, j] = adj[j, i] = sqrt(sum((vecA .- vecB) .^ 2))
    #         end
    #     end
    # end
    for row in eachrow(networkModel)
        i, j = row[:response], row[:referent]
        adj[i, j] = adj[j, i] = rotation.f(row[:density])
    end

    positions = [[0, 0] for i in 1:N]
    try
        positions = rotation.alg.layout(adj, 2)
    catch e
        positions = rotation.alg.layout(adj)
    end

    for row in eachrow(networkModel)
        # vecR = Vector{Float64}(unitModel[!, row[:relationship]])
        # row[:weight_x] = sum(first.(positions) .* vecR) / sum(vecR)
        # row[:weight_y] = sum(last.(positions) .* vecR) / sum(vecR)
        i, j = row[:response], row[:referent]
        row[:weight_x] = (positions[i][1] + positions[j][1]) / 2
        row[:weight_y] = (positions[i][2] + positions[j][2]) / 2
    end

    help_two_vectors(networkModel)

    # # check assumptions
    # if nrow(networkModel) != nrow(rotation.ena.networkModel)
    #     error("Cannot copy rotation from an ENA model with a different network size.")
    # end

    # # copy the rotation
    # networkModel[!, :weight_x] = copy(rotation.ena.networkModel[!, :weight_x])
    # networkModel[!, :weight_y] = copy(rotation.ena.networkModel[!, :weight_y])
end