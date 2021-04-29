struct CodeNetworkRotation <: AbstractCodeNetworkRotation
    alg::Module # any module that defines a layout(adj_matrix[, n_dim]) function that returns an array of n_dim vectors
end

# Implement rotation
function rotate!(rotation::AbstractCodeNetworkRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame, codeModel::DataFrame)

    N = max(maximum(networkModel[!, :response]), maximum(networkModel[!, :referent]))
    M = median(networkModel[!, :density])
    adj = zeros(N, N)
    for row in eachrow(networkModel)
        i, j = row[:response], row[:referent]
        s = row[:density]
        adj[i, j] = adj[j, i] = s
        # if s > M
        #     adj[i, j] = adj[j, i] = 1
        # end
        # if s > 0
        #     adj[i, j] = adj[j, i] = 1/s
        # end
    end

    positions = [[0, 0] for row in eachrow(codeModel)]
    try
        positions = rotation.alg.layout(adj, 2)
    catch e
        positions = rotation.alg.layout(adj)
    end

    for row in eachrow(networkModel)
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