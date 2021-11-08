struct UMAPRotation <: AbstractUMAPRotation
    randomState::Int
    numNeighbors::Int
    minDistance::Float64
    weight::Float64
end

# Simplified constructor
function UMAPRotation(; numNeighbors=35, minDistance=0.0000000001, randomState=42, weight=1)
    
    if randomState == 42
        @warn "The UMAP algorithm, which is being used in your nonlinear ENA model, is non-deterministic. The random state has been set to 42 to help reproducibility, but the onus is still on the researcher to demonstrate that their results reflect the data, not the particular random state chosen."
    end
    
    return UMAPRotation(randomState, numNeighbors, minDistance, weight)
end

# Implement rotate
function rotate!(rotation::UMAPRotation, networkModel::DataFrame, codeModel::DataFrame, metadata::DataFrame, subspaceModel::DataFrame)
    # Seed
    Random.seed!(rotation.randomState)
    
    # Prepare Data
    allCols = setdiff(Symbol.(names(subspaceModel)), [:ENA_UNIT, :pos_x, :pos_y])
    relCols = networkModel[!, :relationship]
    conCols = setdiff(allCols, relCols)
    weight = 1 / length(allCols)
    for col in conCols
        colVals = Vector{Float64}(subspaceModel[!, col])
        colVals = colVals .- minimum(colVals)
        colVals /= maximum(colVals)
        colVals *= weight
        subspaceModel[!, col] = colVals
    end
    
    # Run Model
    X = Matrix{Float64}(transpose(Matrix{Float64}(subspaceModel[!, allCols])))
    weights = [(col in relCols ? 1 : rotation.weight) for col in allCols]
    metric = WeightedEuclidean(weights/sum(weights))
    model = UMAP_(X, 2; n_neighbors=rotation.numNeighbors, min_dist=rotation.minDistance, metric=metric)
    subspaceModel[!, :pos_x] = model.embedding[1, :]
    subspaceModel[!, :pos_y] = model.embedding[2, :]
    
    # Prepare one-hot encodings for network "elbows"
    X = nothing
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        point = zeros(length(allCols))
        for (i, col) in enumerate(allCols)
            if col == r
                point[i] = 1
            end
        end
        
        if isnothing(X)
            X = point
        else
            X = hcat(X, point)
        end
    end

    # Project network components into the embedding
    weights = [(col in relCols ? 1 : 0) for col in allCols]
    metric = WeightedEuclidean(weights/sum(weights))
    networkEmbedding = UMAP.transform(model, X; n_neighbors=1, min_dist=rotation.minDistance, metric=metric)
    networkModel[!, :weight_x] = networkEmbedding[1, :]
    networkModel[!, :weight_y] = networkEmbedding[2, :]
    
    # Unseed
    Random.default_rng()
end
