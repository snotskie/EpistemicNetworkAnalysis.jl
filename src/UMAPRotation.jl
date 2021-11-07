struct UMAPRotation <: AbstractUMAPRotation
    randomState::Int
end

# Simplified constructor
function UMAPRotation()
    @warn "The UMAP algorithm, which is being used in your nonlinear ENA model, is non-deterministic. The random state has been set to 42 to help reproducibility, but the onus is still on the researcher to demonstrate that their results reflect the data, not the particular random state chosen."
    
    return UMAPRotation(42)
end

# Implement rotate
function rotate!(rotation::UMAPRotation, networkModel::DataFrame, codeModel::DataFrame, metadata::DataFrame, subspaceModel::DataFrame)
    # Seed
    Random.seed!(rotation.randomState)
    
    # Prepare Data
    allCols = setdiff(Symbol.(names(subspaceModel)), [:ENA_UNIT])
    relCols = networkModel[!, :relationship]
    conCols = setdiff(allCols, relCols)
    display(allCols)
    display(relCols)
    display(conCols)
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
    weights = ones(length(allCols))
    metric = WeightedEuclidean(weights/sum(weights))

    ## Params
    n_neighbors = 35 # TODO parameterize this
    min_dist = 0.0000000001 # TODO parameterize this
    n_components = 2
    n_epochs = 300
    learning_rate = 1
    init = :spectral
    spread = 1
    set_operation_ratio = 1
    local_connectivity = 1
    repulsion_strength = 1
    neg_sample_rate = 5
    
#     ## From the UMAP source, to bypass an ambiguity bug
#     knns, dists = UMAP.knn_search(X, n_neighbors, metric)
#     graph = UMAP.fuzzy_simplicial_set(knns, dists, n_neighbors, size(X, 2), local_connectivity, set_operation_ratio)
#     embedding = UMAP.initialize_embedding(graph, n_components, Val(init))
#     embedding = UMAP.optimize_embedding(graph, embedding, embedding, n_epochs, learning_rate, min_dist, spread, repulsion_strength, neg_sample_rate, move_ref=true)

#     ## Done! Copy into the model we actually care about
#     model = UMAP_(graph, hcat(embedding...), X, knns, dists)
    
    model = UMAP_(accumX; n_neighbors=knn, min_dist=0.0000000001)
    
    subspaceModel[!, :pos_x] = model.embedding[1, :]
    subspaceModel[!, :pos_y] = model.embedding[2, :]
    
    # Prepare one-hot encodings for network "elbows"
    X = nothing
    for networkRow in eachrow(ena.networkModel)
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
    networkEmbedding = UMAP.transform(model, X; n_neighbors=1, min_dist=0.0000000001, metric=metric) # TODO param just min_dist here
    networkModel[!, :weight_x] = networkEmbedding[1, :]
    networkModel[!, :weight_y] = networkEmbedding[2, :]
    
    # Unseed
    Random.seed!(Dates.value(now()))
end