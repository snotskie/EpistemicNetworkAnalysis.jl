abstract type AbstractMulticlassRotation <: AbstractGroupDifferenceRotation end
struct MulticlassRotation <: AbstractMulticlassRotation
    groupVar::Symbol
end

"""
    MulticlassRotation(
        groupVar::Symbol
    )

Define a rotation for comparing multiple groups, by maximizing between-group variance

See also: `MeansRotation` and `LDARotation`

## Example

```julia
rotation = MulticlassRotation(:Act)
```

## Statistical Tests

Models using an `MulticlassRotation` will run the following statistical tests:

- `KruskalWallisTest` for each dimension
"""
MulticlassRotation

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractMulticlassRotation, M<:AbstractLinearENAModel{R}}

    edgeIDs = model.edges.edgeID
    groupNames = sort(unique(model.metadata[!, model.rotation.groupVar]))

    ## Mean Vectors
    all_rows = map(eachrow(model.accum)) do row
        return false
    end

    group_ns = []
    group_vecs = []
    for g in groupNames
        group_rows = model.metadata[!, model.rotation.groupVar] .== g
        all_rows = all_rows .| group_rows
        group_vec = Vector{Float64}(mean.(model.accum[group_rows, col] for col in edgeIDs))
        push!(group_vecs, group_vec)
        push!(group_ns, sum(group_rows))
    end

    mean_vec = Vector{Float64}(mean.(model.accum[all_rows, col] for col in edgeIDs))

    ## Offset Vectors
    offset_vecs = []
    for vec in group_vecs
        push!(offset_vecs, vec - mean_vec)
    end

    # Run Computations
    ## Between-Group Scatter Matrix
    Sb = sum(
        n * vec * transpose(vec)
        for (n, vec) in zip(group_ns, offset_vecs) 
    )

    ## Solve Eigenvalue Problem
    vals = eigvals(Sb)
    vecs = eigvecs(Sb)

    ## Add to the model
    ns = size(vecs, 2)
    embedding = similar(model.embedding, ns)
    for i in 1:ns
        embedding[i, :label] = "MCMR$(i)"
        embedding[i, :eigen_value] = vals[i]
        # vecs are stored column-major, from least to most discrimination
        axis = real.(vecs[:, end-i+1])
        axis /= sqrt(sum(axis .^ 2))
        embedding[i, model.edges.edgeID] = axis
    end

    append!(model.embedding, embedding)

    # let parent handle the rest
    super = rotationsupertype(M, AbstractMulticlassRotation)
    rotate!(super, model)
end