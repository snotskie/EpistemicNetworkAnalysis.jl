abstract type AbstractLDARotation <: AbstractGroupDifferenceRotation end
struct LDARotation <: AbstractLDARotation
    groupVar::Symbol
end

"""
    LDARotation(
        groupVar::Symbol
    )

Define a rotation for comparing multiple groups, via [linear discriminant analysis](https://en.wikipedia.org/wiki/Linear_discriminant_analysis)

See also: `MeansRotation` and `MulticlassRotation`

## Example

```julia
rotation = LDARotation(:Act)
```

## Statistical Tests

Models using an `LDARotation` will run the following statistical tests:

- `KruskalWallisTest` for each dimension
"""
LDARotation

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLDARotation, M<:AbstractLinearENAModel{R}}

    # Prepare the datalinear discriminant analysis
    groups = sort(unique(model.metadata[!, model.rotation.groupVar]))
    groupMap = Dict(group => i for (i, group) in enumerate(groups))
    nc = length(groups)
    X = Matrix{Float64}(transpose(Matrix{Float64}(model.accum[!, model.edges.edgeID])))
    for j in 1:size(X, 2)
        X[:, j] = X[:, j] .- mean(X[:, j])
    end

    y = map(model.metadata[!, model.rotation.groupVar]) do group
        return groupMap[group]
    end

    ## Run the LDA
    ldaModel = projection(fit(MulticlassLDA, X, y))

    ## Add to the model
    ns = size(ldaModel)[2]
    embedding = similar(model.embedding, ns)
    for i in 1:ns
        embedding[i, :label] = "LDA$(i)"
        embedding[i, model.edges.edgeID] = ldaModel[:, i]
    end

    append!(model.embedding, embedding)

    # let parent handle the rest
    super = rotationsupertype(M, AbstractLDARotation)
    rotate!(super, model)
end