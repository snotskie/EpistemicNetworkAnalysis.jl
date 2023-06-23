abstract type AbstractLDARotation <: AbstractLinearENARotation end
struct LDARotation <: AbstractLDARotation
    groupVar::Symbol
end

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLDARotation, M<:AbstractLinearENAModel{R}}

    # Prepare the data
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

function test!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLDARotation, M<:AbstractLinearENAModel{R}}

    super = rotationsupertype(M, AbstractLDARotation)
    test!(super, model)

    groups = sort(unique(model.metadata[!, model.rotation.groupVar]))
    for i in 1:nrow(model.embedding)
        if model.embedding[i, :label] == "LDA$(i)"
            test!(M, model, KruskalWallisTest, dim=i, groupVar=model.rotation.groupVar, groups=groups)
        end
    end
end

# insert groups
function defaultplotkwargs(
        ::Type{M},
        model::AbstractLinearENAModel;
        groupBy::Union{Symbol,Nothing}=model.rotation.groupVar,
        kwargs...
    ) where {R<:AbstractLDARotation, M<:AbstractLinearENAModel{R}}

    kwargs = NamedTuple(kwargs)
    defaults = (
        groupBy=groupBy,
        kwargs...
    )

    super = rotationsupertype(M, AbstractLDARotation)
    return defaultplotkwargs(super, model, merge(defaults, kwargs))
end