abstract type AbstractManualRotation <: AbstractLinearENARotation end
struct ManualRotation <: AbstractManualRotation
    embedding::DataFrame
end

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractManualRotation, M<:AbstractLinearENAModel{R}}

    embedding = similar(model.embedding, nrow(model.rotation.embedding))
    embedding.label .= model.rotation.embedding.label
    manualEdges = Symbol.(names(model.rotation.embedding))
    for edgeID in model.edges.edgeID
        if edgeID in manualEdges
            embedding[!, edgeID] .= model.rotation.embedding[!, edgeID]
        else
            @warn "Edge $(edgeID) does not exist in the manual embedding, setting loading to zero instead"
            embedding[!, edgeID] .= 0
        end
    end

    append!(model.embedding, embedding)

    # let parent handle the rest
    super = rotationsupertype(M, AbstractManualRotation)
    rotate!(super, model)
end