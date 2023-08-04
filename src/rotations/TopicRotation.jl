abstract type AbstractTopicRotation <: AbstractLinearENARotation end
struct TopicRotation <: AbstractTopicRotation
    topicName::AbstractString
    controlNodes::Array{Symbol}
    treatmentNodes::Array{Symbol}
end

"""
    TopicRotation(
        topicName::AbstractString,
        controlNodes::Array{Symbol},
        treatmentNodes::Array{Symbol}
    )

Define a rotation that places its x-axis through the mean of `controlNodes` on the left and the mean of `treatmentNodes` on the right, ie., through an *a priori* defined topic

## Example
```julia
rotation = TopicRotation(
    "Gendered Language",
    [:Women, :Love],
    [:Men, :Honor]
)
```
"""
TopicRotation

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractTopicRotation, M<:AbstractLinearENAModel{R}}

    embedding = similar(model.embedding, 1)
    embedding[1, :label] = model.rotation.topicName

    controlRows = map(model.nodes.nodeID) do nodeID
        return nodeID in model.rotation.controlNodes
    end

    treatmentRows = map(model.nodes.nodeID) do nodeID
        return nodeID in model.rotation.treatmentNodes
    end

    edgeIDs = model.edges.edgeID
    for edgeID in edgeIDs
        muControl = mean(model.nodes[controlRows, edgeID])
        muTreatment = mean(model.nodes[treatmentRows, edgeID])
        embedding[1, edgeID] = muTreatment - muControl
    end

    append!(model.embedding, embedding)

    # let parent handle the rest
    super = rotationsupertype(M, AbstractTopicRotation)
    rotate!(super, model)
end