abstract type AbstractTopicRotation <: AbstractLinearENARotation end
struct TopicRotation <: AbstractTopicRotation
    themeName::AbstractString
    controlNodes::Array{Symbol}
    treatmentNodes::Array{Symbol}
end

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractTopicRotation, M<:AbstractLinearENAModel{R}}

    embedding = similar(model.embedding, 1)
    embedding[1, :label] = model.rotation.themeName

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
 
#  # Implement rotation
#  function rotate!(rotation::AbstractTopicRotation, networkModel::DataFrame, codeModel::DataFrame, metadata::DataFrame, subspaceModel::DataFrame)
    
#     help_one_vector(networkModel, subspaceModel)
#  end
 
# # Override plotting pieces
# ## Base - Inject a groupBy and some labels when none are given
# function plot(ena::AbstractENAModel{<:AbstractTopicRotation};
#     xlabel=nothing, ylabel=nothing,
#     kwargs...)

#     if isnothing(xlabel)
#         left = join(ena.rotation.controlNodes, "/")
#         if length(left) > 15
#             left = join(map(x->string(x)[1:min(5,end)], ena.rotation.controlNodes), "/")
#         end

#         right = join(ena.rotation.treatmentNodes, "/")
#         if length(right) > 15
#             right = join(map(x->string(x)[1:min(5,end)], ena.rotation.treatmentNodes), "/")
#         end

#         xlabel = "$(left) to $(right)"
#         # xlabel = "$( to $(join(ena.rotation.treatmentNodes, '/'))"
#     end

#     if isnothing(ylabel)
#         ylabel = "SVD"
#     end

#     return invoke(plot, Tuple{AbstractENAModel{<:AbstractLinearENARotation}}, ena;
#                 xlabel=xlabel, ylabel=ylabel, kwargs...)
# end