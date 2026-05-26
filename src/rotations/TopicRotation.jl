abstract type AbstractTopicRotation <: AbstractLinearENARotation end
struct TopicRotation <: AbstractTopicRotation
    topicName::AbstractString
    controlNodes::Array{Symbol}
    treatmentNodes::Array{Symbol}
    symmetricNormalization::Bool
    function TopicRotation(
            topicName::AbstractString,
            controlNodes::Array{<:Any},
            treatmentNodes::Array{<:Any}=[]
        )

        @assert length(controlNodes) + length(treatmentNodes) > 0 "At least one coded required for TopicRotation"
        if length(treatmentNodes) == 0
            return new(
                topicName,
                Symbol[],
                convert(Array{Symbol}, Symbol.(controlNodes))
            )
        else
            return new(
                topicName,
                convert(Array{Symbol}, Symbol.(controlNodes)),
                convert(Array{Symbol}, Symbol.(treatmentNodes))
            )
        end
    end
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

    edgeIDs = model.edges.edgeID
    embedding = similar(model.embedding, 1)
    embedding[1, :label] = model.rotation.topicName

    controlRows = map(model.nodes.nodeID) do nodeID
        return nodeID in model.rotation.controlNodes
    end

    treatmentRows = map(model.nodes.nodeID) do nodeID
        return nodeID in model.rotation.treatmentNodes
    end

    for edgeID in edgeIDs
        if sum(controlRows) > 0
            muControl = mean(model.nodes[controlRows, edgeID])
            if sum(treatmentRows) > 0
                muTreatment = mean(model.nodes[treatmentRows, edgeID])
                embedding[1, edgeID] = muTreatment - muControl
            else
                embedding[1, edgeID] = -muControl
            end
        else
            if sum(treatmentRows) > 0
                muTreatment = mean(model.nodes[treatmentRows, edgeID])
                embedding[1, edgeID] = muTreatment
            else
                @error "No control or treatment codes given for TopicRotation"
                # abort, let parent handle the rest
                super = rotationsupertype(M, AbstractTopicRotation)
                rotate!(super, model)
                return
            end
        end
    end

    append!(model.embedding, embedding)

    offTopicNodes = setdiff(Symbol.(model.nodes.nodeID), Symbol.(model.rotation.controlNodes), Symbol.(model.rotation.treatmentNodes))
    if length(offTopicNodes) > 0
        offTopicRows = map(model.nodes.nodeID) do nodeID
            return nodeID in offTopicNodes
        end

        embedding = similar(model.embedding, 1)
        if length(offTopicNodes) == 1
            embedding[1, :label] = string(offTopicNodes[1])
            for edgeID in edgeIDs
                muOffTopic = mean(model.nodes[offTopicRows, edgeID])
                embedding[1, edgeID] = muOffTopic
            end
        else
            # make a copy of node equivalent to accum and mean center
            X = Matrix{Float64}(model.nodes[offTopicRows, edgeIDs])
            for i in axes(X, 2)
                X[:, i] .-= mean(X[:, i])
            end

            # run SVD on the nodes data, then add to the model
            pca = fit(PCA, X', pratio=1.0, method=:svd) # BUGFIX force use svd, https://github.com/snotskie/EpistemicNetworkAnalysis.jl/issues/56#issuecomment-1910540698
            @debug "eigvals(pca) = $(eigvals(pca))"
            svd = transpose(projection(pca))
            @debug svd
            # BUGFIX https://github.com/snotskie/EpistemicNetworkAnalysis.jl/issues/56#issuecomment-1910540698
            # Prevent SVD from adding more dimensions than are possible
            numSVDDims = min(size(svd, 1), length(edgeIDs))
            @debug numSVDDims
            embedding = similar(model.embedding, numSVDDims)
            embedding[1:numSVDDims, edgeIDs] = svd[1:numSVDDims, :]
            embedding[1:numSVDDims, :eigen_value] = eigvals(pca)[1:numSVDDims]
            embedding[1:numSVDDims, :label] = ["OffTopic$(i)" for i in 1:numSVDDims]
        end

        append!(model.embedding, embedding)
    end

    # let parent handle the rest
    super = rotationsupertype(M, AbstractTopicRotation)
    rotate!(super, model)
end

# override default edge filter
function defaultedgefilter(
        ::Type{M},
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1},
        rotation::AbstractTopicRotation,
        config::NamedTuple
    ) where {R<:AbstractTopicRotation, M<:Union{
        AbstractPlainENAModel{R},
        AbstractDigraphENAModel{R},
        AbstractBiplotENAModel{R},
        AbstractCodewiseENAModel{R}
    }}
    super = rotationsupertype(M, AbstractTopicRotation)
    edgeFilter = defaultedgefilter(super, data, codes, conversations, units, rotation, config)
    return (row)->(
        edgeFilter(row) && (
            Symbol(row[:ground]  ) in Symbol.(rotation.controlNodes  ) ||
            Symbol(row[:ground]  ) in Symbol.(rotation.treatmentNodes) ||
            Symbol(row[:response]) in Symbol.(rotation.controlNodes  ) ||
            Symbol(row[:response]) in Symbol.(rotation.treatmentNodes)
        )
    )
end

# flip x and y, sometimes
# function defaultplotkwargs(
#         ::Type{M},
#         model::AbstractLinearENAModel;
#         x::Int=findfirst(model.embedding.label .== model.rotation.topicName),
#         y::Int=findfirst(model.embedding.label .== model.rotation.topicName) > 1 ? 1 : 2,
#         kwargs...
#     ) where {R<:AbstractTopicRotation, M<:AbstractLinearENAModel{R}}
#     kwargs = NamedTuple(kwargs)
#     defaults = (
#         x=x,
#         y=y,
#         kwargs...
#     )

#     super = rotationsupertype(M, AbstractMeansRotation)
#     return defaultplotkwargs(super, model, merge(defaults, kwargs))
# end