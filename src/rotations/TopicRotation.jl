abstract type AbstractTopicRotation <: AbstractLinearENARotation end
struct TopicRotation <: AbstractTopicRotation
    topicName::AbstractString
    controlNodes::Array{Symbol}
    treatmentNodes::Array{Symbol}
    offTopicNodes::Array{Symbol}
end

"""
    TopicRotation(
        topicName::AbstractString,
        controlNodes::Array{Symbol},
        treatmentNodes::Array{Symbol}
        offTopicNodes::Array{Symbol}=[]
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

function TopicRotation(
        topicName::AbstractString,
        controlNodes::Array{Symbol},
        treatmentNodes::Array{Symbol}
    )

    return TopicRotation(topicName, controlNodes, treatmentNodes, [])
end

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractTopicRotation, M<:AbstractLinearENAModel{R}}

    edgeIDs = model.edges.edgeID
    embedding = similar(model.embedding, 1)
    if length(model.rotation.offTopicNodes) > 0
        offTopicRows = map(model.nodes.nodeID) do nodeID
            return nodeID in model.rotation.offTopicNodes
        end

        if length(model.rotation.offTopicNodes) == 1
            embedding = similar(model.embedding, 2)
            embedding[1, :label] = string(model.rotation.offTopicNodes[1])
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
            embedding = similar(model.embedding, numSVDDims+length(model.rotation.offTopicNodes)+1)
            embedding[1:numSVDDims, edgeIDs] = svd[1:numSVDDims, :]
            embedding[1:numSVDDims, :label] = ["OffTopic$(i)" for i in 1:numSVDDims]
            embedding[1:numSVDDims, :eigen_value] = eigvals(pca)[1:numSVDDims]
            for i in 1:length(model.rotation.offTopicNodes)
                embedding[i+numSVDDims, :label] = string(model.rotation.offTopicNodes[i])
                for edgeID in edgeIDs
                    muOffTopic = mean(model.nodes[offTopicRows, edgeID])
                    embedding[i+numSVDDims, edgeID] = muOffTopic
                end
            end
        end
    end

    embedding[end, :label] = model.rotation.topicName

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
                embedding[end, edgeID] = muTreatment - muControl
            else
                embedding[end, edgeID] = -muControl
            end
        else
            if sum(treatmentRows) > 0
                muTreatment = mean(model.nodes[treatmentRows, edgeID])
                embedding[end, edgeID] = muTreatment
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

    # let parent handle the rest
    super = rotationsupertype(M, AbstractTopicRotation)
    rotate!(super, model)
end

# flip x and y, sometimes
function defaultplotkwargs(
        ::Type{M},
        model::AbstractLinearENAModel;
        x::Int=1,
        y::Int=2,
        kwargs...
    ) where {R<:AbstractTopicRotation, M<:AbstractLinearENAModel{R}}

    if length(model.rotation.offTopicNodes) == 1
        x = 2
        y = 1
    elseif length(model.rotation.offTopicNodes) > 1
        x = 1 + length(filter(label -> startswith(label, "OffTopic"), model.embedding.label)) + length(model.rotation.offTopicNodes)
        y = 1
    end

    kwargs = NamedTuple(kwargs)
    defaults = (
        x=x,
        y=y,
        kwargs...
    )

    super = rotationsupertype(M, AbstractMeansRotation)
    return defaultplotkwargs(super, model, merge(defaults, kwargs))
end