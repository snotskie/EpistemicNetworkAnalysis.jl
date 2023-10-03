abstract type AbstractTrainedRotation{T<:AbstractLinearENAModel} <: AbstractLinearENARotation end
struct TrainedRotation{T<:AbstractLinearENAModel} <: AbstractTrainedRotation{T}
    trainmodel::T
end

"""
    TrainedRotation(
        trainmodel::AbstractLinearENAModel
    )

Define a rotation that copies the embedding of an existing "training" model for use in a second "testing" model

## Example

```julia
# Randomly assign units to train or test sets
data[!, :rand] = rand(1:100, nrow(data))
trainFilter(unit) = unit.rand < 80
testFilter(unit) = !trainFilter(unit)

# Run model on the training set
trainmodel = ENAModel(
    data, codes, conversations, units,
    unitFilter=trainFilter,
    rotation=LDARotation(:Act)
)

# Run new model on the test set, but using the embedding of the trained model
testmodel = ENAModel(
    data, codes, conversations, units,
    unitFilter=testFilter,
    rotation=TrainedRotation(trainmodel)
)
```

## Statistical Tests

Models using a `TrainedRotation` will run the same statistical tests as the `trainmodel`.
"""
TrainedRotation

# BUGFIX for to_xlsx
NamedTupleTools.fieldvalues(rotation::AbstractTrainedRotation) = [string(nameof(typeof(rotation.trainmodel)))]

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractTrainedRotation, M<:AbstractLinearENAModel{R}}

    embedding = similar(model.embedding, nrow(model.rotation.trainmodel.embedding))
    embedding.label .= model.rotation.trainmodel.embedding.label
    copiedEdges = Symbol.(names(model.rotation.trainmodel.embedding))
    for edgeID in model.edges.edgeID
        if edgeID in copiedEdges
            embedding[!, edgeID] .= model.rotation.trainmodel.embedding[!, edgeID]
        else
            @warn "Edge $(edgeID) does not exist in the copied model's embedding, setting loading to zero instead"
            embedding[!, edgeID] .= 0
        end
    end

    append!(model.embedding, embedding)

    # let parent handle the rest
    super = rotationsupertype(M, AbstractTrainedRotation)
    rotate!(super, model)
end

function substantiate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractTrainedRotation, M<:AbstractLinearENAModel{R}}

    # just copy what the trained model already did
    model.nodes = copy(model.rotation.trainmodel.nodes)
end

function test!(
        ::Type{M}, ::AbstractLinearENAModel, model::AbstractLinearENAModel
    ) where {CR<:AbstractLinearENARotation, CM<:AbstractLinearENAModel{CR}, R<:AbstractTrainedRotation{CM}, M<:AbstractLinearENAModel{R}}

    test!(CM, model.rotation.trainmodel, model)
end

function defaultplotkwargs(
        ::Type{M},
        model::AbstractLinearENAModel;
        x::Int=1,
        y::Int=2,
        meanCenter::Bool=model.config.sphereNormalize,
        origin::Array{<:Real}=(meanCenter ?  [mean(model.points[x, :]), mean(model.points[y, :])] : [0,0]),
        kwargs...
    ) where {R<:AbstractTrainedRotation, M<:AbstractLinearENAModel{R}}
    
    copiedDefaults = defaultplotkwargs(
        typeof(model.rotation.trainmodel),
        model.rotation.trainmodel;
        x=x,
        y=y,
        meanCenter=meanCenter,
        origin=origin,
        kwargs...
    )

    defaults = (
        xlabel="$(copiedDefaults.xlabel)*",
        ylabel="$(copiedDefaults.ylabel)*",
    )
    
    kwargs = NamedTuple(kwargs)
    return merge(copiedDefaults, defaults, kwargs)
end