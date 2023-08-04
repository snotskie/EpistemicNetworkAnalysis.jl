abstract type AbstractGroupDifferenceRotation <: AbstractLinearENARotation
    # groupVar::Symbol
end

function test!(
        ::Type{M}, trainmodel::AbstractLinearENAModel, testmodel::AbstractLinearENAModel
    ) where {R<:AbstractGroupDifferenceRotation, M<:AbstractLinearENAModel{R}}

    super = rotationsupertype(M, AbstractGroupDifferenceRotation)
    test!(super, trainmodel, testmodel)

    groups = sort(unique(trainmodel.metadata[!, trainmodel.rotation.groupVar]))
    for i in 1:nrow(testmodel.embedding)
        test!(M, trainmodel, testmodel, KruskalWallisTest, dim=i, groupVar=trainmodel.rotation.groupVar, groups=groups)
    end
end

# insert groups
function defaultplotkwargs(
        ::Type{M},
        model::AbstractLinearENAModel;
        groupBy::Union{Symbol,Nothing}=model.rotation.groupVar,
        kwargs...
    ) where {R<:AbstractGroupDifferenceRotation, M<:AbstractLinearENAModel{R}}

    kwargs = NamedTuple(kwargs)
    defaults = (
        groupBy=groupBy,
        kwargs...
    )

    super = rotationsupertype(M, AbstractGroupDifferenceRotation)
    return defaultplotkwargs(super, model, merge(defaults, kwargs))
end