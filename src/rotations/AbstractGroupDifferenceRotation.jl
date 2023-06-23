abstract type AbstractGroupDifferenceRotation <: AbstractLinearENARotation
    # groupVar::Symbol
end

function test!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractGroupDifferenceRotation, M<:AbstractLinearENAModel{R}}

    super = rotationsupertype(M, AbstractGroupDifferenceRotation)
    test!(super, model)

    groups = sort(unique(model.metadata[!, model.rotation.groupVar]))
    for i in 1:nrow(model.embedding)
        test!(M, model, KruskalWallisTest, dim=i, groupVar=model.rotation.groupVar, groups=groups)
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