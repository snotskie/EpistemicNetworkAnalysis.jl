# DONE

# use macro helper to define a standard ENA struct with all the bells
@enamodel BiplotENAModel AbstractLinearENAModel

# Documentation
"""
    BiplotENAModel(
        # Required
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1};

        # Optional
        rotation::AbstractLinearENARotation=SVDRotation(),
        unitFilter::Function=unit->true,
        sphereNormalize::Bool=true,
        dropEmpty::Bool=false,
        recenterEmpty::Bool=false
    )

Construct a biplot model of unit-wise counts of code occurences, without measuring connections between codes. Model will have perfect goodness of fit between `points` and `pointsHat`, will be much simpler than other model types, but will lose most information compared to other model types.

`BiplotENAModel` follows the same argument and field structure as `ENAModel`, except `edgeFilter` and `windowSize` are in effect ignored.
"""
BiplotENAModel

# override default model constructor kwargs
function defaultmodelkwargs(
        ::Type{M};
        prev_config::NamedTuple=NamedTuple(),
        kwargs...
    ) where {R<:AbstractLinearENARotation, M<:AbstractBiplotENAModel{R}}

    kwargs = NamedTuple(kwargs)
    super = modelsupertype(M, AbstractBiplotENAModel)
    parentdefaults = defaultmodelkwargs(super)
    definitivedefaults = (
        edgeFilter=(row)->(
            row[:kind] == :count
        ),
        windowSize=1
    )

    return merge(parentdefaults, prev_config, definitivedefaults, kwargs)
end

# let the parent handle it from there