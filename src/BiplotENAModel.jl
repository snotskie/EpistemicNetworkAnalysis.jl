# DONE

# use macro helper to define a standard ENA struct with all the bells
@enamodel(BiplotENAModel, AbstractLinearENAModel)

# override default model constructor kwargs
function defaultmodelkwargs(
        T::Type{M{R}};
        kwargs...
    ) where {M<:AbstractBiplotENAModel,R<:AbstractLinearENARotation}

    super = supertype(AbstractBiplotENAModel){R}
    parentdefaults = defaultmodelkwargs(super; kwargs...)
    defaults = (
        edgeFilter=(row)->(row[:kind] == :count),
        windowSize=1
    )

    return merge(parentdefaults, defaults, kwargs)
end

# let the parent handle it from there