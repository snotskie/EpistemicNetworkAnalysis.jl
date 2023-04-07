# DONE

# use macro helper to define a standard ENA struct with bells
@enamodel BiplotENAModel LinearENAModel SVDRotation

# override default model constructor kwargs
function defaultmodelkwargs(AbstractBiplotENAModel; kwargs...)
    parentdefaults = defaultmodelkwargs(supertype(AbstractBiplotENAModel); kwargs...)
    defaults = (
        relationshipFilter=(i, j, ci, cj)->(i == j),
        windowSize=1
    )

    return merge(parentdefaults, defaults, kwargs)
end

# let the parent handle it from there