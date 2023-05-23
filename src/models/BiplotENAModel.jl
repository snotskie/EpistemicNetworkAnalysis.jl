# DONE

# use macro helper to define a standard ENA struct with all the bells
@enamodel BiplotENAModel AbstractLinearENAModel

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