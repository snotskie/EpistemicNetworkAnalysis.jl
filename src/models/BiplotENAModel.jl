# DONE

# use macro helper to define a standard ENA struct with all the bells
@enamodel BiplotENAModel AbstractLinearENAModel

# override default model constructor kwargs
function defaultmodelkwargs(
        ::Type{M};
        kwargs...
    ) where {R<:AbstractLinearENARotation, M<:AbstractBiplotENAModel{R}}

    kwargs = NamedTuple(kwargs)
    super = modelsupertype(M, AbstractBiplotENAModel)
    parentdefaults = defaultmodelkwargs(super; kwargs...)
    definitivedefaults = (
        edgeFilter=(row)->(
            row[:kind] == :count #&& parentdefaults.edgeFilter(row)
        ),
        windowSize=1
    )

    return merge(parentdefaults, kwargs, definitivedefaults)
end

# let the parent handle it from there