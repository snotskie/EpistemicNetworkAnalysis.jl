struct SVDRotation <: AbstractSVDRotation
    dim1::Integer
    dim2::Integer
end

# Simplified constructor
function SVDRotation(dim1::Integer=1)
    return SVDRotation(dim1, dim1+1)
end

# Implement rotation
function rotate!(rotation::AbstractSVDRotation, networkModel::DataFrame, codeModel::DataFrame, metadata::DataFrame, subspaceModel::DataFrame)

    ## Run an ortho svd and use those values as the axis weights
    pcaModel = projection(help_deflating_svd(networkModel, subspaceModel))
    networkModel[!, :weight_x] = pcaModel[:, rotation.dim1]
    networkModel[!, :weight_y] = pcaModel[:, rotation.dim2]
end

# Override plotting pieces
## Inject labels when none are given
function plot(ena::AbstractENAModel{<:AbstractSVDRotation};
    xlabel=nothing, ylabel=nothing,
    kwargs...)

    if isnothing(xlabel)
        xlabel = string("SVD", ena.rotation.dim1)
    end

    if isnothing(ylabel)
        ylabel = string("SVD", ena.rotation.dim2)
    end

    return invoke(plot, Tuple{AbstractENAModel{<:AbstractENARotation}}, ena;
                  xlabel=xlabel, ylabel=ylabel, kwargs...)
end