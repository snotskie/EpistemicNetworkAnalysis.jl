struct SVDRotation <: AbstractSVDRotation
    # no fields
end

# Implement rotation
function rotate!(rotation::AbstractSVDRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)

    ## Run an ortho svd and use those values as the axis weights
    pcaModel = projection(help_deflating_svd(networkModel, unitModel))
    networkModel[!, :weight_x] = pcaModel[:, 1]
    networkModel[!, :weight_y] = pcaModel[:, 2]
end

# Override plotting pieces
## Inject labels when none are given
function plot(ena::AbstractENAModel{<:AbstractSVDRotation};
    xlabel=nothing, ylabel=nothing,
    kwargs...)

    if isnothing(xlabel)
        xlabel = "SVD1"
    end

    if isnothing(ylabel)
        ylabel = "SVD2"
    end

    return invoke(plot, Tuple{AbstractENAModel{<:AbstractENARotation}}, ena;
                  xlabel=xlabel, ylabel=ylabel, kwargs...)
end