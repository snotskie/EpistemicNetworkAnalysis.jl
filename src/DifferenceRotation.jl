struct DifferenceRotation <: AbstractDifferenceRotation
   i::Int
   j::Int
end

# Implement rotation
function rotate!(rotation::AbstractDifferenceRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)
    vecA = Vector{Float64}(unitModel[rotation.i, networkModel[!, :relationship]])
    vecB = Vector{Float64}(unitModel[rotation.j, networkModel[!, :relationship]])
    networkModel[!, :weight_x] .= vecB - vecA
    help_one_vector(networkModel, unitModel)
end

# Override plotting pieces
## Base - Inject a groupBy and some labels when none are given
function plot(ena::AbstractENAModel{<:AbstractDifferenceRotation};
    xlabel=nothing, ylabel=nothing,
    kwargs...)

    if isnothing(xlabel)
        if ena.rotateOn == :codeModel
            xlabel = "$(ena.codeModel[ena.rotation.i, :code]) to $(ena.codeModel[ena.rotation.j, :code])"
        else
            xlabel = "$(ena.metadata[ena.rotation.i, :ENA_UNIT]) to $(ena.metadata[ena.rotation.j, :ENA_UNIT])"
        end
    end

    if isnothing(ylabel)
        ylabel = "SVD"
    end

    return invoke(plot, Tuple{AbstractENAModel{<:AbstractENARotation}}, ena;
                  xlabel=xlabel, ylabel=ylabel, kwargs...)
end