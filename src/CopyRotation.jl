struct CopyRotation <: AbstractCopyRotation
    ena::AbstractENAModel
end

# Accessor pass through
function Base.getproperty(rotation::AbstractCopyRotation, field::Symbol)
    if field == :ena
        return getfield(rotation, field)
    else
        return getfield(rotation.ena.rotation, field)
    end
end

# Implement rotation
function rotate!(rotation::AbstractCopyRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)

    # check assumptions
    if nrow(networkModel) != nrow(rotation.ena.networkModel)
        error("Cannot copy rotation from an ENA model with a different network size.")
    end

    # copy the rotation
    networkModel[!, :weight_x] = copy(rotation.ena.networkModel[!, :weight_x])
    networkModel[!, :weight_y] = copy(rotation.ena.networkModel[!, :weight_y])
end

# Plot pass through
# TODO the type conversion is throwing an error
# function plot(ena::AbstractENAModel{<:AbstractCopyRotation}; kwargs...)
#     return invoke(plot, Tuple{typeof(ena.rotation.ena)}, ena.rotation.ena, ena); kwargs)
# end