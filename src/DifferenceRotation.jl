struct DifferenceRotation <: AbstractDifferenceRotation
   i::Int
   j::Int
end

# Implement rotation
function rotate!(rotation::AbstractDifferenceRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)
    vecA = Vector{Float64}(unitModel[rotation.i, networkModel[!, :relationship]])
    vecB = Vector{Float64}(unitModel[rotation.j, networkModel[!, :relationship]])
    networkModel[!, :weight_x] .= vecA - vecB
    help_one_vector(networkModel, unitModel)
end