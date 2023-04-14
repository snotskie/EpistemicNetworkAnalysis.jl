struct DirectionRotation <: AbstractDirectionRotation
    i::Int
 end
 
 # Implement rotation
 function rotate!(rotation::AbstractDirectionRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)
     vecA = Vector{Float64}(unitModel[rotation.i, networkModel[!, :relationship]])
     networkModel[!, :weight_x] .= vecA
     help_one_vector(networkModel, unitModel)
 end