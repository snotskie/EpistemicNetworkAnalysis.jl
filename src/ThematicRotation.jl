struct ThematicRotation <: AbstractThematicRotation
    controlThemes::Array{Symbol}
    treatmentThemes::Array{Symbol}
 end
 
 # Implement rotation
 function rotate!(rotation::AbstractThematicRotation, networkModel::DataFrame, codeModel::DataFrame, metadata::DataFrame, subspaceModel::DataFrame)
    controlRows = map(eachrow(codeModel)) do codeRow
        if codeRow[:code] in rotation.controlThemes
            return true
        else
            return false
        end
    end
    
    treatmentRows = map(eachrow(codeModel)) do codeRow
        if codeRow[:code] in rotation.treatmentThemes
            return true
        else
            return false
        end
    end

    controlMean = combine(codeModel[controlRows, :], (r => mean => r for r in networkModel[!, :relationship])...)
    treatmentMean = combine(codeModel[treatmentRows, :], (r => mean => r for r in networkModel[!, :relationship])...)
    controlVec = Vector{Float64}(controlMean[1, networkModel[!, :relationship]])
    treatmentVec = Vector{Float64}(treatmentMean[1, networkModel[!, :relationship]])
    networkModel[!, :weight_x] .= treatmentVec .- controlVec
    help_one_vector(networkModel, subspaceModel)
 end
 
# Override plotting pieces
## Base - Inject a groupBy and some labels when none are given
function plot(ena::AbstractENAModel{<:AbstractThematicRotation};
    xlabel=nothing, ylabel=nothing,
    kwargs...)

    if isnothing(xlabel)
        left = join(ena.rotation.controlThemes, "/")
        if length(left) > 15
            left = join(map(x->string(x)[1:min(5,end)], ena.rotation.controlThemes), "/")
        end

        right = join(ena.rotation.treatmentThemes, "/")
        if length(right) > 15
            right = join(map(x->string(x)[1:min(5,end)], ena.rotation.treatmentThemes), "/")
        end

        xlabel = "$(left) to $(right)"
        # xlabel = "$( to $(join(ena.rotation.treatmentThemes, '/'))"
    end

    if isnothing(ylabel)
        ylabel = "SVD"
    end

    return invoke(plot, Tuple{AbstractENAModel{<:AbstractLinearENARotation}}, ena;
                xlabel=xlabel, ylabel=ylabel, kwargs...)
end