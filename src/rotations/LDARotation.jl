struct LDARotation <: AbstractLDARotation
    groupVar::Symbol
    dim1::Integer
    dim2::Integer
end

# Simplified constructor
function LDARotation(groupVar::Symbol, dim1::Integer=1)
    return LDARotation(groupVar, dim1, dim1+1)
end

# Implement rotation
function rotate!(rotation::AbstractLDARotation, networkModel::DataFrame, codeModel::DataFrame, metadata::DataFrame, subspaceModel::DataFrame)

    # Check assumptions
    # if nrow(subspaceModel) != nrow(metadata)
    #     error("Cannot perform LDA rotation when rotateOn=:codeModel")
    # end

    # Prepare the data
    groups = sort(unique(metadata[!, rotation.groupVar]))
    groupMap = Dict(group => i for (i, group) in enumerate(groups))
    nc = length(groups)
    X = Matrix{Float64}(transpose(Matrix{Float64}(subspaceModel[!, networkModel[!, :relationship]])))
    for j in 1:size(X, 2)
        X[:, j] = X[:, j] .- mean(X[:, j])
    end

    y = map(metadata[!, rotation.groupVar]) do group
        return groupMap[group]
    end

    ## Run the LDA
    ldaModel = projection(fit(MulticlassLDA, nc, X, y))
    networkModel[!, :weight_x] = ldaModel[:, rotation.dim1]
    if rotation.dim2 == 0
        help_one_vector(networkModel, subspaceModel)
    else
        networkModel[!, :weight_y] = ldaModel[:, rotation.dim2]
        help_two_vectors(networkModel)
    end

    # if size(ldaModel, 2) >= 2
    #     networkModel[!, :weight_x] = ldaModel[:, rotation.dim1]
    #     networkModel[!, :weight_y] = ldaModel[:, rotation.dim2]

    #     ## Normalize the axis weights
    #     s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    #     if s != 0
    #         networkModel[!, :weight_x] /= s
    #     end

    #     s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
    #     if s != 0
    #         networkModel[!, :weight_y] /= s
    #     end
    # else
    #     ## Try to use MR1's x-axis as my approximate y-axis
    #     groups = unique(metadata[!, rotation.groupVar])
    #     rotate!(MeansRotation(rotation.groupVar, groups[1], groups[2]), networkModel, subspaceModel, metadata)
    #     networkModel[!, :weight_y] = networkModel[!, :weight_x]
    #     networkModel[!, :weight_x] = ldaModel[:, 1]
    #     help_two_vectors(networkModel)
    # end
end

# Override plotting pieces
## Base - Inject a groupBy and some labels when none are given
function plot(ena::AbstractENAModel{<:AbstractLDARotation};
    negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
    extraColors::Array{<:Colorant,1}=DEFAULT_EXTRA_COLORS,
    groupBy=nothing,
    xlabel=nothing, ylabel=nothing,
    kwargs...)

    if isnothing(groupBy) || groupBy == ena.rotation.groupVar
        groupBy = ena.rotation.groupVar
    end

    if isnothing(xlabel)
        xlabel = string("LDA", ena.rotation.dim1)
    end

    if isnothing(ylabel)
        if ena.rotation.dim2 == 0
            ylabel = "SVD"
        else
            ylabel = string("LDA", ena.rotation.dim2)
        end
    end

    return invoke(plot, Tuple{AbstractENAModel{<:AbstractLinearENARotation}}, ena;
                  negColor=negColor, posColor=posColor, extraColors=extraColors,
                  groupBy=groupBy, xlabel=xlabel, ylabel=ylabel, kwargs...)
end

# function test(ena::AbstractENAModel{<:AbstractLDARotation})
#     results = invoke(test, Tuple{AbstractENAModel}, ena)

#     ## 99% copy pasta from rotation function, since I don't have the best way to de-dup this work just yet
#     groups = sort(unique(ena.metadata[!, ena.rotation.groupVar]))
#     groupMap = Dict(group => i for (i, group) in enumerate(groups))
#     nc = length(groups)
#     subspaceModel = ena.centroidModel
#     if ena.rotateOn == :accumModel
#         subspaceModel = ena.accumModel
#     end

#     X = Matrix{Float64}(transpose(Matrix{Float64}(subspaceModel[!, ena.networkModel[!, :relationship]])))
#     for j in 1:size(X, 2)
#         X[:, j] = X[:, j] .- mean(X[:, j])
#     end

#     y = map(ena.metadata[!, ena.rotation.groupVar]) do group
#         return groupMap[group]
#     end


#     ## Run SNR test
#     try
#         ldaModel = fit(MulticlassLDA, nc, X, y)
#         W = MultivariateStats.withclass_scatter(ldaModel)
#         B = MultivariateStats.betweenclass_scatter(ldaModel)
#         snr = sum(diag(inv(W) * B))
#         results[:signal_to_noise_ratio] = snr
#     catch e
#         # do nothing
#     end

#     return results
# end
