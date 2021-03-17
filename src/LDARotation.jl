struct LDARotation <: AbstractLDARotation
    groupVar::Symbol
end

# Implement rotation
function rotate!(rotation::AbstractLDARotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)

    # Check assumptions
    if nrow(unitModel) != nrow(metadata)
        error("Cannot perform LDA rotation when rotateOn=:codeModel")
    end

    # Prepare the data
    groups = sort(unique(metadata[!, rotation.groupVar]))
    groupMap = Dict(group => i for (i, group) in enumerate(groups))
    nc = length(groups)
    X = Matrix{Float64}(transpose(Matrix{Float64}(unitModel[!, networkModel[!, :relationship]])))
    y = map(metadata[!, rotation.groupVar]) do group
        return groupMap[group]
    end

    ## Run the LDA
    ldaModel = projection(fit(MulticlassLDA, nc, X, y))
    if size(ldaModel, 2) >= 2
        networkModel[!, :weight_x] = ldaModel[:, 1]
        networkModel[!, :weight_y] = ldaModel[:, 2]

        ## Normalize the axis weights
        s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
        if s != 0
            networkModel[!, :weight_x] /= s
        end

        s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
        if s != 0
            networkModel[!, :weight_y] /= s
        end
    else
        ## Use MR1's x-axis as my approximate y-axis
        groups = unique(metadata[!, rotation.groupVar])
        rotate!(MeansRotation(rotation.groupVar, groups[1], groups[2]), networkModel, unitModel, metadata)
        networkModel[!, :weight_y] = networkModel[!, :weight_x]
        networkModel[!, :weight_x] = ldaModel[:, 1]
        help_two_vectors(networkModel)
    end
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
        xlabel = "LDA1"
    end

    if isnothing(ylabel)
        ylabel = "LDA2"
    end

    return invoke(plot, Tuple{AbstractENAModel{<:AbstractENARotation}}, ena;
                  negColor=negColor, posColor=posColor, extraColors=extraColors,
                  groupBy=groupBy, xlabel=xlabel, ylabel=ylabel, kwargs...)
end