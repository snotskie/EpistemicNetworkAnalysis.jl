"""
TODO document
"""
struct LDARotation <: AbstractLDARotation
    groupVar::Symbol
end

# Implement rotation
function rotate!(rotation::AbstractLDARotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)

    # Prepare the data
    groups = sort(unique(metadata[!, rotation.groupVar]))
    groupMap = Dict(group => i for (i, group) in enumerate(groups))
    nc = length(groups)
    X = Matrix{Float64}(transpose(Matrix{Float64}(unitModel[!, networkModel[!, :relationship]])))
    y = map(metadata[!, rotation.groupVar]) do group
        return groupMap[group]
    end

    # y = Vector{Float64}(y)

    ## Run the LDA
    ldaModel = projection(fit(MulticlassLDA, nc, X, y))
    networkModel[!, :weight_x] = ldaModel[:, 1]
    if size(ldaModel, 2) >= 2
        networkModel[!, :weight_y] = ldaModel[:, 2]
    else
        ## Find the first svd dim of the data orthogonal to the x weights, use these as the y weights
        xAxis = Matrix{Float64}(unitModel[!, networkModel[!, :relationship]]) *
        Matrix{Float64}(networkModel[!, [:weight_x]])
        xAxis = xAxis .- mean(xAxis)
        controlModel = DataFrame(xAxis)
        pcaModel = projection(help_deflating_svd(networkModel, unitModel, controlModel))
        networkModel[!, :weight_y] = pcaModel[:, 1]
    end

    ## Normalize the axis weights
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
    if s != 0
        networkModel[!, :weight_y] /= s
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