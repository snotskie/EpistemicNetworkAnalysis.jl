struct MulticlassRotation <: AbstractMulticlassRotation
    groupVar::Symbol
    dim1::Integer
    dim2::Integer
end

# Simplified constructor
function MulticlassRotation(groupVar::Symbol, dim1::Integer=1)
    return MulticlassRotation(groupVar, dim1, dim1+1)
end

# Implement rotation
function rotate!(rotation::AbstractMulticlassRotation, networkModel::DataFrame, codeModel::DataFrame, metadata::DataFrame, subspaceModel::DataFrame)

    # Prepare the data
    ## Helpers
    relationshipIds = networkModel[!, :relationship]
    
    ## Groups
    groupNames = sort(unique(metadata[!, rotation.groupVar]))

    ## Mean Vectors
    all_rows = map(eachrow(subspaceModel)) do row
        return false
    end

    group_ns = []
    group_vecs = []
    for g in groupNames
        group_rows = metadata[!, rotation.groupVar] .== g
        all_rows = all_rows .| group_rows
        group_vec = Vector{Float64}(mean.(subspaceModel[group_rows, col] for col in relationshipIds))
        push!(group_vecs, group_vec)
        push!(group_ns, sum(group_rows))
    end

    mean_vec = Vector{Float64}(mean.(subspaceModel[all_rows, col] for col in relationshipIds))

    ## Offset Vectors
    offset_vecs = []
    for vec in group_vecs
        push!(offset_vecs, vec - mean_vec)
    end

    # Run Computations
    ## Between-Group Scatter Matrix
    Sb = sum(
        n * vec * transpose(vec)
        for (n, vec) in zip(group_ns, offset_vecs) 
    )

    ## Solve Eigenvalue Problem
    vals = eigvals(Sb)
    vecs = eigvecs(Sb)

    ## Store Axes
    R = real.(vecs[:, end-rotation.dim1+1]) # vecs are stored column-major, not row major
    R /= sqrt(sum(R .^ 2))
    networkModel[!, :weight_x] .= R
    if rotation.dim2 == 0
        help_one_vector(networkModel, subspaceModel)
    else
        R = real.(vecs[:, end-rotation.dim2+1]) # vecs are stored column-major, not row major
        R /= sqrt(sum(R .^ 2))
        networkModel[!, :weight_y] .= R
        help_two_vectors(networkModel)
    end
        
end

# Override plotting pieces
## Base - Inject a groupBy and some labels when none are given
function plot(ena::AbstractENAModel{<:AbstractMulticlassRotation};
    negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
    extraColors::Array{<:Colorant,1}=DEFAULT_EXTRA_COLORS,
    groupBy=nothing,
    xlabel=nothing, ylabel=nothing,
    kwargs...)

    if isnothing(groupBy) || groupBy == ena.rotation.groupVar
        groupBy = ena.rotation.groupVar
    end

    if isnothing(xlabel)
        xlabel = string("MCR", ena.rotation.dim1)
    end

    if isnothing(ylabel)
        if ena.rotation.dim2 == 0
            ylabel = "SVD"
        else
            ylabel = string("MCR", ena.rotation.dim2)
        end
    end

    return invoke(plot, Tuple{AbstractENAModel{<:AbstractLinearENARotation}}, ena;
                  negColor=negColor, posColor=posColor, extraColors=extraColors,
                  groupBy=groupBy, xlabel=xlabel, ylabel=ylabel, kwargs...)
end
