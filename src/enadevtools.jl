# Type Tree

## Rotations
abstract type AbstractENARotation end
abstract type AbstractLinearENARotation <: AbstractENARotation end
# abstract type AbstractNonlinearENARotation <: AbstractENARotation end

## Models
abstract type AbstractENAModel{T<:AbstractENARotation} end
abstract type AbstractLinearENAModel{T<:AbstractLinearENARotation} <: AbstractENAModel{T} end
# abstract type AbstractNonlinearENAModel{T<:AbstractNonlinearENARotation} <: AbstractENAModel{T} end

## EdgePainters
abstract type AbstractEdgePainter end

## Type Helpers
function modelsupertype(::Type{M}, ::Type{N}) where {
        R<:AbstractENARotation,
        M<:AbstractENAModel{R},
        N<:AbstractENAModel
    }
    
    return supertype(N){R}
end

function rotationsupertype(::Type{M}, ::Type{S}) where {
        R<:AbstractENARotation,
        M<:AbstractENAModel{R},
        S<:AbstractENARotation
    }
    
    return M.name.wrapper{supertype(S)}
end

## Type Macros
macro enamodel(
        self, parent,
        defaultrotation=:SVDRotation,
        rotationtype=:AbstractLinearENARotation,
        abstractname=Symbol(string(:Abstract, self))
    )

    return quote
        # make abstract type
        abstract type $(esc(abstractname)){T<:$(esc(rotationtype))} <: $(esc(parent)){T} end

        # make struct
        mutable struct $(esc(self)){T<:$(esc(rotationtype))} <: $(esc(abstractname)){T}
            # required arguments
            data::DataFrame
            codes::Array{Symbol,1}
            conversations::Array{Symbol,1}
            units::Array{Symbol,1}
        
            # rotation
            rotation::T
        
            # model
            metadata::DataFrame
            points::DataFrame
            pointsHat::DataFrame
            pointsNodes::DataFrame
            accum::DataFrame
            accumHat::DataFrame
            edges::DataFrame
            nodes::DataFrame
            embedding::DataFrame
            config::NamedTuple
        end

        # make constructor
        function $(esc(self))(
                data::DataFrame,
                codes::Array{<:Any,1},
                conversations::Array{<:Any,1},
                units::Array{<:Any,1};
                rotateBy::R=$(esc(defaultrotation))(),
                kwargs...
            ) where {R<:$(esc(rotationtype))}

            # call common ENA constructor
            return constructENA(
                $(esc(self)){R},
                data, codes, conversations, units,
                rotateBy; kwargs...
            )
        end
    end
end

# Wrapper constructor, accepts any type for codes etc.
function constructENA(
        ::Type{M},
        data::DataFrame,
        codes::Array{<:Any,1},
        conversations::Array{<:Any,1},
        units::Array{<:Any,1},
        rotation::AbstractENARotation;
        kwargs...
    ) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    return constructENA(
        M,
        data,
        Symbol.(codes),
        Symbol.(conversations),
        Symbol.(units),
        rotation;
        kwargs...
    )
end

# Base constructor, requires symbols
function constructENA(
        ::Type{M},
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1},
        rotation::AbstractENARotation;
        kwargs...
    ) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}

    kwargs = defaultmodelkwargs(M; kwargs...)
    model = M(
        populateENAfields(
            M, data, codes, conversations, units, rotation;
            kwargs...
        )...
    )

    accumulate!(M, model)
    approximate!(M, model)
    rotate!(M, model)
    return model
end

## Unimplemented Functions
function populateENAfields(
        ::Type{M},
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1},
        rotation::AbstractENARotation;
        config...
    ) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    
    error("Unimplemented")
end

function defaultmodelkwargs(::Type{M}, model::AbstractENAModel; kwargs...) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function defaultplotkwargs(::Type{M}, model::AbstractENAModel; kwargs...) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

# NOTE: when implementing these functions elsewhere, M should be the *most* specific
# type that the function applies to, while model should be the *least* specific.
# Also, accumulate! and approximate! should be generic to the rotation type when possible,
# while rotate! should be generic to the ENA type when possible.
function accumulate!(::Type{M}, model::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function approximate!(::Type{M}, model::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function rotate!(::Type{M}, model::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function summarize!(::Type{M}, model::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

# in linear, do plot like the consruct helper, override its components under there
function plot(::Type{M}, model::AbstractENAModel, plotconfig::NamedTuple) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

## Wrapper Functions
function plot(model::AbstractENAModel; kwargs...)
    plotconfig = defaultplotkwargs(typeof(model), model; kwargs...)
    return plot(typeof(model), model, plotconfig)
end

# Helpers

function computeNetworkDensities(model, rows=!; normalize=false)

    # sum densities for each edge
    edgeIDs = model.edges.edgeID
    edgeDensities = Dict(
        edgeID => sum(model.accum[rows, edgeID])
        for edgeID in edgeIDs
    )

    # find the density of each code row dot, by "splitting" the density of each line between its two codes
    nodeIDs = model.nodes.nodeID
    nodeDensities = Dict(
        nodeID => 0.0
        for nodeID in nodeIDs
    )

    for edge in eachrow(model.edges)
        nodeDensities[edge.ground] += edgeDensities[edge.edgeID]
        nodeDensities[edge.response] += edgeDensities[edge.edgeID]
    end

    # optionally normalize each
    if normalize
        s = maximum(values(edgeDensities))
        if s != 0
            for edgeID in edgeIDs
                edgeDensities[edge] /= s
            end
        end

        s = maximum(values(nodeDensities))
        if s != 0
            for nodeID in nodeIDs
                nodeDensities[node] /= s
            end
        end
    end

    return edgeDensities, nodeDensities
end

function addPointsToModelFromDim(model, dim)
    edgeIDs = model.edges.edgeID
    unitIDs = model.accum.unitID
    nodeIDs = model.nodes.nodeID
    points = Matrix{Float64}(model.accum[!, edgeIDs]) * Vector{Float64}(dim[edgeIDs])
    df = similar(model.points, 1)
    df[1, unitIDs] = points
    append!(model.points, df)
    pointsHat = Matrix{Float64}(model.accumHat[!, edgeIDs]) * Vector{Float64}(dim[edgeIDs])
    df = similar(model.pointsHat, 1)
    df[1, unitIDs] = pointsHat
    append!(model.pointsHat, df)
    pointsNodes = Matrix{Float64}(model.nodes[!, edgeIDs]) * Vector{Float64}(dim[edgeIDs])
    df = similar(model.pointsNodes, 1)
    df[1, nodeIDs] = pointsNodes
    append!(model.pointsNodes, df)
end

function layoutSubplots(ps::Array{Plot}, plotconfig::NamedTuple)
    for p in ps
        xticks!(p, plotconfig.xticks)
        yticks!(p, plotconfig.yticks)
        if plotconfig.lims > 0
            xlims!(p, plotconfig.xlims...)
            ylims!(p, plotconfig.ylims...)
        end
    end

    N = ceil(Int, sqrt(length(ps)))
    M = ceil(Int, length(ps)/N)
    layout = grid(N, M)
    while length(ps) < N*M
        push!(ps, plot(legend=false,grid=false,foreground_color_subplot=:white))
    end

    return plot(ps..., size=(plotconfig.size*M, plotconfig.size*N), layout=layout)
end

# function help_deflating_svd(networkModel::DataFrame, subspaceModel::DataFrame, controlModel::Union{Nothing,DataFrame}=nothing)
#     X = Matrix{Float64}(subspaceModel[!, networkModel[!, :relationship]])
#     for i in 1:size(X)[2]
#         xcol = X[:, i]
#         xcol = xcol .- mean(xcol) # mean center
#         X[:, i] = xcol
#     end

#     if !isnothing(controlModel)
#         C = Matrix{Float64}(controlModel)
#         for i in 1:size(C)[2]
#             ccol = C[:, i]
#             ccol = ccol .- mean(ccol) # mean center
#             C[:, i] = ccol
#             for j in 1:size(X)[2] # deflate
#                 xcol = X[:, j]
#                 scalar = dot(xcol, ccol) / dot(ccol, ccol)
#                 xcol -= scalar * ccol
#                 X[:, j] = xcol
#             end
#         end
#     end

#     # then, once we've deflated or not, we run an SVD on the data
#     pcaModel = fit(PCA, X', pratio=1.0)
#     return pcaModel
# end

# function help_one_vector(networkModel::DataFrame, subspaceModel::DataFrame)
#     ## Normalize the axis weights
#     s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
#     if s != 0
#         networkModel[!, :weight_x] /= s
#     end

#     ## Find the first svd dim of the data orthogonal to the x weights, use these as the y weights
#     xAxis = Matrix{Float64}(subspaceModel[!, networkModel[!, :relationship]]) *
#             Matrix{Float64}(networkModel[!, [:weight_x]])
#     xAxis = xAxis .- mean(xAxis)
#     controlModel = DataFrame(xAxis, :auto)
#     pcaModel = projection(help_deflating_svd(networkModel, subspaceModel, controlModel))
#     networkModel[!, :weight_y] = pcaModel[:, 1]
# end

# function help_two_vectors(networkModel::DataFrame)
#     ## Normalize the weights for both axes
#     s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
#     if s != 0
#         networkModel[!, :weight_x] /= s
#     end

#     s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
#     if s != 0
#         networkModel[!, :weight_y] /= s
#     end

#     ## Orthogonalization: replace y weights with their rejection from the x weights
#     before = copy(networkModel[!, :weight_y])
#     scalar = dot(networkModel[!, :weight_y], networkModel[!, :weight_x]) / dot(networkModel[!, :weight_x], networkModel[!, :weight_x])
#     networkModel[!, :weight_y] -= scalar * networkModel[!, :weight_x]
#     after = copy(networkModel[!, :weight_y])

#     ## Raise a warning about interpreting the y-axis when before and after have a large angle between them
#     theta = dot(before, after)
#     theta /= sqrt(dot(before, before))
#     theta /= sqrt(dot(after, after))
#     if theta > 1 # bugfix for rounding error
#         theta = 1
#     elseif theta < -1
#         theta = -1
#     end
#     angle = acos(theta) * 180 / pi
#     if abs(angle) > 5
#         @warn """The angle between the y-axis and the direction of the requested effect is larger than 5 degrees ($angle degrees).
# This can undermine interpreting the y-axis in terms of the requested effect."""
#     end

#     ## Re-normalize the weights for the y-axis
#     s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
#     if s < 0.05
#         networkModel[!, :weight_y] .= 0
#         @warn "During the rotation step, the y axis was deflated to zero due to close correlation with the x axis."
#     elseif s != 0
#         networkModel[!, :weight_y] /= s
#     end
# end

# function help_xs_and_ys(ena, displayRows, flipX::Bool, flipY::Bool)
#     unitModel = ena.accumModel
#     xs = unitModel[displayRows, :pos_x] * (flipX ? -1 : 1)
#     ys = unitModel[displayRows, :pos_y] * (flipY ? -1 : 1)
#     return (xs, ys)
# end

function nonlinearGradientMap(lo, mid, hi; grains=100, curve=1.5)
    return vcat(
        [weighted_color_mean((100-i)^curve/grains^curve, lo, mid) for i in 1:grains],
        [weighted_color_mean(1-i^curve/grains^curve, mid, hi) for i in 1:grains]
    )
end

# function rotatedLabel(label, x, y)
#     angle = atan(y, x) * 180 / pi
#     if x < 0
#         angle = atan(-y, -x) * 180 / pi
#     else
    
#     return text(label, :top, default(:xtickfontsize), rotation=angle)
# end