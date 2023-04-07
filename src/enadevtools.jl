module enadevtools

# Type Tree

## Rotations
abstract type AbstractENARotation end
abstract type AbstractLinearENARotation <: AbstractENARotation end
# abstract type AbstractNonlinearENARotation <: AbstractENARotation end

## Models
abstract type AbstractENAModel{T<:AbstractENARotation} end
abstract type AbstractLinearENAModel{T<:AbstractLinearENARotation} <: AbstractENAModel{T} end
# abstract type AbstractNonlinearENAModel{T<:AbstractNonlinearENARotation} <: AbstractENAModel{T} end

## Type Macros
macro enamodel(
        self, parent;
        defaultrotation=:SVDRotation,
        rotationtype=:AbstractLinearENARotation
        abstractname=Symbol(string(:Abstract, self))
    )

    return quote
        # make abstract type
        abstract type $abstractname{T<:$rotationtype} <: $parent{T<:$rotationtype} end

        # make struct
        struct $self{T<:$rotationtype} <: $abstractname{T<:$rotationtype}
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
        function $self(
                data::DataFrame,
                codes::Array{Symbol,1},
                conversations::Array{Symbol,1},
                units::Array{Symbol,1};
                rotateBy::T=$defaultrotation(),
                kwargs...
            ) where T <: $rotationtype

            # call common ENA constructor
            return constructENA(
                $self{T},
                data, codes, conversations, units,
                rotateBy, kwargs...
            )
        end
    end
end

return constructENA(
        T::Type{M{R}},
        data, codes, conversations, units, rotation;
        kwargs...
    ) where {M<:AbstractENAModel,R<:AbstractENARotation}

    kwargs = defaultmodelkwargs(T; kwargs...)
    model = T(
        data, codes, conversations, units, rotation,
        populateENAfields(
            T, data, codes, conversations, units, rotation;
            kwargs...
        )...
    )

    accumulate!(T, model)
    approximate!(T, model)
    rotate!(T, model)
    return model
end

## Unimplemented Functions
function populateENAfields(
        T::Type{M{R}},
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1},
        rotation::R;
        config...
    ) where {M<:AbstractENAModel,R<:AbstractENARotation}
    
    error("Unimplemented")
end

function accumulate!(T::Type{M{R}}, model::T) where {M<:AbstractENAModel,R<:AbstractENARotation}
    error("Unimplemented")
end

function approximate!(T::Type{M{R}}, model::T) where {M<:AbstractENAModel,R<:AbstractENARotation}
    error("Unimplemented")
end

function rotate!(T::Type{M{R}}, model::T) where {M<:AbstractENAModel,R<:AbstractENARotation}
    error("Unimplemented")
end

function tests(T::Type{M{R}}, model::T) where {M<:AbstractENAModel,R<:AbstractENARotation}
    error("Unimplemented")
end

# in linear, do plot like the consruct helper, override its components under there
function plot(T::Type{M{R}}, model::T; kwargs...) where {M<:AbstractENAModel,R<:AbstractENARotation}
    error("Unimplemented")
end

## Wrapper Functions
function plot(model::AbstractENAModel; kwargs...)
    return plot(type(model), model; kwargs...)
end

# TODO below here

# Helpers

function computeNetworkDensities(model, rows=!)

    # sum densities for each edge
    edgeNames = model.edges.edge
    edgeDensities = sum(eachcol(model.accum[rows, edgeNames]))

    # find the density of each code row dot, by "splitting" the density of each line between its two codes
    nodeNames = model.nodes[!, :node]
    nodeDensities = zeros(length(nodeNames))
    for (k, edge) in eachrow(model.edges)
        i, j = edge.i, edge.j
        nodeDensities[i] += edgeDensities[k]
        nodeDensities[j] += edgeDensities[k]
    end

    # normalize each
    s = maximum(edgeDensities)
    if s != 0
        edgeDensities /= s
    end

    s = maximum(nodeDensities)
    if s != 0
        nodeDensities /= s
    end

    edgeDensityDict = Dict(
        edge => density
        for (edge, density) in zip(edgeNames, edgeDensities)
    )

    nodeDensityDict = Dict(
        node => density
        for (node, density) in zip(nodeNames, nodeDensities)
    )

    return edgeDensityDict, nodeDensityDict
end


function help_deflating_svd(networkModel::DataFrame, subspaceModel::DataFrame, controlModel::Union{Nothing,DataFrame}=nothing)
    X = Matrix{Float64}(subspaceModel[!, networkModel[!, :relationship]])
    for i in 1:size(X)[2]
        xcol = X[:, i]
        xcol = xcol .- mean(xcol) # mean center
        X[:, i] = xcol
    end

    if !isnothing(controlModel)
        C = Matrix{Float64}(controlModel)
        for i in 1:size(C)[2]
            ccol = C[:, i]
            ccol = ccol .- mean(ccol) # mean center
            C[:, i] = ccol
            for j in 1:size(X)[2] # deflate
                xcol = X[:, j]
                scalar = dot(xcol, ccol) / dot(ccol, ccol)
                xcol -= scalar * ccol
                X[:, j] = xcol
            end
        end
    end

    # then, once we've deflated or not, we run an SVD on the data
    pcaModel = fit(PCA, X', pratio=1.0)
    return pcaModel
end

function help_one_vector(networkModel::DataFrame, subspaceModel::DataFrame)
    ## Normalize the axis weights
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    ## Find the first svd dim of the data orthogonal to the x weights, use these as the y weights
    xAxis = Matrix{Float64}(subspaceModel[!, networkModel[!, :relationship]]) *
            Matrix{Float64}(networkModel[!, [:weight_x]])
    xAxis = xAxis .- mean(xAxis)
    controlModel = DataFrame(xAxis, :auto)
    pcaModel = projection(help_deflating_svd(networkModel, subspaceModel, controlModel))
    networkModel[!, :weight_y] = pcaModel[:, 1]
end

function help_two_vectors(networkModel::DataFrame)
    ## Normalize the weights for both axes
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
    if s != 0
        networkModel[!, :weight_y] /= s
    end

    ## Orthogonalization: replace y weights with their rejection from the x weights
    before = copy(networkModel[!, :weight_y])
    scalar = dot(networkModel[!, :weight_y], networkModel[!, :weight_x]) / dot(networkModel[!, :weight_x], networkModel[!, :weight_x])
    networkModel[!, :weight_y] -= scalar * networkModel[!, :weight_x]
    after = copy(networkModel[!, :weight_y])

    ## Raise a warning about interpreting the y-axis when before and after have a large angle between them
    theta = dot(before, after)
    theta /= sqrt(dot(before, before))
    theta /= sqrt(dot(after, after))
    if theta > 1 # bugfix for rounding error
        theta = 1
    elseif theta < -1
        theta = -1
    end
    angle = acos(theta) * 180 / pi
    if abs(angle) > 5
        @warn """The angle between the y-axis and the direction of the requested effect is larger than 5 degrees ($angle degrees).
This can undermine interpreting the y-axis in terms of the requested effect."""
    end

    ## Re-normalize the weights for the y-axis
    s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
    if s < 0.05
        networkModel[!, :weight_y] .= 0
        @warn "During the rotation step, the y axis was deflated to zero due to close correlation with the x axis."
    elseif s != 0
        networkModel[!, :weight_y] /= s
    end
end

function help_xs_and_ys(ena, displayRows, flipX::Bool, flipY::Bool)
    unitModel = ena.accumModel
    xs = unitModel[displayRows, :pos_x] * (flipX ? -1 : 1)
    ys = unitModel[displayRows, :pos_y] * (flipY ? -1 : 1)
    return (xs, ys)
end

function help_plot_ci(p, xs, ys, color, shape, label, showCIs::Bool=true)
    if length(xs) > 0
        x = mean(xs)
        y = mean(ys)
        Plots.plot!(p, [x], [y],
            label=label,
            seriestype=:scatter,
            markersize=4,
            markershape=shape,
            markercolor=color,
            markerstrokecolor=color)
    end

    if length(xs) > 1 && showCIs
        ci_x = collect(confint(OneSampleTTest(xs)))
        ci_y = collect(confint(OneSampleTTest(ys)))
        Plots.plot!(p, [ci_x[1], ci_x[2]], [ci_y[1], ci_y[1]],
            label=nothing,
            seriestype=:line,
            linewidth=1,
            linecolor=color)

        Plots.plot!(p, [ci_x[1], ci_x[2]], [ci_y[2], ci_y[2]],
            label=nothing,
            seriestype=:line,
            linewidth=1,
            linecolor=color)

        Plots.plot!(p, [ci_x[1], ci_x[1]], [ci_y[1], ci_y[2]],
            label=nothing,
            seriestype=:line,
            linewidth=1,
            linecolor=color)

        Plots.plot!(p, [ci_x[2], ci_x[2]], [ci_y[1], ci_y[2]],
            label=nothing,
            seriestype=:line,
            linewidth=1,
            linecolor=color)
    end
end

function help_nonlinear_gradient(lo, mid, hi; grains=100, curve=1.5)
    return vcat(
        [weighted_color_mean((100-i)^curve/grains^curve, lo, mid) for i in 1:grains],
        [weighted_color_mean(1-i^curve/grains^curve, mid, hi) for i in 1:grains]
    )
end

function rotatedLabel(label, x, y)
    angle = atan(y, x) * 180 / pi
    if x < 0
        angle = atan(-y, -x) * 180 / pi
    else
    
    return text(label, :top, default(:xtickfontsize), rotation=angle)
end

# TODO exports

end # module