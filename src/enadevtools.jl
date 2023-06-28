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

        # make default constructor
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

        # make re-modeler constructor
        function $(esc(self))(
                prev_model::AbstractENAModel;
                rotateBy::R=prev_model.rotation, #$(esc(defaultrotation))(),
                kwargs...
            ) where {R<:$(esc(rotationtype))}

            # call common ENA re-modeler
            return remodelENA(
                $(esc(self)){R},
                prev_model,
                rotateBy;
                kwargs...
            )
        end

        # make re-rotator constructor
        function $(esc(self))(
                prev_model::$(esc(self));
                rotateBy::R=$(esc(defaultrotation))()
            ) where {R<:$(esc(rotationtype))}

            # call common ENA re-rerotator
            return rerotateENA(
                $(esc(self)){R},
                prev_model,
                rotateBy
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
    test!(M, model)
    return model
end

# Re-model constructor
function remodelENA(
        ::Type{M},
        prev_model::AbstractENAModel,
        rotation::AbstractENARotation;
        kwargs...
    ) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}

    kwargs = defaultmodelkwargs(M; prev_config=prev_model.config, kwargs...)
    model = M(
        populateENAfields(
            M,
            copy(prev_model.data),
            copy(prev_model.codes),
            copy(prev_model.conversations),
            copy(prev_model.units),
            rotation;
            kwargs...
        )...
    )

    accumulate!(M, model)
    approximate!(M, model)
    rotate!(M, model)
    test!(M, model)
    return model
end

# Re-rotate constructor
function rerotateENA(
        ::Type{M},
        prev_model::AbstractENAModel,
        rotation::AbstractENARotation
    ) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}

    model = M(
        copy(prev_model.data),
        copy(prev_model.codes),
        copy(prev_model.conversations),
        copy(prev_model.units),
        rotation,
        copy(prev_model.metadata),
        similar(prev_model.points, 0),
        similar(prev_model.pointsHat, 0),
        similar(prev_model.pointsNodes, 0),
        copy(prev_model.accum),
        copy(prev_model.accumHat),
        copy(prev_model.edges),
        copy(prev_model.nodes),
        similar(prev_model.embedding, 0),
        prev_model.config
    )

    rotate!(M, model)
    test!(M, model)
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

function defaultplotkwargs(::Type{M}, model::AbstractENAModel, config::NamedTuple) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    # wrapper for the above, to simplify logic of children that override their parent plotkwargs
    defaultplotkwargs(M, model; Dict(zip(keys(config), values(config)))...)
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

function test!(::Type{M}, model::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function test!(::Type{M}, model::AbstractENAModel, test::Type{<:HypothesisTests.HypothesisTest}; kwargs...) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function summary(::Type{M}, model::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    columns = [
        :label,
        setdiff(Symbol.(names(model.embedding)), [:label, model.edges.edgeID...])...
    ]

    return copy(model.embedding[!, columns])
end

function show(io::IO, model::AbstractENAModel)
    details = (
        ENATool="EpistemicNetworkAnalysis.jl",
        ToolVersion=PkgVersion.Version(EpistemicNetworkAnalysis),
        ToolAuthor=PkgVersion.Author(EpistemicNetworkAnalysis),
        ModelConfig=(
            codes=model.codes,
            conversations=model.conversations,
            units=model.units,
            model.config...
        ),
        NumberOfUnits=nrow(model.metadata),
        RotationType=string(nameof(typeof(model.rotation))),
        RotationConfig=
            length(propertynames(model.rotation)) > 0 ?
            namedtuple(propertynames(model.rotation), fieldvalues(model.rotation)) :
            (),
        Dimensions=Tables.rowtable(summary(model))
    )

    pprint(io, details)
end

# in linear, do plot like the consruct helper, override its components under there
function plot(::Type{M}, model::AbstractENAModel, plotconfig::NamedTuple) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

## Wrapper Functions
function summary(model::AbstractENAModel)
    return summary(typeof(model), model)
end

# function show(model::AbstractENAModel)
#     return show(typeof(model), model)
# end

function plot(model::AbstractENAModel; kwargs...)
    plotconfig = defaultplotkwargs(typeof(model), model; kwargs...)
    return plot(typeof(model), model, plotconfig)
end

# Serialization
function to_xlsx(filename::AbstractString, model::AbstractENAModel)
    i = 1
    namesheet(xf, sheetname) = begin
        if i == 1
            XLSX.rename!(xf[1], sheetname)
        else
            XLSX.addsheet!(xf, sheetname)
        end
    end

    fixcelltype(x) = begin
        if typeof(x) <: Union{Missing, Bool, Float64, Int64, Dates.Date, Dates.DateTime, Dates.Time, String}
            return x
        elseif typeof(x) <: Function
            return join(string.(code_lowered(x)[1].code), "; ")
        else
            return string(x)
        end
    end

    fixcelltypes(col) = fixcelltype.(col)
    prepare(df) = fixcelltypes.(eachcol(df))
    putsheet(xf, df, sheetname) = begin
        namesheet(xf, sheetname)
        XLSX.writetable!(xf[i], prepare(df), names(df))
        i += 1
    end

    putvector(xf, pair, sheetname) = begin
        namesheet(xf, sheetname)
        XLSX.writetable!(xf[i], [fixcelltypes(last(pair))], [first(pair)])
        i += 1
    end

    puttuple(xf, tp, sheetname) = begin
        namesheet(xf, sheetname)
        XLSX.writetable!(xf[i], fixcelltypes.([[x] for x in values(tp)]), keys(tp))
        i += 1
    end

    putrotation(xf, rotation, sheetname) = begin
        namesheet(xf, sheetname)
        cols = ["RotationType", propertynames(rotation)...]
        vals = [string(nameof(typeof(rotation))), fixcelltypes.(fieldvalues(rotation))...]
        XLSX.writetable!(xf[i], [[x] for x in vals], cols)
        i += 1
    end

    XLSX.openxlsx(filename, mode="w") do xf
        putsheet(xf, DataFrame([
            "Data"          "The original data given to the ENA model";
            "Codes"         "A list of the column names defining the qualitative codes used in the model";
            "Conversations" "A list of the column names used to distinguish conversations";
            "Units"         "A list of the column names used to distinguish units of analysis";
            "Model"         "A row recording the configuration of the model used";
            "Rotation"      "A row recording the configuration of the rotation used";
            "Metadata"      "The metadata associated with each unit of analysis";
            "Accum"         "The normalized accumulated connection weights for each unit of analysis"
            "Points"        "The plotted unit points, one column per unit of analysis, one row per dimension. By default, the first row is the X-axis and the second row is the Y-axis";
            "Edges"         "The edges used in the model";
            "Embedding"     "The embedding arrived at by the model, one column per edge, one row per dimension. Rows are sorted in parallel to the Points sheet";
            "Nodes"         "The pseudo-accumulated weights for each node, found by fitting nodes to each column of Accum, and used to calculate PointsNodes";
            "PointsNodes"   "The plotted node points, one column per node, one row per dimension, found by Embedding * Nodes. Rows are sorted in parallel to the Points sheet";
            "AccumHat"      "The approximated unit weights, used to calculate PointsHat";
            "PointsHat"     "The approximated plotted unit points (aka 'centroids'), found by Embedding * AccumHat, used to calculate model goodness of fit"
        ], ["Sheet", "Description"]), "Documentation")
        putsheet(xf, model.data, "Data")
        putvector(xf, "Codes" => model.codes, "Codes")
        putvector(xf, "Conversations" => model.conversations, "Conversations")
        putvector(xf, "Units" => model.units, "Units")
        puttuple(xf, (
            ModelType=string(nameof(typeof(model))),
            model.config...
        ), "Model")
        putrotation(xf, model.rotation, "Rotation")
        putsheet(xf, model.metadata, "Metadata")
        putsheet(xf, model.accum, "Accum")
        putsheet(xf, model.points, "Points")
        putsheet(xf, model.edges, "Edges")
        putsheet(xf, model.embedding, "Embedding")
        putsheet(xf, model.nodes, "Nodes")
        putsheet(xf, model.pointsNodes, "PointsNodes")
        putsheet(xf, model.accumHat, "AccumHat")
        putsheet(xf, model.pointsHat, "PointsHat")
    end
end

# function from_xlsx(filename::AbstractString)
#     xf = XLSX.readxlsx(filename)
#     sheetnames = XLSX.sheetnames(xf)
#     expected = [
#         "Data",
#         "Codes",
#         "Conversations",
#         "Units",
#         "Model",
#         "Rotation",
#         "Metadata",
#         "Accum",
#         "Points",
#         "Edges",
#         "Embedding",
#         "Nodes",
#         "PointsNodes",
#         "AccumHat",
#         "PointsHat"
#     ]

#     for sheetname in expected
#         @assert(sheetname in sheetnames,
#             "Expected sheet named '$sheetname' not found in the spreadsheet: $filename")
#     end

#     # Helper
#     df(sheetname) = DataFrame(XLSX.gettable(xf[sheetname]))

#     # Basic config
#     data = df("Data")
#     conversations = df("Conversations").Conversations
#     units = df("Units").Units
#     codes = df("Codes").Codes

#     # Rotation
#     rotationTypeName = df("Rotation").RotationType[1]
#     rotationType = try
#         getfield(@__MODULE__, Symbol(rotationTypeName))
#     catch UndefVarError
#         @error "Unknown rotation type: $rotationTypeName"
#     end

#     @assert(rotationType <: AbstractENARotation,
#         "Invalid rotation type: $rotationTypeName")
    
#     rotationArgs = df("Rotation")[1, 2:end]
#     rotation = rotationType(rotationArgs...)
#     # This is the tricky part, zero-knowledge type conversion
#     # from whatever type XLSX and DataFrame comes up with,
#     # to whatever type (or array type) the rotation requires
    
#     # Module components
#     # Should be mostly just loading in dataframes

#     # Module config
#     # Same trickiness as the rotation

#     # Model re-construction
#     # The idea is to just call the underlying struct directly

#     return xf # TODO return the reconstructed model
# end

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

function addPointsToModelFromDim(model, i)
    edgeIDs = model.edges.edgeID
    unitIDs = model.accum.unitID
    nodeIDs = model.nodes.nodeID
    axis = Vector{Float64}(model.embedding[i, edgeIDs])
    points = Matrix{Float64}(model.accum[!, edgeIDs]) * axis
    df = similar(model.points, 1)
    df[1, unitIDs] = points
    append!(model.points, df)
    pointsHat = Matrix{Float64}(model.accumHat[!, edgeIDs]) * axis
    df = similar(model.pointsHat, 1)
    df[1, unitIDs] = pointsHat
    append!(model.pointsHat, df)
    pointsNodes = Matrix{Float64}(model.nodes[!, edgeIDs]) * axis
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

    # BUGFIX: divide the default arrow size by ceil(sqrt(number of plots in the multi plot)) = number of plots on the first row
    GR.setarrowsize(1/N)

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