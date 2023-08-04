# DONE

# use macro helper to define a standard ENA struct with all the bells
@enamodel DigraphENAModel AbstractLinearENAModel

# Documentation
"""
    DigraphENAModel(
        # Required
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1};

        # Optional
        rotation::AbstractLinearENARotation=SVDRotation(),
        unitFilter::Function=unit->true,
        edgeFilter::Function=edge->edge.kind == :directed,
        windowSize::Real=Inf,
        sphereNormalize::Bool=true,
        dropEmpty::Bool=false,
        recenterEmpty::Bool=false
    )

Construct a directed ENA model. Nodes are positioned to maximize goodness of fit between plotted points and units' weighted average of edge vectors.

`DigraphENAModel` follows the same argument and field structure as `ENAModel`.

Ensure that `edgeFilter` only includes `:directed` edges.
"""
DigraphENAModel

# override default model constructor kwargs
function defaultmodelkwargs(
        ::Type{M};
        prev_config::NamedTuple=NamedTuple(),
        kwargs...
    ) where {R<:AbstractLinearENARotation, M<:AbstractDigraphENAModel{R}}

    kwargs = NamedTuple(kwargs)
    super = modelsupertype(M, AbstractDigraphENAModel)
    parentdefaults = defaultmodelkwargs(super)
    definitivedefaults = (
        edgeFilter=(row)->(
            row[:kind] == :directed
        ),
        windowSize=1
    )

    return merge(parentdefaults, prev_config, definitivedefaults, kwargs)
end

function approximate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLinearENARotation, M<:AbstractDigraphENAModel{R}}

    # Regression model for placing the code dots into the approximated high-dimensional space
    ## start with a small amount of noise to prevent colinearity issues
    # X = Matrix{Float64}(rand(nrow(model.accum), nrow(model.nodes)) / 1000000000)
    X = Matrix{Float64}(zeros(nrow(model.accum), nrow(model.nodes)))
    nodeIndexMap = Dict(
        nodeID => i
        for (i, nodeID) in enumerate(model.nodes.nodeID)
    )

#     ## Regression model for placing the code dots into the approximated high-dimensional space
#     X = Matrix{Float64}(zeros(nrow(accumModel), nrow(codeModel)))
#     for (i, unitRow) in enumerate(eachrow(accumModel))
#         for r in relationshipIds
#             a, b = relationshipMap[r]
#             X[i, a] -= unitRow[r]
#             X[i, b] += unitRow[r]
#         end
#     end

    ## for each unit's edges, split the edge between nodes on either end
    for (k, unitRow) in enumerate(eachrow(model.accum))
        for edge in eachrow(model.edges)
            i, j = nodeIndexMap[edge.ground], nodeIndexMap[edge.response]
            X[k, i] -= unitRow[edge.edgeID]
            X[k, j] += unitRow[edge.edgeID]
        end
    end

#     X = (transpose(X) * X)^-1 * transpose(X)
    ## run partial regression
    X = (transpose(X) * X)^-1 * transpose(X)

#     ## Fit each dimension of the original high-dimensional space
#     for networkRow in eachrow(networkModel)
#         r = networkRow[:relationship]
#         y = Vector{Float64}(accumModel[:, r])
#         r_coefs = X * y
#         codeModel[!, r] = r_coefs .- mean(r_coefs) # NOTE in dENA, these MUST be mean centered, else the plot will be illegible
#     end

    ## fit each dimension (edge) of the original high-dimensional space
    for edge in model.edges.edgeID
        y = Vector{Float64}(model.accum[:, edge])
        coefs = X * y # regress on this edge
        model.nodes[!, edge] = coefs[1:end] .- mean(coefs) # NOTE in dENA, these MUST be mean centered, else the plot will be illegible
    end

#     ## Normalize the code positions so they fit on the plot
#     s = maximum(abs.(Matrix{Float64}(codeModel[:, relationshipIds])))
#     if s != 0
#         for codeRow in eachrow(codeModel)
#             codeRow[relationshipIds] = Vector{Float64}(codeRow[relationshipIds]) / s * 0.8
#         end
#     end

# TODO do I need to do the above?

#     ## Refit the units: in high-d space, the refit units are as close as possible to their center of mass wrt the network
#     ## This is by-definition the refit space
#     for networkRow in eachrow(networkModel)
#         r = networkRow[:relationship]
#         for (i, unitRow) in enumerate(eachrow(centroidModel))
#             unitRow[r] =
#                 sum(
#                     codeModel[relationshipMap[k][2], r] * accumModel[i, k] -
#                     codeModel[relationshipMap[k][1], r] * accumModel[i, k]
#                     for k in relationshipIds
#                 )
#         end
#     end

    # find accumHat
    ## Refit the units: in high-d space, the refit units are as close as possible to their
    ## center of mass wrt the network
    ## This is by-definition the refit space
    for unitEdgeID in model.edges.edgeID
        for (k, unitRow) in enumerate(eachrow(model.accumHat))
            unitRow[unitEdgeID] =
                sum(
                    model.nodes[nodeIndexMap[nodeEdge.response], unitEdgeID] * model.accum[k, nodeEdge.edgeID] -
                    model.nodes[nodeIndexMap[nodeEdge.ground  ], unitEdgeID] * model.accum[k, nodeEdge.edgeID]
                    for nodeEdge in eachrow(model.edges)
                )
        end
    end
end

# let the parent handle it from there


# struct DigraphENAModel{T} <: AbstractDigraphENAModel{T}
#     # Inherits:
#     codes::Array{Symbol,1}
#     conversations::Array{Symbol,1}
#     units::Array{Symbol,1}
#     rotation::T
#     accumModel::DataFrame # all the unit-level data we compute
#     centroidModel::DataFrame # accumModel with re-approximated relationship columns
#     metadata::DataFrame
#     codeModel::DataFrame # all the code-level data we compute
#     networkModel::DataFrame # all the connections-level data we compute
#     relationshipMap::Dict{Symbol,Tuple{Int,Int}}

#     # Adds:
#     windowSize::Int
# end

# function DigraphENAModel(
#         # required
#         data::DataFrame, codes::Array{Symbol,1}, conversations::Array{Symbol,1}, units::Array{Symbol,1};
#         # optional
#         windowSize::Int=4, rotateBy::AbstractENARotation=SVDRotation(),
#         sphereNormalize::Bool=true, dropEmpty::Bool=false, recenterEmpty::Bool=false,
#         deflateEmpty::Bool=false, meanCenter::Bool=true, subspaces::Int=0, fitNodesToCircle=false,
#         subsetFilter::Function=x->true, relationshipFilter::Function=(i,j,ci,cj)->(i!=j)
#     )

#     # Checking that the options are sane
#     if windowSize < 2
#         error("The windowSize must be greater than 1")
#     elseif !sphereNormalize && deflateEmpty
#         error("When sphereNormalize=false, deflateEmpty=true has no interpretive validity")
#     elseif dropEmpty && recenterEmpty
#         error("When dropEmpty=true, recenterEmpty=true has no interpretive validity")
#     elseif length(codes) < 3
#         error("There must be at least 3 codes in the model")
#     elseif length(conversations) < 1
#         error("Need at least one column to define conversations")
#     elseif length(units) < 1
#         error("Need at least one column to define units")
#     end

#     # Preparing model structures
#     ## Relationships between codes
#     relationshipMap = Dict{Symbol,Tuple{Int,Int}}()
#     for (i, code1) in enumerate(codes)
#         for (j, code2) in enumerate(codes)
#             if relationshipFilter(i, j, code1, code2)
#                 key = Symbol(string(code1, "_", code2))
#                 relationshipMap[key] = (i, j)
#             end
#         end
#     end

#     relationshipIds = collect(keys(relationshipMap))

#     ## Adding a new column to the raw data that labels each unit properly
#     ena_ids = map(eachrow(data)) do dataRow
#         return join(dataRow[units], ".")
#     end

#     data = hcat(data, DataFrame(:ENA_UNIT => ena_ids))

#     ## Unit model placeholders
#     accumModel = combine(first, groupby(data, :ENA_UNIT))
#     emptyValues = DataFrame(Dict(
#         r => Real[0 for i in 1:nrow(accumModel)]
#         for r in relationshipIds
#     ))

#     emptyPositions = DataFrame(
#         pos_x=Real[0 for i in 1:nrow(accumModel)],
#         pos_y=Real[0 for i in 1:nrow(accumModel)]
#     )

#     accumModel = hcat(accumModel, emptyValues)
#     accumModel = hcat(accumModel, emptyPositions) 
    
#     ## Replace codes "first" metadata with sums instead (ie, code and count)
#     sums = combine(
#         groupby(data[!, [codes..., :ENA_UNIT]], :ENA_UNIT),
#         (code => sum for code in codes)...
#     )
    
#     for code in codes
#         accumModel[!, code] = sums[!, Symbol(string(code, "_sum"))]
#     end

#     ## Network model placeholders
#     networkModel = DataFrame(
#         relationship=relationshipIds,
#         response=[i for (i, j) in values(relationshipMap)],
#         referent=[j for (i, j) in values(relationshipMap)],
#         density=Real[0 for r in relationshipIds], # how thick to make the line
#         weight_x=Real[0 for r in relationshipIds], # the weight I contribute to dim_x's
#         weight_y=Real[0 for r in relationshipIds] # the weight I contribute to dim_y's
#     )
    
#     ## Code model placeholders
#     codeModel = DataFrame(
#         code=codes,
#         density=Real[0 for c in codes], # how thick to make the dot
#         pos_x=Real[0 for c in codes], # where to plot this code on the fitted plot's x-axis
#         pos_y=Real[0 for c in codes] # where to plot this code on the fitted plot's y-axis
#     )

#     emptyValues2 = DataFrame(Dict(
#         r => Real[0 for c in codes] # my position on the approximated x-axis if we use a relationship as the x-axis
#         for r in relationshipIds
#     ))

#     codeModel = hcat(codeModel, emptyValues2)

#     # Accumulation step
#     ## Raw counts for all the units
#     counts = Dict(
#         unit => [[0 for j in codes] for i in codes]
#         for unit in accumModel[!, :ENA_UNIT]
#     )

#     prev_convo = data[1, conversations]
#     howrecents = [Inf for c in codes]
#     for line in eachrow(data)
#         if prev_convo != line[conversations]
#             prev_convo = line[conversations]
#             howrecents .= Inf
#         end

#         for (i, code) in enumerate(codes)
#             if line[code] > 0
#                 howrecents[i] = 0
#             else
#                 howrecents[i] += 1
#             end
#         end

#         unit = line[:ENA_UNIT]
#         for r in relationshipIds
#             i, j = relationshipMap[r]
#             if howrecents[i] == 0 && 0 < howrecents[j] < windowSize
#                 counts[unit][i][j] += 1
#             end
#         end
#     end

#     ## Overwrite the accum model's placeholders
#     for unitRow in eachrow(accumModel)
#         unit = unitRow[:ENA_UNIT]
#         vector = [counts[unit][i][j] for (i,j) in values(relationshipMap)]
#         for (k, r) in enumerate(relationshipIds)
#             unitRow[r] = vector[k]
#         end
#     end

#     # Filtering
#     ## User-defined unit subsetting - we kept all the data in the count model up to this point so the user can define filters as they please
#     filter!(subsetFilter, accumModel)

#     ## Maybe drop rows with empty values or resposition them
#     if dropEmpty
#         ### Find values to keep
#         nonZeroRows = map(eachrow(accumModel)) do unitRow
#             return !all(iszero.(values(unitRow[relationshipIds])))
#         end

#         ### reassign
#         accumModel = accumModel[nonZeroRows, :]
#     elseif recenterEmpty
#         ### Find values that actually have value
#         nonZeroRows = map(eachrow(accumModel)) do unitRow
#             return !all(iszero.(values(unitRow[relationshipIds])))
#         end

#         ### Count them and find the mean point
#         N = sum(nonZeroRows)
#         meanPoint = map(relationshipIds) do r
#             return sum(accumModel[!, r]) / N
#         end

#         ### For all rows, if they are empty, replace with the mean point
#         for unitRow in eachrow(accumModel)
#             if all(iszero.(values(unitRow[relationshipIds])))
#                 for (k, r) in enumerate(relationshipIds)
#                     unitRow[r] = meanPoint[k]
#                 end
#             end
#         end
#     end

#     ## Normalize each unit, if requested
#     if sphereNormalize
#         for unitRow in eachrow(accumModel)
#             vector = Vector{Float64}(unitRow[relationshipIds])
#             s = sqrt(sum(vector .^ 2))
#             if s != 0
#                 unitRow[relationshipIds] = vector / s
#             end
#         end
#     end

#     ## If we didn't sphere normalize, then still scale it down make plots easier to read
#     if !sphereNormalize
#         s = maximum(maximum(accumModel[!, r]) for r in relationshipIds)
#         for r in relationshipIds
#             accumModel[!, r] /= s
#         end
#     end

#     # Now that we have all the data counted, divvy and copy it between accumModel, centroidModel, and metadata
#     modelColumns = [:pos_x, :pos_y, relationshipIds...]
#     nonModelColumns = setdiff(Symbol.(names(accumModel)), modelColumns)
#     metadata = accumModel[!, nonModelColumns]
#     centroidModel = copy(accumModel[!, [:ENA_UNIT, modelColumns...]])
#     accumModel = accumModel[!, [:ENA_UNIT, modelColumns...]]

#     # Compute the position of the codes in the approximated high-dimensional space
#     ## Compute the density of each network row line
#     for networkRow in eachrow(networkModel)
#         r = networkRow[:relationship]
#         networkRow[:density] = sum(accumModel[!, r])
#     end

#     ## Normalize the network densities
#     s = maximum(networkModel[!, :density])
#     networkModel[!, :density] /= s

#     ## Compute the density of each code row dot, by "splitting" the density of each line between its two codes
#     for networkRow in eachrow(networkModel)
#         r = networkRow[:relationship]
#         i, j = relationshipMap[r]
#         codeModel[i, :density] += networkRow[:density]
#         # codeModel[j, :density] += networkRow[:density]
#     end

#     ## Normalize the code densities
#     s = maximum(codeModel[!, :density])
#     codeModel[!, :density] /= s

#     ## Regression model for placing the code dots into the approximated high-dimensional space
#     X = Matrix{Float64}(zeros(nrow(accumModel), nrow(codeModel)))
#     for (i, unitRow) in enumerate(eachrow(accumModel))
#         for r in relationshipIds
#             a, b = relationshipMap[r]
#             X[i, a] -= unitRow[r]
#             X[i, b] += unitRow[r]
#         end
#     end

#     X = (transpose(X) * X)^-1 * transpose(X)

#     ## Fit each dimension of the original high-dimensional space
#     for networkRow in eachrow(networkModel)
#         r = networkRow[:relationship]
#         y = Vector{Float64}(accumModel[:, r])
#         r_coefs = X * y
#         codeModel[!, r] = r_coefs .- mean(r_coefs) # NOTE in dENA, these MUST be mean centered, else the plot will be illegible
#     end

#     ## Normalize the code positions so they fit on the plot
#     s = maximum(abs.(Matrix{Float64}(codeModel[:, relationshipIds])))
#     if s != 0
#         for codeRow in eachrow(codeModel)
#             codeRow[relationshipIds] = Vector{Float64}(codeRow[relationshipIds]) / s * 0.8
#         end
#     end

#     # Rotation step
#     ## Use the given rotation method, probably one of the out-of-the-box ENARotations, but could be anything user defined
#     ## But first, maybe deflate the rotation model
#     subspaceModel = accumModel
#     if deflateEmpty
#         deflatedModel = copy(subspaceModel)
#         zAxis = map(eachrow(networkModel)) do networkRow
#             r = networkRow[:relationship]
#             return sum(subspaceModel[!, r])
#         end

#         s = sqrt(sum(zAxis .^ 2))
#         if s != 0
#             zAxis /= s
#         end

#         deflatedModel[!, :pos_z] = Matrix{Float64}(subspaceModel[!, relationshipIds]) * Vector{Float64}(zAxis)
#         for r in relationshipIds
#             rAxis = map(relationshipIds) do rp
#                 if r == rp
#                     return 1
#                 else
#                     return 0
#                 end
#             end

#             scalar = dot(rAxis, zAxis) / dot(zAxis, zAxis)
#             deflatedModel[!, r] = subspaceModel[!, r] .- (scalar * deflatedModel[!, :pos_z])
#         end

#         subspaceModel = deflatedModel
#     end

#     if subspaces >= 2
#         deflatedModel = copy(subspaceModel)
#         pcaModel = projection(help_deflating_svd(networkModel, subspaceModel))
#         for r in relationshipIds
#             rAxis = map(relationshipIds) do rp
#                 if r == rp
#                     return 1
#                 else
#                     return 0
#                 end
#             end

#             # rAxisNew = rAxis * 0
#             deflatedModel[!, r] = deflatedModel[!, r] * 0
#             for i in 1:min(subspaces, size(pcaModel, 2))
#                 zAxis = pcaModel[:, i]
#                 deflatedModel[!, :pos_z] = Matrix{Float64}(subspaceModel[!, relationshipIds]) * Vector{Float64}(zAxis)

#                 scalar = dot(rAxis, zAxis) / dot(zAxis, zAxis)
#                 deflatedModel[!, r] = deflatedModel[!, r] + scalar * deflatedModel[!, :pos_z]
#             end
#         end

#         subspaceModel = deflatedModel
#     end

#     rotate!(rotateBy, networkModel, codeModel, metadata, subspaceModel)

#     # Layout step
#     ## Project the pos_x and pos_y for the original units onto the plane, now that we have the rotation
#     ## These are what are really drawn (unlike in rENA, where accumModel is what is drawn)
#     accumModel[!, :pos_x] = Matrix{Float64}(accumModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_x])
#     accumModel[!, :pos_y] = Matrix{Float64}(accumModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_y])

#     ## Same for the codes
#     ## These aren't used to compute what's really drawn, they are labels floating around
#     ## in the same high-d space as what is really drawn, that we interpret in terms of the
#     ## center of mass. If we *were* to use these to compute what's really drawn,
#     ## it should actually give us the same result as the projection we used for centroidRow's above,
#     ## since that's the property we defined the refit space to have. (ignoring the intercept)
#     codeModel[!, :pos_x] = Matrix{Float64}(codeModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_x])
#     codeModel[!, :pos_y] = Matrix{Float64}(codeModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_y])

#     # Centroid placement step
#     ## If requested, refit the nodes onto a unit circle. This lowers the goodness of fit, but improves readablity
#     if fitNodesToCircle
#         for codeRow in eachrow(codeModel)
#             vector = Vector{Float64}(codeRow[[:pos_x, :pos_y]])
#             s = sqrt(sum(vector .^ 2))
#             if s != 0
#                 codeRow[[:pos_x, :pos_y]] = vector / s
#                 codeRow[relationshipIds] = Vector{Float64}(codeRow[relationshipIds]) / s
#             end
#         end
#     end

#     ## Refit the units: in high-d space, the refit units are as close as possible to their center of mass wrt the network
#     ## This is by-definition the refit space
#     for networkRow in eachrow(networkModel)
#         r = networkRow[:relationship]
#         for (i, unitRow) in enumerate(eachrow(centroidModel))
#             unitRow[r] =
#                 sum(
#                     codeModel[relationshipMap[k][2], r] * accumModel[i, k] -
#                     codeModel[relationshipMap[k][1], r] * accumModel[i, k]
#                     for k in relationshipIds
#                 )
#         end
#     end

#     ## Same for the refit units
#     ## This is really only for testing the goodness of fit
#     centroidModel[!, :pos_x] = Matrix{Float64}(centroidModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_x])
#     centroidModel[!, :pos_y] = Matrix{Float64}(centroidModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_y])

#     ## Maybe translate everything so overall mean lies at the origin
#     if meanCenter
#         mu_x = mean(centroidModel[!, :pos_x])
#         mu_y = mean(centroidModel[!, :pos_y])
#         centroidModel[!, :pos_x] = centroidModel[!, :pos_x] .- mu_x
#         centroidModel[!, :pos_y] = centroidModel[!, :pos_y] .- mu_y
#         mu_x = mean(accumModel[!, :pos_x])
#         mu_y = mean(accumModel[!, :pos_y])
#         accumModel[!, :pos_x] = accumModel[!, :pos_x] .- mu_x
#         accumModel[!, :pos_y] = accumModel[!, :pos_y] .- mu_y
#     end

#     # Testing step
#     ## Test that the angle between the dimensions is 90 degrees
#     theta = dot(networkModel[!, :weight_x], networkModel[!, :weight_y])
#     theta /= sqrt(dot(networkModel[!, :weight_x], networkModel[!, :weight_x]))
#     theta /= sqrt(dot(networkModel[!, :weight_y], networkModel[!, :weight_y]))
#     angle = acos(theta) * 180 / pi
#     if abs(angle-90) > 0.0001 # allow for a little approximation error
#         @warn """The angle between the axes of this model is $(angle) degrees, when it should be 90.
# This can lead to strange visual effects when plotting on orthogonal axes.
# This can undermine interpreting betweenness between units.
# This can undermind interpreting variance explained by the axes.
# And this can cause problems with ENA's optimization algorithm fitting the codes and the lines."""
#     end

#     # Done!
#     return DigraphENAModel(
#         codes, conversations, units, rotateBy,
#         accumModel, centroidModel, metadata, codeModel, networkModel,
#         relationshipMap,
#         windowSize
#     )
# end

# # Override plotting pieces
# ## Base - Inject a groupBy and some labels when none are given
# function plot(ena::AbstractDigraphENAModel{T};
#     showArrows=nothing, reverseLineSort=nothing,
#     kwargs...) where T <: AbstractMeansRotation

#     if isnothing(showArrows)
#         showArrows = true
#     else
#         showArrows = false
#     end

#     if isnothing(reverseLineSort)
#         reverseLineSort = true
#     else
#         reverseLineSort = false
#     end

#     return invoke(plot, Tuple{AbstractENAModel{<:T}}, ena;
#                   showArrows=showArrows, reverseLineSort=reverseLineSort, kwargs...)
# end