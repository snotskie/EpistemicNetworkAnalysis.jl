# Helpers
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

# Implementation
function defaultmodelkwargs(
        ::Type{M};
        prev_config::NamedTuple=NamedTuple(),
        unitFilter::Function=x->true,
        edgeFilter::Function=x->true,
        windowSize::Real=Inf,
        sphereNormalize::Bool=true,
        dropEmpty::Bool=false,
        recenterEmpty::Bool=false,
        # deflateEmpty::Bool=false,
        kwargs...
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    kwargs = NamedTuple(kwargs)
    defaults = (
        unitFilter=unitFilter,
        edgeFilter=edgeFilter,
        windowSize=windowSize,
        sphereNormalize=sphereNormalize,
        dropEmpty=dropEmpty,
        recenterEmpty=recenterEmpty,
        # deflateEmpty=deflateEmpty,
        kwargs...
    )

    return merge(defaults, prev_config, kwargs)
end

function populateENAfields(
        ::Type{M},
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1},
        rotation::AbstractLinearENARotation;
        config...
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    # sanity checks
    if nrow(data) < 3
        error("The data::DataFrame parameter should have at least 3 rows")
    elseif length(codes) < 3
        error("The codes::Array parameter should have at least 3 items")
    elseif length(units) == 0
        error("The units::Array parameter should have at least 1 item") 
    end

    # note, empty conversations is allowed. amounts to a whole conversation model

    config = NamedTuple(config)

    # edges: empty starting point
    edges = DataFrame(
        edgeID=Symbol[],
        kind=Symbol[],
        ground=Symbol[],
        response=Symbol[]
    )

    # edges: fill in all possible combinations
    for (i, ground) in enumerate(codes)
        for (j, response) in enumerate(codes)
            if i == j # single code count edges
                ckey = Symbol(string("_", ground, "_"))
                ekey = Symbol(string(ground, "_to_", response))
                newEdge = similar(edges, 2)
                newEdge[1, :edgeID] = ckey
                newEdge[1, :kind] = :count
                newEdge[1, :ground] = ground
                newEdge[1, :response] = response
                newEdge[2, :edgeID] = ekey
                newEdge[2, :kind] = :echo
                newEdge[2, :ground] = ground
                newEdge[2, :response] = response
                append!(edges, newEdge)
            elseif i < j # directed and undirected edges
                ukey = Symbol(string(ground, "_", response))
                dkey = Symbol(string(ground, "_to_", response))
                newEdges = similar(edges, 2)
                newEdges[1, :edgeID] = ukey
                newEdges[1, :kind] = :undirected
                newEdges[1, :ground] = ground
                newEdges[1, :response] = response
                newEdges[2, :edgeID] = dkey
                newEdges[2, :kind] = :directed
                newEdges[2, :ground] = ground
                newEdges[2, :response] = response
                append!(edges, newEdges)
            else # remaining directed edges
                dkey = Symbol(string(ground, "_to_", response))
                newEdges = similar(edges, 1)
                newEdges[1, :edgeID] = dkey
                newEdges[1, :kind] = :directed
                newEdges[1, :ground] = ground
                newEdges[1, :response] = response
                append!(edges, newEdges)
            end
        end
    end

    # edges: filter out unused combinations
    filter!(config.edgeFilter, edges)
    edgeIDs = edges.edgeID

    # nodes: zero'd starting point
    nodes = DataFrame(Dict(
        :nodeID=>codes,
        (name=>zeros(Real, length(codes)) for name in edgeIDs)...
    ))

    # accum, accumHat, metadata: generate IDs for each unit
    unitIDs = map(eachrow(data)) do dataRow
        return Symbol(join(dataRow[units], "."))
    end

    # accum, accumHat, metadata: add to the data and group by unit ID
    if :unitID in Symbol.(names(data))
        data.unitIDs .= unitIDs
    else
        data = hcat(data, DataFrame(:unitID => unitIDs))
    end

    tempAccum = combine(first, groupby(data, :unitID))

    # accum, accumHat, metadata: filter unused units
    filter!(config.unitFilter, tempAccum)

    # accum, accumHat, metadata: placeholder zeros for all unit/edge pairs
    tempValues = DataFrame(Dict(
        edgeID=>zeros(Real, nrow(tempAccum))
        for edgeID in edgeIDs
    ))

    # accum, accumHat, metadata: add placeholders and split into multiple dfs
    # makeunique=true because 
    tempAccum = hcat(tempAccum, tempValues)
    accum = copy(tempAccum[!, [:unitID, edgeIDs...]])
    accumHat = copy(accum)
    metaNames = setdiff(Symbol.(names(tempAccum)), edgeIDs)
    metadata = copy(tempAccum[!, metaNames])

    # embedding: empty starting point
    embedding = DataFrame(Dict(
        :label=>String[],
        :variance_explained=>Real[],
        # :pearson=>Real[],
        :coregistration=>Real[],
        (edgeID=>Real[] for edgeID in edgeIDs)...
    ))

    # points: empty starting point
    points = DataFrame(Dict(
        unitID=>Real[] for unitID in accum.unitID
    ))

    pointsHat = similar(points, 0)
    pointsNodes = DataFrame(Dict(
        nodeID=>Real[] for nodeID in nodes.nodeID
    ))

    # done!
    return data,
        codes,
        conversations,
        units,
        rotation,
        metadata,
        points,
        pointsHat,
        pointsNodes,
        accum,
        accumHat,
        edges,
        nodes,
        embedding,
        config
end

function accumulate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    prev_convo = model.data[1, model.conversations]
    howrecents = Dict(
        node.nodeID => Inf
        for node in eachrow(model.nodes)
    )

    unitIndexMap = Dict(
        unit.unitID => i
        for (i, unit) in enumerate(eachrow(model.accum))
    )

    for line in eachrow(model.data)

        # reset on new conversations
        if prev_convo != line[model.conversations]
            prev_convo = line[model.conversations]
            for key in keys(howrecents)
                howrecents[key] = Inf
            end
        end

        # count how recently codes occured: increment everyone
        for nodeID in model.nodes.nodeID
            howrecents[nodeID] += 1
        end

        # update the accum counts of the unit on this line
        if haskey(unitIndexMap, line.unitID) # in case the unit was filtered out above
            unit = unitIndexMap[line.unitID]
            for edge in eachrow(model.edges)
                if edge.kind == :count
                    if line[edge.response] > 0
                        model.accum[unit, edge.edgeID] += 1
                    end
                elseif edge.kind == :echo || edge.kind == :directed
                    if (line[edge.response] > 0 && howrecents[edge.ground] < model.config.windowSize) ||
                        (line[edge.ground] > 0 && line[edge.response] > 0)
                        model.accum[unit, edge.edgeID] += 1
                    end
                elseif edge.kind == :undirected
                    if (line[edge.response] > 0 && howrecents[edge.ground] < model.config.windowSize) ||
                    (line[edge.ground] > 0 && howrecents[edge.response] < model.config.windowSize) ||
                    (line[edge.ground] > 0 && line[edge.response] > 0)
                        model.accum[unit, edge.edgeID] += 1
                    end
                end
            end
        end

        # count how recently codes occured: reset the found ones
        for nodeID in model.nodes.nodeID
            if line[nodeID] > 0
                howrecents[nodeID] = 0
            end
        end
    end

    # normalize each unit, if requested
    edgeIDs = model.edges.edgeID
    if model.config.sphereNormalize
        for i in 1:nrow(model.accum)
            vector = Vector{Float64}(model.accum[i, edgeIDs])
            s = sqrt(sum(vector .^ 2))
            if s != 0
                model.accum[i, edgeIDs] = vector / s
            end
        end
    # else, still scale it down to make plots easier to read
    else
        s = maximum(maximum(model.accum[!, r]) for r in edgeIDs)
        for r in edgeIDs
            model.accum[!, r] /= s
        end
    end

    # maybe drop rows with empty values
    if model.config.dropEmpty
        droppedRows = map(eachrow(model.accum)) do unitRow
            return all(iszero.(values(unitRow[edgeIDs])))
        end

        droppedUnits = model.accum[droppedRows, :]
        model.accum = model.accum[.!droppedRows, :]
        model.accumHat = model.accumHat[.!droppedRows, :]
        model.metadata = model.metadata[.!droppedRows, :]
        select!(model.points, Not(droppedUnits.unitID))
        select!(model.pointsHat, Not(droppedUnits.unitID))
    # else, maybe recenter the empty rows to the mean
    elseif model.config.recenterEmpty
        zeroRows = map(eachrow(model.accum)) do unitRow
            return all(iszero.(values(unitRow[edgeIDs])))
        end

        nonZeroRows = .!zeroRows
        N = sum(nonZeroRows)
        meanPoint = transpose(sum.(eachcol(model.accum[nonZeroRows, edgeIDs])) / N)
        model.accum[zeroRows, edgeIDs] .= meanPoint
    # else, maybe deflate the model so empty and the mean always align
    # elseif model.config.deflateEmpty
    #     meanAxis = sum.(eachcol(model.accum[!, edgeIDs]))
    #     s = sqrt(sum(meanAxis .^ 2))
    #     if s != 0
    #         meanAxis /= s
    #     end

    #     meanPoints = Matrix{Float64}(model.accum[!, edgeIDs]) * meanAxis
    #     for edge in edgeIDs
    #         edgeAxis = [edge == edgep ? 1 : 0 for edgep in edgeIDs]
    #         scalar = dot(edgeAxis, meanAxis) / dot(meanAxis, meanAxis)
    #         model.accum[!, edge] .-= scalar * meanPoints
    #     end
    end

    # renormalize if some of the above happened
    # if model.config.deflateEmpty
    #     if model.config.sphereNormalize
    #         for i in 1:nrow(model.accum)
    #             vector = Vector{Float64}(model.accum[i, edgeIDs])
    #             s = sqrt(sum(vector .^ 2))
    #             if s != 0
    #                 model.accum[i, edgeIDs] = vector / s
    #             end
    #         end
    #     else
    #         s = maximum(maximum(model.accum[!, r]) for r in edgeIDs)
    #         for r in edgeIDs
    #             model.accum[!, r] /= s
    #         end
    #     end
    # end
end

function substantiate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    # Regression model for placing the code dots into the approximated high-dimensional space
    ## start with a small amount of noise to prevent colinearity issues
    # X = Matrix{Float64}(rand(nrow(model.accum), nrow(model.nodes)) / 1000000000)
    X = Matrix{Float64}(zeros(nrow(model.accum), nrow(model.nodes)))
    nodeIndexMap = Dict(
        nodeID => i
        for (i, nodeID) in enumerate(model.nodes.nodeID)
    )

    ## for each unit's edges, split the edge between nodes on either end
    for (k, unitRow) in enumerate(eachrow(model.accum))
        for edge in eachrow(model.edges)
            i, j = nodeIndexMap[edge.ground], nodeIndexMap[edge.response]
            X[k, i] += unitRow[edge.edgeID] / 2
            X[k, j] += unitRow[edge.edgeID] / 2
        end
    end

    ## run partial regression
    X = (transpose(X) * X)^-1 * transpose(X)

    ## fit each dimension (edge) of the original high-dimensional space
    for edge in model.edges.edgeID
        y = Vector{Float64}(model.accum[:, edge])
        coefs = X * y # regress on this edge
        model.nodes[!, edge] = coefs[1:end]
    end
end

function approximate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    # find accumHat
    ## Refit the units: in high-d space, the refit units are as close as possible to their
    ## center of mass wrt the network
    ## This is by-definition the refit space
    nodeIndexMap = Dict(
        nodeID => i
        for (i, nodeID) in enumerate(model.nodes.nodeID)
    )
    
    for unitEdgeID in model.edges.edgeID
        for (k, unitRow) in enumerate(eachrow(model.accumHat))
            unitRow[unitEdgeID] =
                sum(
                    model.nodes[nodeIndexMap[nodeEdge.ground  ], unitEdgeID] * model.accum[k, nodeEdge.edgeID] +
                    model.nodes[nodeIndexMap[nodeEdge.response], unitEdgeID] * model.accum[k, nodeEdge.edgeID]
                    for nodeEdge in eachrow(model.edges)
                ) / 2
        end
    end
end

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    # for existing dimensions given to us by the child class,
    # normalize them, orthogonalize them from each other,
    # and add point positions to the model
    edgeIDs = model.edges.edgeID
    numExistingDims = nrow(model.embedding)
    if numExistingDims > 0
        # normalize the first dimension (the rest are normalized below)
        v1 = Vector(model.embedding[1, edgeIDs])
        s = sqrt(sum(v1 .^ 2))
        if s != 0
            model.embedding[1, edgeIDs] .= v1 / s
        end
    end

    for i in 1:numExistingDims
        addPointsToModelFromDim(model, i)

        # orthogonalize, by rejection, each existing dimension from each previous dimension.
        # note, a dimension will be orthogonalize before it is used to add points to the model
        vi = Vector(model.embedding[i, edgeIDs])
        denom = dot(vi, vi)
        for j in (i+1):numExistingDims
            vj = Vector(model.embedding[j, edgeIDs])
            s = sqrt(sum(vj .^ 2))
            if s != 0
                model.embedding[j, edgeIDs] .= vj / s
            end

            scalar = dot(vj, vi) / denom
            model.embedding[j, edgeIDs] .= vj - scalar * vi
            s = sqrt(sum(vj .^ 2))
            if s < 0.05
                model.embedding[j, edgeIDs] .= 0
                @warn "During the rotation step, axis $j was deflated to zero due to close correlation with axis $i."
            elseif s != 0
                model.embedding[j, edgeIDs] .= vj / s
            end
        end
    end

    # make a mean centered copy of accum
    edgeIDs = model.edges.edgeID
    X = Matrix{Float64}(model.accum[!, edgeIDs])
    for i in axes(X, 2)
        X[:, i] .-= mean(X[:, i])
    end

    # if needed, deflate those points' dimensions from X before we run SVD on it
    unitIDs = model.accum.unitID
    for i in 1:numExistingDims
        col = Vector{Float64}(model.points[i, unitIDs])
        col .-= mean(col)
        denom = dot(col, col)
        for j in axes(X, 2)
            scalar = dot(X[:, j], col) / denom
            X[:, j] -= scalar * col
        end
    end

    # then, once we've deflated or not, we run SVD on the data, then add to the model
    svd = transpose(projection(fit(PCA, X', pratio=1.0)))
    df = similar(model.embedding, size(svd, 1))
    df[!, edgeIDs] = svd
    df.label = ["SVD$(i)" for i in 1:nrow(df)]
    append!(model.embedding, df)
    for i in (numExistingDims+1):nrow(model.embedding)
        addPointsToModelFromDim(model, i)
    end
end

# function tests(
#         ::Type{M}, model::AbstractLinearENAModel
#     ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

#     # TODO maybe cut this, think about tests as summary stats dict
# end