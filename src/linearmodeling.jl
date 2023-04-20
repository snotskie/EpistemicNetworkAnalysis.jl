function defaultmodelkwargs(
        ::Type{M};
        unitFilter::Function=x->true,
        edgeFilter::Function=x->true,
        windowSize::Real=Inf,
        sphereNormalize::Bool=true,
        dropEmpty::Bool=false,
        recenterEmpty::Bool=false,
        deflateEmpty::Bool=false,
        rotateBy::AbstractLinearENARotation=SVDRotation(),
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
        deflateEmpty=deflateEmpty,
        rotateBy=rotateBy,
        kwargs...
    )

    return merge(defaults, kwargs)
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

    config = NamedTuple(config)

    # edges: empty starting point
    edges = DataFrame(
        edge=Symbol[],
        kind=Symbol[],
        i=Int[],
        j=Int[],
        ground=Symbol[],
        response=Symbol[]
    )

    # edges: fill in all possible combinations
    for (i, ground) in enumerate(codes)
        for (j, response) in enumerate(codes)
            if i == j # single code count edges
                newEdge = similar(edges, 1)
                newEdge[1, :edge] = Symbol(string("_", ground, "_"))
                newEdge[1, :kind] = :count
                newEdge[1, :i] = i
                newEdge[1, :j] = j
                newEdge[1, :ground] = ground
                newEdge[1, :response] = response
                append!(edges, newEdge)
            elseif i < j # directed and undirected edges
                ukey = Symbol(string(ground, "_", response))
                dkey = Symbol(string(ground, "_to_", response))
                newEdges = similar(edges, 2)
                newEdges[1, :edge] = ukey
                newEdges[1, :kind] = :undirected
                newEdges[1, :i] = i
                newEdges[1, :j] = j
                newEdges[1, :ground] = ground
                newEdges[1, :response] = response
                newEdges[2, :edge] = dkey
                newEdges[2, :kind] = :directed
                newEdges[2, :i] = i
                newEdges[2, :j] = j
                newEdges[2, :ground] = ground
                newEdges[2, :response] = response
                append!(edges, newEdges)
            else # remaining directed edges
                dkey = Symbol(string(ground, "_to_", response))
                newEdges = similar(edges, 1)
                newEdges[1, :edge] = dkey
                newEdges[1, :kind] = :directed
                newEdges[1, :i] = i
                newEdges[1, :j] = j
                newEdges[1, :ground] = ground
                newEdges[1, :response] = response
                append!(edges, newEdges)
            end
        end
    end

    # edges: filter out unused combinations
    filter!(config.edgeFilter, edges)
    edgeNames = edges.edge

    # nodes: zero'd starting point
    nodes = DataFrame(Dict(
        :node=>codes,
        (name=>zeros(Real, length(codes)) for name in edgeNames)...
    ))

    # accum, accumHat, metadata: generate IDs for each unit
    unitIDs = map(eachrow(data)) do dataRow
        return Symbol(join(dataRow[units], "."))
    end

    # accum, accumHat, metadata: add to the data and group by unit ID
    data = hcat(data, DataFrame(:unitID => unitIDs))
    tempAccum = combine(first, groupby(data, :unitID))

    # accum, accumHat, metadata: filter unused units
    filter!(config.unitFilter, tempAccum)

    # accum, accumHat, metadata: placeholder zeros for all unit/edge pairs
    tempValues = DataFrame(Dict(
        edge=>zeros(Real, nrow(tempAccum))
        for edge in edgeNames
    ))

    # accum, accumHat, metadata: add placeholders and split into multiple dfs
    # makeunique=true because 
    tempAccum = hcat(tempAccum, tempValues)
    accum = copy(tempAccum[!, [:unitID, edgeNames...]])
    accumHat = copy(accum)
    metaNames = setdiff(Symbol.(names(tempAccum)), edgeNames)
    metadata = copy(tempAccum[!, metaNames])

    # embedding: empty starting point
    embedding = DataFrame(Dict(
        :label=>String[],
        :variance_explained=>Real[],
        :pearson=>Real[],
        :coregistration=>Real[],
        (name=>Real[] for name in edgeNames)...
    ))

    # points: empty starting point
    points = DataFrame(Dict(
        unit=>Real[] for unit in accum.unitID
    ))

    pointsHat = similar(points, 0)
    pointsNodes = DataFrame(Dict(
        node=>Real[] for node in nodes.node
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
    howrecents = [Inf for c in model.codes]
    unitIndexMap = Dict(
        unitRow.unitID => i
        for (i, unitRow) in enumerate(eachrow(model.accum))
    )

    for line in eachrow(model.data)

        # reset on new conversations
        if prev_convo != line[model.conversations]
            prev_convo = line[model.conversations]
            howrecents .= Inf
        end

        # count how recently codes occured
        for (i, code) in enumerate(model.codes)
            if line[code] > 0
                howrecents[i] = 0
            else
                howrecents[i] += 1
            end
        end

        # update the accum counts of the unit on this line
        unit = unitIndexMap[line.unitID]
        for edge in eachrow(model.edges)
            if edge.kind == :count
                i = edge.i
                if howrecents[i] == 0
                    model.accum[unit, edge.edge] += 1
                end
            elseif edge.kind == :directed
                i, j = edge.i, edge.j
                if howrecents[j] == 0 && howrecents[i] < model.config.windowSize
                    model.accum[unit, edge.edge] += 1
                end
            elseif edge.kind == :undirected               
                i, j = edge.i, edge.j
                if howrecents[i] == 0 && howrecents[j] < model.config.windowSize
                    model.accum[unit, edge.edge] += 1
                elseif howrecents[j] == 0 && howrecents[i] < model.config.windowSize
                    model.accum[unit, edge.edge] += 1
                end
            end
        end
    end

    # maybe drop rows with empty values
    edgeNames = model.edges.edge
    if model.config.dropEmpty
        droppedRows = map(eachrow(model.accum)) do unitRow
            return all(iszero.(values(unitRow[edgeNames])))
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
            return all(iszero.(values(unitRow[edgeNames])))
        end

        nonZeroRows = .!zeroRows
        N = sum(nonZeroRows)
        meanPoint = transpose(sum(eachcol(model.accum[nonZeroRows, edgeNames])) / N)
        model.accum[zeroRows, edgeNames] .= meanPoint
    # else, maybe deflate the model so empty and the mean always align
    elseif model.config.deflateEmpty
        meanAxis = transpose(sum(eachcol(model.accum[!, edgeNames])))
        s = sqrt(sum(meanAxis .^ 2))
        if s != 0
            meanAxis /= s
        end

        meanPoints = Matrix{Float64}(model.accum[!, edgeNames]) * Vector{Float64}(meanAxis)
        for edge in edgeNames
            edgeAxis = [edge == edgep ? 1 : 0 for edgep in edgeNames]
            scalar = dot(edgeAxis, meanAxis) / dot(meanAxis, meanAxis)
            model.accum[!, edge] .-= scalar * meanPoints
        end
    end

    # normalize each unit, if requested
    if model.config.sphereNormalize
        for i in 1:nrow(model.accum)
            vector = Vector{Float64}(model.accum[i, edgeNames])
            s = sqrt(sum(vector .^ 2))
            if s != 0
                model.accum[i, edgeNames] = vector / s
            end
        end
    # else, still scale it down to make plots easier to read
    else
        s = maximum(maximum(model.accum[!, r]) for r in edgeNames)
        for r in edgeNames
            model.accum[!, r] /= s
        end
    end
end

function approximate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    # compute densities
    edgeDensityDict, nodeDensityDict = computeNetworkDensities(model)

    # Regression model for placing the code dots into the approximated high-dimensional space
    ## start with a small amount of noise to prevent colinearity issues
    X = Matrix{Float64}(rand(nrow(model.accum), nrow(model.nodes)) / 1000000000)

    ## for each unit's edges, split the edge between nodes on either end
    for (k, unitRow) in enumerate(eachrow(model.accum))
        for edge in eachrow(model.edges)
            i, j = edge.i, edge.j
            X[k, i] += unitRow[edge.edge] / 2
            X[k, j] += unitRow[edge.edge] / 2
        end
    end

    ## run partial regression
    X = (transpose(X) * X)^-1 * transpose(X)

    ## fit each dimension (edge) of the original high-dimensional space
    for edge in model.edges.edge
        y = Vector{Float64}(model.accum[:, edge])
        coefs = X * y # regress on this edge
        model.nodes[!, edge] = coefs[1:end]
    end

    # find accumHat
    ## Refit the units: in high-d space, the refit units are as close as possible to their
    ## center of mass wrt the network
    ## This is by-definition the refit space
    for unitEdge in model.edges.edge
        for (k, unitRow) in enumerate(eachrow(model.accumHat))
            unitRow[unitEdge] =
                sum(
                    model.nodes[nodeEdge.i, unitEdge] * model.accum[k, nodeEdge.edge] +
                    model.nodes[nodeEdge.j, unitEdge] * model.accum[k, nodeEdge.edge]
                    for nodeEdge in eachrow(model.edges)
                ) / 2
        end
    end
end

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    # for existing dimensions given to us by the child class,
    # add normalize them, orthogonalize them from each other,
    # and add point positions to the model
    edgeNames = model.edges.edge
    numExistingDims = nrow(model.embedding)
    if numExistingDims > 0
        # normalize the first dimension (the rest are normalized below)
        v1 = Vector(model.embedding[1, edgeNames])
        s = sqrt(sum(v1 .^ 2))
        if s != 0
            model.embedding[1, edgeNames] .= v1 / s
        end
    end

    for i in 1:numExistingDims
        addPointsToModelFromDim(model, model.embedding[i, edgeNames])

        # orthogonalize, by rejection, each existing dimension from each previous dimension.
        # note, a dimension will be orthogonalize before it is used to add points to the model
        vi = Vector(model.embedding[i, edgeNames])
        denom = dot(vi, vi)
        for j in (i+1):numExistingDims
            vj = Vector(model.embedding[j, edgeNames])
            s = sqrt(sum(vj .^ 2))
            if s != 0
                model.embedding[j, edgeNames] .= vj / s
            end

            scalar = dot(vj, vi) / denom
            model.embedding[j, edgeNames] -= scalar * vi
            s = sqrt(sum(vj .^ 2))
            if s < 0.05
                model.embedding[j, edgeNames] .= 0
                @warn "During the rotation step, axis $j was deflated to zero due to close correlation with axis $i."
            elseif s != 0
                model.embedding[j, edgeNames] .= vj / s
            end
        end
    end

    # make a mean centered copy of accum
    edgeNames = model.edges.edge
    X = Matrix{Float64}(model.accum[!, edgeNames])
    for i in 1:size(X)[2]
        X[:, i] .-= mean(X[:, i])
    end

    # if needed, deflate those points' dimensions from X before we run SVD on it
    if numExistingDims > 0
        P = transpose(Matrix{Float64}(model.points))
        for i in 1:size(P)[2]
            col = P[:, i] .- mean(P[:, i])
            denom = dot(col, col)
            for j in 1:size(X)[2]
                scalar = dot(X[:, j], col) / denom
                X[:, j] -= scalar * col
            end
        end
    end

    # then, once we've deflated or not, we run SVD on the data, then add to the model
    svd = transpose(projection(fit(PCA, X', pratio=1.0)))
    df = similar(model.embedding, size(svd)[1])
    df[!, edgeNames] = svd
    df.label = ["SVD$(i)" for i in 1:nrow(df)]
    append!(model.embedding, df)
    for i in (numExistingDims+1):nrow(model.embedding)
        addPointsToModelFromDim(model, model.embedding[i, :])
    end

    # Now that we have the full embedding ready, run basic stats on it
    unitIDs = model.accum.unitID
    total_variance = sum(var.(eachrow(model.points[!, unitIDs])))
    for i in 1:nrow(model.embedding)
        points = Vector(model.points[i, unitIDs])
        pointsHat = Vector(model.pointsHat[i, unitIDs])
        model.embedding[i, :variance_explained] = var(points) / total_variance
        model.embedding[i, :pearson] = cor(points, pointsHat)
        pointsDiffs = [
            a - b
            for a in points
            for b in points
        ]

        pointsHatDiffs = [
            a - b
            for a in pointsHat
            for b in pointsHat
        ]

        model.embedding[i, :coregistration] = cor(pointsDiffs, pointsHatDiffs)
    end
end

# function tests(
#         ::Type{M}, model::AbstractLinearENAModel
#     ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

#     # TODO maybe cut this, think about tests as summary stats dict
# end