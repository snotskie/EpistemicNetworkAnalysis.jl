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
        lineNormalize::Bool=false,
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
        lineNormalize=lineNormalize,
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

    # sanity checks: errors
    if nrow(data) < 3
        error("The data::DataFrame parameter should have at least 3 rows")
    elseif length(codes) < 3
        error("The codes::Array parameter should have at least 3 items")
    elseif length(units) == 0
        error("The units::Array parameter should have at least 1 item") 
    end

    # note, empty conversations is allowed. amounts to a whole conversation model

    config = NamedTuple(config)
    # sanity checks: warnings
    if config.dropEmpty && config.recenterEmpty
        @warn "dropEmpty and recenterEmpty were both set to true. In this case, recenterEmpty is moot, as empty units are dropped, so they cannot be recentered."
    end
    if config.lineNormalize && !config.dropEmpty && !config.recenterEmpty
        @warn "When setting lineNormalize to true, you should also set either dropEmpty or recenterEmpty to true."
    end

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
        :eigen_value=>Union{Number,Missing}[],
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
    if model.config.sphereNormalize || model.config.lineNormalize
        for i in 1:nrow(model.accum)
            vector = Vector{Float64}(model.accum[i, edgeIDs])
            s = sqrt(sum(vector .^ 2))
            if s != 0
                model.accum[i, edgeIDs] = vector / s
            end
        end

        if model.config.lineNormalize
            # fulcrum = ones(length(edgeIDs)) ./ sqrt(length(edgeIDs)) # unit length in the direction of 1
            fulcrum = mean.(eachcol(model.accum[!, edgeIDs]))
            fulcrum /= sqrt(sum(fulcrum .^ 2))
            numer = dot(fulcrum, fulcrum)
            for i in 1:nrow(model.accum)
                vector = Vector{Float64}(model.accum[i, edgeIDs])
                denom = dot(fulcrum, vector)
                if denom != 0
                    vec_ext = vector * numer / denom # project fulcrum onto vector, "extending" vector in line to the fulcrum
                    new_dist = acos(denom) # new distance = arc length = angle in radians
                    old_dist = sqrt(sum((fulcrum - vec_ext) .^ 2))
                    model.accum[i, edgeIDs] = fulcrum + new_dist/old_dist*(vec_ext - fulcrum) # fulcrum, moved in the direction of vec_ext, by a distance equal to original arc length
                end
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
    @debug "numExistingDims = $(numExistingDims)"
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
        @debug "vi = $(vi)"
        denom = dot(vi, vi)
        if denom != 0
            @debug "denom = $(denom)"
            for j in (i+1):numExistingDims
                vj = Vector(model.embedding[j, edgeIDs])
                @debug "vj = $(vj)"
                s = sqrt(sum(vj .^ 2))
                if s != 0
                    model.embedding[j, edgeIDs] .= vj / s
                    vj = Vector(model.embedding[j, edgeIDs])
                end
                @debug "vj after scaling = $(vj)"

                scalar = dot(vj, vi) / denom
                @debug "scalar = $(scalar)"
                model.embedding[j, edgeIDs] .= vj - scalar * vi
                vj = Vector(model.embedding[j, edgeIDs])
                s = sqrt(sum(vj .^ 2))
                if s < 0.05
                    model.embedding[j, edgeIDs] .= 0
                    vj = Vector(model.embedding[j, edgeIDs])
                    @warn "During the rotation step, axis $j was deflated to zero due to close correlation with axis $i."
                elseif s != 0
                    model.embedding[j, edgeIDs] .= vj / s
                    vj = Vector(model.embedding[j, edgeIDs])
                end
                @debug "vj after rejection and rescaling = $(vj)"
                @debug "dot(vi, vj) = $(dot(vi, vj))"
            end
        end
    end

    # make a copy of accum, and if needed, deflate it before we run SVD on it
    X = Matrix{Float64}(model.accum[!, edgeIDs])
    for i in 1:numExistingDims
        # for j in axes(X, 2)
        #     X[:, j] .-= mean(X[:, j])
        # end

        vi = X * Vector{Float64}(model.embedding[i, edgeIDs])
        vi .-= mean(vi)
        denom = dot(vi, vi)
        if denom != 0
            for j in axes(X, 2)
                vj = X[:, j]
                vj .-= mean(vj)
                scalar = dot(vj, vi) / denom
                X[:, j] -= scalar * vi
            end
        end
    end

    # mean center X before SVD
    for i in axes(X, 2)
        X[:, i] .-= mean(X[:, i])
    end

    # then, once we've deflated or not, we run SVD on the data, then add to the model
    pca = fit(PCA, X', pratio=1.0, method=:svd) # BUGFIX force use svd, https://github.com/snotskie/EpistemicNetworkAnalysis.jl/issues/56#issuecomment-1910540698
    @debug "eigvals(pca) = $(eigvals(pca))"
    svd = transpose(projection(pca))
    # BUGFIX https://github.com/snotskie/EpistemicNetworkAnalysis.jl/issues/56#issuecomment-1910540698
    # Prevent SVD from adding more dimensions than are possible
    numSVDDims = min(size(svd, 1), length(edgeIDs) - numExistingDims)
    df = similar(model.embedding, numSVDDims)
    df[!, edgeIDs] = svd[1:numSVDDims, :]
    df.label = ["SVD$(i)" for i in 1:nrow(df)]
    df.eigen_value = eigvals(pca)[1:numSVDDims]
    append!(model.embedding, df)
    for i in (numExistingDims+1):nrow(model.embedding)
        addPointsToModelFromDim(model, i)
    end

    # ensure all values are defined in all dimensions
    for i in 1:nrow(model.embedding)
        for col in names(model.embedding)
            try
                model.embedding[i, col] = model.embedding[i, col]
            catch e
                if e isa UndefRefError
                    if eltype(model.embedding[:, col]) >: Missing
                        model.embedding[i, col] = missing
                    elseif eltype(model.embedding[:, col]) >: Nothing
                        model.embedding[i, col] = nothing
                    elseif eltype(model.embedding[:, col]) <: AbstractString
                        model.embedding[i, col] = "UNDEF"
                    elseif eltype(model.embedding[:, col]) >: Symbol
                        model.embedding[i, col] = :UNDEF
                    elseif eltype(model.embedding[:, col]) <: Number
                        model.embedding[i, col] = 0.0
                    elseif eltype(model.embedding[:, col]) >: Bool
                        model.embedding[i, col] = false
                    else
                        @error "An embedding value was undefined in a column for which a sensible 'missing' value does not exist. The type is $(eltype(model.embedding[:, col]))"
                    end
                else
                    rethrow(e)
                end
            end
        end
    end
end