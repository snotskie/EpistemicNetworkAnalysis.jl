struct ENAModel{T} <: AbstractENAModel{T}
    # Inherits:
    codes::Array{Symbol,1}
    conversations::Array{Symbol,1}
    units::Array{Symbol,1}
    rotation::T
    rotateOn::Symbol
    accumModel::DataFrame # all the unit-level data we compute
    centroidModel::DataFrame # accumModel with re-approximated relationship columns
    metadata::DataFrame
    codeModel::DataFrame # all the code-level data we compute
    networkModel::DataFrame # all the connections-level data we compute
    relationshipMap::Dict{Symbol,Tuple{Int,Int}}

    # Adds:
    windowSize::Int
end

function ENAModel(
        # required
        data::DataFrame, codes::Array{Symbol,1}, conversations::Array{Symbol,1}, units::Array{Symbol,1};
        # optional
        windowSize::Int=4, rotateBy::AbstractENARotation=SVDRotation(), rotateOn::Symbol=:accumModel,
        sphereNormalize::Bool=true, dropEmpty::Bool=false, recenterEmpty::Bool=false,
        deflateEmpty::Bool=false, meanCenter::Bool=true, subspaces::Int=0,
        subsetFilter::Function=x->true, relationshipFilter::Function=(i,j,ci,cj)->(i<j)
    )

    # Checking that the options are sane
    if windowSize < 1
        error("The windowSize must be positive")
    elseif !sphereNormalize && deflateEmpty
        error("When sphereNormalize=false, deflateEmpty=true has no interpretive validity")
    elseif dropEmpty && recenterEmpty
        error("When dropEmpty=true, recenterEmpty=true has no interpretive validity")
    elseif length(codes) < 3
        error("There must be at least 3 codes in the model")
    elseif length(conversations) < 1
        error("Need at least one column to define conversations")
    elseif length(units) < 1
        error("Need at least one column to define units")
    elseif !(rotateOn in [:centroidModel, :accumModel])
        error("rotateOn must be either :centroidModel or :accumModel")
    end

    # Preparing model structures
    ## Relationships between codes
    relationshipMap = Dict{Symbol,Tuple{Int,Int}}()
    for (i, code1) in enumerate(codes)
        for (j, code2) in enumerate(codes)
            if relationshipFilter(i, j, code1, code2)
                key = Symbol(string(code1, "_", code2))
                relationshipMap[key] = (i, j)
            end
        end
    end

    relationshipIds = collect(keys(relationshipMap))

    ## Adding a new column to the raw data that labels each unit properly
    ena_ids = map(eachrow(data)) do dataRow
        return join(dataRow[units], ".")
    end

    data = hcat(data, DataFrame(:ENA_UNIT => ena_ids))

    ## Unit model placeholders
    accumModel = combine(first, groupby(data, :ENA_UNIT))
    emptyValues = DataFrame(Dict(
        r => Real[0 for i in 1:nrow(accumModel)]
        for r in relationshipIds
    ))

    emptyPositions = DataFrame(
        pos_x=Real[0 for i in 1:nrow(accumModel)],
        pos_y=Real[0 for i in 1:nrow(accumModel)]
    )

    accumModel = hcat(accumModel, emptyValues)
    accumModel = hcat(accumModel, emptyPositions) 
    
    ## Replace codes "first" metadata with sums instead (ie, code and count)
    sums = combine(
        groupby(data[!, [codes..., :ENA_UNIT]], :ENA_UNIT),
        (code => sum for code in codes)...
    )
    
    for code in codes
        accumModel[!, code] = sums[!, Symbol(string(code, "_sum"))]
    end

    ## Network model placeholders
    networkModel = DataFrame(
        relationship=relationshipIds,
        response=[i for (i, j) in values(relationshipMap)],
        referent=[j for (i, j) in values(relationshipMap)],
        density=Real[0 for r in relationshipIds], # how thick to make the line
        weight_x=Real[0 for r in relationshipIds], # the weight I contribute to dim_x's
        weight_y=Real[0 for r in relationshipIds] # the weight I contribute to dim_y's
    )
    
    ## Code model placeholders
    codeModel = DataFrame(
        code=codes,
        density=Real[0 for c in codes], # how thick to make the dot
        pos_x=Real[0 for c in codes], # where to plot this code on the fitted plot's x-axis
        pos_y=Real[0 for c in codes] # where to plot this code on the fitted plot's y-axis
    )

    emptyValues2 = DataFrame(Dict(
        r => Real[0 for c in codes] # my position on the approximated x-axis if we use a relationship as the x-axis
        for r in relationshipIds
    ))

    codeModel = hcat(codeModel, emptyValues2)

    # Accumulation step
    ## Raw counts for all the units
    counts = Dict(
        unit => [[0 for j in codes] for i in codes]
        for unit in accumModel[!, :ENA_UNIT]
    )

    prev_convo = data[1, conversations]
    howrecents = [Inf for c in codes]
    for line in eachrow(data)
        if prev_convo != line[conversations]
            prev_convo = line[conversations]
            howrecents .= Inf
        end

        for (i, code) in enumerate(codes)
            if line[code] > 0
                howrecents[i] = 0
            else
                howrecents[i] += 1
            end
        end

        unit = line[:ENA_UNIT]
        for r in relationshipIds
            i, j = relationshipMap[r]
            if howrecents[i] == 0 && howrecents[j] < windowSize
                counts[unit][i][j] += 1
            elseif howrecents[j] == 0 && howrecents[i] < windowSize
                counts[unit][i][j] += 1
            end
        end
    end

    ## Overwrite the accum model's placeholders
    for unitRow in eachrow(accumModel)
        unit = unitRow[:ENA_UNIT]
        vector = [counts[unit][i][j] for (i,j) in values(relationshipMap)]
        for (k, r) in enumerate(relationshipIds)
            unitRow[r] = vector[k]
        end
    end

    # Filtering
    ## User-defined unit subsetting - we kept all the data in the count model up to this point so the user can define filters as they please
    filter!(subsetFilter, accumModel)

    ## Maybe drop rows with empty values or resposition them
    if dropEmpty
        ### Find values to keep
        nonZeroRows = map(eachrow(accumModel)) do unitRow
            return !all(iszero.(values(unitRow[relationshipIds])))
        end

        ### reassign
        accumModel = accumModel[nonZeroRows, :]
    elseif recenterEmpty
        ### Find values that actually have value
        nonZeroRows = map(eachrow(accumModel)) do unitRow
            return !all(iszero.(values(unitRow[relationshipIds])))
        end

        ### Count them and find the mean point
        N = sum(nonZeroRows)
        meanPoint = map(relationshipIds) do r
            return sum(accumModel[!, r]) / N
        end

        ### For all rows, if they are empty, replace with the mean point
        for unitRow in eachrow(accumModel)
            if all(iszero.(values(unitRow[relationshipIds])))
                for (k, r) in enumerate(relationshipIds)
                    unitRow[r] = meanPoint[k]
                end
            end
        end
    end

    ## Normalize each unit, if requested
    if sphereNormalize
        for unitRow in eachrow(accumModel)
            vector = Vector{Float64}(unitRow[relationshipIds])
            s = sqrt(sum(vector .^ 2))
            if s != 0
                unitRow[relationshipIds] = vector / s
            end
        end
    end

    ## If we didn't sphere normalize, then still scale it down make plots easier to read
    if !sphereNormalize
        s = maximum(maximum(accumModel[!, r]) for r in relationshipIds)
        for r in relationshipIds
            accumModel[!, r] /= s
        end
    end

    # Now that we have all the data counted, divvy and copy it between accumModel, centroidModel, and metadata
    modelColumns = [:pos_x, :pos_y, relationshipIds...]
    nonModelColumns = setdiff(Symbol.(names(accumModel)), modelColumns)
    metadata = accumModel[!, nonModelColumns]
    centroidModel = copy(accumModel[!, [:ENA_UNIT, modelColumns...]])
    accumModel = accumModel[!, [:ENA_UNIT, modelColumns...]]

    # Compute the position of the codes in the approximated high-dimensional space
    ## Compute the density of each network row line
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        networkRow[:density] = sum(accumModel[!, r])
    end

    ## Normalize the network densities
    s = maximum(networkModel[!, :density])
    networkModel[!, :density] /= s

    ## Compute the density of each code row dot, by "splitting" the density of each line between its two codes
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        i, j = relationshipMap[r]
        codeModel[i, :density] += networkRow[:density]
        codeModel[j, :density] += networkRow[:density]
    end

    ## Normalize the code densities
    s = maximum(codeModel[!, :density])
    codeModel[!, :density] /= s

    ## Regression model for placing the code dots into the approximated high-dimensional space
    X = Matrix{Float64}(rand(nrow(accumModel), nrow(codeModel)) / 1000000000)
    for (i, unitRow) in enumerate(eachrow(accumModel))
        for r in relationshipIds
            a, b = relationshipMap[r]
            X[i, a] += unitRow[r] / 2
            X[i, b] += unitRow[r] / 2
        end
    end

    X = (transpose(X) * X)^-1 * transpose(X)

    ## Fit each dimension of the original high-dimensional space
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        y = Vector{Float64}(accumModel[:, r])
        r_coefs = X * y
        codeModel[:, r] = r_coefs[1:end]
    end

    ## Refit the units: in high-d space, the refit units are as close as possible to their center of mass wrt the network
    ## This is by-definition the refit space
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        for (i, unitRow) in enumerate(eachrow(centroidModel))
            unitRow[r] =
                sum(
                    codeModel[relationshipMap[k][1], r] * accumModel[i, k] +
                    codeModel[relationshipMap[k][2], r] * accumModel[i, k]
                    for k in relationshipIds
                ) / 2
        end
    end

    # Rotation step
    ## Use the given rotation method, probably one of the out-of-the-box ENARotations, but could be anything user defined
    ## But first, maybe deflate the rotation model
    subspaceModel = accumModel
    if rotateOn == :centroidModel
        subspaceModel = centroidModel
    end

    if deflateEmpty
        deflatedModel = copy(subspaceModel)
        zAxis = map(eachrow(networkModel)) do networkRow
            r = networkRow[:relationship]
            return sum(subspaceModel[!, r])
        end

        s = sqrt(sum(zAxis .^ 2))
        if s != 0
            zAxis /= s
        end

        deflatedModel[!, :pos_z] = Matrix{Float64}(subspaceModel[!, relationshipIds]) * Vector{Float64}(zAxis)
        for r in relationshipIds
            rAxis = map(relationshipIds) do rp
                if r == rp
                    return 1
                else
                    return 0
                end
            end

            scalar = dot(rAxis, zAxis) / dot(zAxis, zAxis)
            deflatedModel[!, r] = subspaceModel[!, r] .- (scalar * deflatedModel[!, :pos_z])
        end

        subspaceModel = deflatedModel
    end

    if subspaces >= 2
        deflatedModel = copy(subspaceModel)
        pcaModel = projection(help_deflating_svd(networkModel, subspaceModel))
        for r in relationshipIds
            rAxis = map(relationshipIds) do rp
                if r == rp
                    return 1
                else
                    return 0
                end
            end

            # rAxisNew = rAxis * 0
            deflatedModel[!, r] = deflatedModel[!, r] * 0
            for i in 1:min(subspaces, size(pcaModel, 2))
                zAxis = pcaModel[:, i]
                deflatedModel[!, :pos_z] = Matrix{Float64}(subspaceModel[!, relationshipIds]) * Vector{Float64}(zAxis)

                scalar = dot(rAxis, zAxis) / dot(zAxis, zAxis)
                deflatedModel[!, r] = deflatedModel[!, r] + scalar * deflatedModel[!, :pos_z]
            end
        end

        subspaceModel = deflatedModel
    end

    rotate!(rotateBy, networkModel, codeModel, metadata, subspaceModel)

    # Layout step
    ## Project the pos_x and pos_y for the original units onto the plane, now that we have the rotation
    ## This is really only for testing the goodness of fit
    accumModel[!, :pos_x] = Matrix{Float64}(accumModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_x])
    accumModel[!, :pos_y] = Matrix{Float64}(accumModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_y])

    ## Same for the refit units
    ## These are what are really drawn (unlike in rENA, where accumModel is what is drawn)
    centroidModel[!, :pos_x] = Matrix{Float64}(centroidModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_x])
    centroidModel[!, :pos_y] = Matrix{Float64}(centroidModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_y])

    ## Same for the codes
    ## These aren't used to compute what's really drawn, they are labels floating around
    ## in the same high-d space as what is really drawn, that we interpret in terms of the
    ## center of mass. If we *were* to use these to compute what's really drawn,
    ## it should actually give us the same result as the projection we used for centroidRow's above,
    ## since that's the property we defined the refit space to have. (ignoring the intercept)
    codeModel[!, :pos_x] = Matrix{Float64}(codeModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_x])
    codeModel[!, :pos_y] = Matrix{Float64}(codeModel[!, relationshipIds]) * Vector{Float64}(networkModel[!, :weight_y])

    ## Maybe translate everything so overall mean lies at the origin
    if meanCenter
        mu_x = mean(centroidModel[!, :pos_x])
        mu_y = mean(centroidModel[!, :pos_y])
        centroidModel[!, :pos_x] = centroidModel[!, :pos_x] .- mu_x
        centroidModel[!, :pos_y] = centroidModel[!, :pos_y] .- mu_y
        mu_x = mean(accumModel[!, :pos_x])
        mu_y = mean(accumModel[!, :pos_y])
        accumModel[!, :pos_x] = accumModel[!, :pos_x] .- mu_x
        accumModel[!, :pos_y] = accumModel[!, :pos_y] .- mu_y
    end

    # Testing step
    ## Test that the angle between the dimensions is 90 degrees
    theta = dot(networkModel[!, :weight_x], networkModel[!, :weight_y])
    theta /= sqrt(dot(networkModel[!, :weight_x], networkModel[!, :weight_x]))
    theta /= sqrt(dot(networkModel[!, :weight_y], networkModel[!, :weight_y]))
    angle = acos(theta) * 180 / pi
    if abs(angle-90) > 0.0001 # allow for a little approximation error
        @warn """The angle between the axes of this model is $(angle) degrees, when it should be 90.
This can lead to strange visual effects when plotting on orthogonal axes.
This can undermine interpreting betweenness between units.
This can undermind interpreting variance explained by the axes.
And this can cause problems with ENA's optimization algorithm fitting the codes and the lines."""
    end

    # Done!
    return ENAModel(
        codes, conversations, units, rotateBy, rotateOn,
        accumModel, centroidModel, metadata, codeModel, networkModel,
        relationshipMap,
        windowSize
    )
end