# TODO make a few different ENAModel overlap wrappers
# TODO make this more memory efficient by adding a unitMetadata df that holds the overlap between the two unit models

"""
TODO: document
"""
struct ENAModel{T} <: AbstractENAModel{T}
    data::DataFrame # the raw data
    codes::Array{Symbol,1}
    conversations::Array{Symbol,1}
    units::Array{Symbol,1}
    windowSize::Int
    rotation::T
    countModel::DataFrame # all the unit-level data we compute
    networkModel::DataFrame # all the connections-level data we compute
    codeModel::DataFrame # all the code-level data we compute
    centroidModel::DataFrame # countModel with re-approximated relationship columns
    relationshipMap::Any
    pearson::Real
end

function ENAModel(data::DataFrame, codes::Array{Symbol,1}, conversations::Array{Symbol,1}, units::Array{Symbol,1};
    windowSize::Int=4, rotateBy::T=SVDRotation(), sphereNormalize::Bool=true, dropEmpty::Bool=true,
    subsetFilter::Function=x->true) where {T<:AbstractENARotation}

    # Preparing model structures
    ## Relationships between codes
    relationshipMap = Dict(Symbol(string(code1, "_", code2)) => (i, j)
                         for (i, code1) in enumerate(codes)
                         for (j, code2) in enumerate(codes)
                         if i < j)

    ## Adding a new column to the raw data that labels each unit properly
    unitJoinedData = hcat(data, DataFrame(:ENA_UNIT => map(eachrow(data)) do dataRow
        return join(dataRow[units], ".")
    end))

    ## Unit model placeholders
    countModel = combine(first, groupby(unitJoinedData, :ENA_UNIT))
    countModel = hcat(countModel, DataFrame(Dict(r => Real[0 for i in 1:nrow(countModel)] # my relative strength for each relationship
                                               for r in keys(relationshipMap))))

    countModel = hcat(countModel, DataFrame(pos_x=Real[0 for i in 1:nrow(countModel)], # my position on the x and y axes
                                          pos_y=Real[0 for i in 1:nrow(countModel)])) 

    ## Network model placeholders
    networkModel = DataFrame(relationship=collect(keys(relationshipMap)),
                             density=Real[0 for r in relationshipMap], # how thick to make the line
                             weight_x=Real[0 for r in relationshipMap], # the weight I contribute to dim_x's
                             weight_y=Real[0 for r in relationshipMap]) # the weight I contribute to dim_y's
    
    ## Code model placeholders
    codeModel = DataFrame(code=codes,
                          density=Real[0 for c in codes], # how thick to make the dot
                          pos_x=Real[0 for c in codes], # where to plot this code on the fitted plot's x-axis
                          pos_y=Real[0 for c in codes]) # where to plot this code on the fitted plot's y-axis

    codeModel = hcat(codeModel, DataFrame(Dict(r => Real[0 for i in 1:nrow(codeModel)]
                     for r in keys(relationshipMap)))) # my position on the approximated x-axis if we use a relationship as the x-axis

    # Accumulation step
    ## Raw counts for all the units
    counts = Dict(unit => [[0 for j in codes] for i in codes]
                  for unit in countModel[!, :ENA_UNIT])

    prev_convo = unitJoinedData[1, conversations]
    howrecents = [Inf for c in codes]
    for line in eachrow(unitJoinedData)
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
        for r in keys(relationshipMap)
            i, j = relationshipMap[r]
            if howrecents[i] == 0 && howrecents[j] < windowSize
                counts[unit][i][j] += 1
            elseif howrecents[j] == 0 && howrecents[i] < windowSize
                counts[unit][i][j] += 1
            end
        end
    end

    ## Normalize and overwrite the unit model's placeholders
    if sphereNormalize
        for unitRow in eachrow(countModel)
            unit = unitRow[:ENA_UNIT]
            vector = [counts[unit][i][j] for (i,j) in values(relationshipMap)]
            s = sqrt(sum(vector .^ 2))
            if s != 0
                vector /= s
            end

            for (k, r) in enumerate(keys(relationshipMap))
                unitRow[r] = vector[k]
            end
        end
    end

    # Filtering
    ## User-defined unit subsetting - we kept all the data in the count model up to this point so the user can define filters as they please
    filter!(subsetFilter, countModel)

    ## Drop empty values
    if dropEmpty
        filter!(countModel) do unitRow
            return !all(iszero.(values(unitRow[networkModel[!, :relationship]])))
        end
    end

    ## Now that we have all the data counted, copy it into the centroid model, then narrow the count model down to just the cols that matter there
    centroidModel = countModel
    countModel = copy(countModel[!, [:ENA_UNIT, :pos_x, :pos_y, networkModel[!, :relationship]...]])

    # Compute the position of the codes in the approximated high-dimensional space
    ## Compute the density of each network row line
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        networkRow[:density] = sum(countModel[!, r])
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
    X = Matrix{Float64}(zeros(nrow(countModel), 1 + nrow(codeModel)))
    X[:, 1] .= 1 # Intercept term
    for (i, unitRow) in enumerate(eachrow(countModel))
        for r in keys(relationshipMap)
            a, b = relationshipMap[r]
            X[i, a+1] += unitRow[r] / 2 # +1 because of intercept
            X[i, b+1] += unitRow[r] / 2 # +1 because of intercept
        end
    end

    X = (transpose(X) * X)^-1 * transpose(X)

    ## Fit each dimension of the original high-dimensional space
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        y = Vector{Float64}(countModel[:, r])
        r_coefs = X * y
        codeModel[:, r] = r_coefs[2:end] # 2:end because of intercept
    end

    ## Refit the units: in high-d space, the refit units are as close as possible to their center of mass wrt the network
    ## This is by-definition the refit space
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        for (i, centroidRow) in enumerate(eachrow(centroidModel))
            centroidRow[r] =
                sum(
                    codeModel[relationshipMap[k][1], r] * countModel[i, k] +
                    codeModel[relationshipMap[k][2], r] * countModel[i, k]
                    for k in keys(relationshipMap)
                ) / 2
        end
    end

    # Rotation step
    ## Use the given lambda, probably one of the out-of-the-box ENARotations, but could be anything user defined
    rotate!(rotateBy, networkModel, centroidModel)

    # Layout step
    ## Project the pos_x and pos_y for the original units onto the plane, now that we have the rotation
    ## This is really only for testing the goodness of fit
    for unitRow in eachrow(countModel)
        unitRow[:pos_x] = sum(
            networkRow[:weight_x] * unitRow[networkRow[:relationship]]
            for networkRow in eachrow(networkModel)
        )

        unitRow[:pos_y] = sum(
            networkRow[:weight_y] * unitRow[networkRow[:relationship]]
            for networkRow in eachrow(networkModel)
        )
    end

    ## Same for the refit units
    ## These are what are really drawn
    for centroidRow in eachrow(centroidModel)
        centroidRow[:pos_x] = sum(
            networkRow[:weight_x] * centroidRow[networkRow[:relationship]]
            for networkRow in eachrow(networkModel)
        )

        centroidRow[:pos_y] = sum(
            networkRow[:weight_y] * centroidRow[networkRow[:relationship]]
            for networkRow in eachrow(networkModel)
        )
    end

    ## Same for the codes
    ## These aren't used to compute what's really drawn, they are labels floating around
    ## in the same high-d space as what is really drawn, that we interpret in terms of the
    ## center of mass. If we *were* to use these to compute what's really drawn,
    ## it should actually give us the same result as the projection we used for centroidRow's above,
    ## since that's the property we defined the refit space to have. (ignoring the intercept)
    for codeRow in eachrow(codeModel)
        codeRow[:pos_x] = sum(
            networkRow[:weight_x] * codeRow[networkRow[:relationship]]
            for networkRow in eachrow(networkModel)
        )

        codeRow[:pos_y] = sum(
            networkRow[:weight_y] * codeRow[networkRow[:relationship]]
            for networkRow in eachrow(networkModel)
        )
    end

    ## Translate everything so overall mean lies at the origin
    centroidModel[!, :pos_x] = centroidModel[!, :pos_x] .- mean(centroidModel[!, :pos_x])
    centroidModel[!, :pos_y] = centroidModel[!, :pos_y] .- mean(centroidModel[!, :pos_y])
    codeModel[!, :pos_x] = codeModel[!, :pos_x] .- mean(centroidModel[!, :pos_x])
    codeModel[!, :pos_y] = codeModel[!, :pos_y] .- mean(centroidModel[!, :pos_y])
    countModel[!, :pos_x] = countModel[!, :pos_x] .- mean(countModel[!, :pos_x])
    countModel[!, :pos_y] = countModel[!, :pos_y] .- mean(countModel[!, :pos_y])

    # Testing step
    ## Run tests for how well the refit space coregisters with the original space
    fitDiffs = Real[]
    dimDiffs = Real[]
    for (i, unitRowA) in enumerate(eachrow(countModel))
        for (j, unitRowB) in enumerate(eachrow(countModel))
            if i < j
                push!(fitDiffs, centroidModel[i, :pos_x] - centroidModel[j, :pos_x])
                push!(fitDiffs, centroidModel[i, :pos_y] - centroidModel[j, :pos_y])
                push!(dimDiffs, unitRowA[:pos_x] - unitRowB[:pos_x])
                push!(dimDiffs, unitRowA[:pos_y] - unitRowB[:pos_y])
            end
        end
    end

    pearson = cor(fitDiffs, dimDiffs)
    # p = pvalue(OneSampleTTest(fitDiffs, dimDiffs))
    # pvalue(EqualVarianceTTest(x, y))
    # pvalue(UnequalVarianceTTest(x, y))
    # pvalue(MannWhitneyUTest(x, y))
    # pvalue(SignedRankTest(x, y))

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
    return ENAModel(unitJoinedData, codes, conversations, units,
                    windowSize, rotateBy, countModel, networkModel, codeModel, centroidModel,
                    relationshipMap, pearson)
end