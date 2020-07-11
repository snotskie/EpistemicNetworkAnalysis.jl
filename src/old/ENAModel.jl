# TODO make a few different ENAModel overlap wrappers
# TODO make this more memory efficient by adding a unitMetadata df that holds the overlap between the two unit models

"""
TODO: document
"""
struct ENAModel
    data::DataFrame # the raw data
    codes::Array{Symbol,1}
    conversations::Array{Symbol,1}
    units::Array{Symbol,1}
    windowSize::Int
    rotation::ENARotation
    unitModel::DataFrame # all the unit-level data we compute
    networkModel::DataFrame # all the connections-level data we compute
    codeModel::DataFrame # all the code-level data we compute
    refitUnitModel::DataFrame # unitModel with re-approximated relationship columns
    relationshipMap::Any
    pvalue::Real
    pearson::Real
    variance_x::Real
    variance_y::Real
end

function ENAModel(data::DataFrame, codes::Array{Symbol,1}, conversations::Array{Symbol,1}, units::Array{Symbol,1};
    windowSize::Int=4, rotateBy::T=SVDRotation(), sphereNormalize::Bool=true, dropEmpty::Bool=true,
    subsetFilter::Function=x->true) where {T<:ENARotation}

    # Preparing model structures
    ## Relationships between codes
    relationships = Dict(Symbol(string(code1, "_", code2)) => (i, j)
                         for (i, code1) in enumerate(codes)
                         for (j, code2) in enumerate(codes)
                         if i < j)

    ## Adding a new column to the raw data that labels each unit properly
    unitVar = :ENA_UNIT
    unitJoinedData = hcat(data, DataFrame(unitVar => map(eachrow(data)) do dataRow
        return join(dataRow[units], ".")
    end))

    ## Unit model placeholders
    unitModel = combine(first, groupby(unitJoinedData, unitVar))
    unitModel = hcat(unitModel, DataFrame(Dict(r => Real[0 for i in 1:nrow(unitModel)] # my relative strength for each relationship
                                               for r in keys(relationships))))

    unitModel = hcat(unitModel, DataFrame(pos_x=Real[0 for i in 1:nrow(unitModel)], # my position on the x and y axes
                                          pos_y=Real[0 for i in 1:nrow(unitModel)])) 

    ## Network model placeholders
    networkModel = DataFrame(relationship=collect(keys(relationships)),
                             density=Real[0 for r in relationships], # how thick to make the line
                             weight_x=Real[0 for r in relationships], # the weight I contribute to dim_x's
                             weight_y=Real[0 for r in relationships]) # the weight I contribute to dim_y's
    
    ## Code model placeholders
    codeModel = DataFrame(code=codes,
                          density=Real[0 for c in codes], # how thick to make the dot
                          pos_x=Real[0 for c in codes], # where to plot this code on the fitted plot's x-axis
                          pos_y=Real[0 for c in codes]) # where to plot this code on the fitted plot's y-axis

    codeModel = hcat(codeModel, DataFrame(Dict(r => Real[0 for i in 1:nrow(codeModel)]
                     for r in keys(relationships)))) # my position on the approximated x-axis if we use a relationship as the x-axis

    # Accumulation step
    ## Raw counts for all the units
    counts = Dict(unit => [[0 for j in codes] for i in codes]
                  for unit in unitModel[!, unitVar])

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

        unit = line[unitVar]
        for r in keys(relationships)
            i, j = relationships[r]
            if howrecents[i] == 0 && howrecents[j] < windowSize
                counts[unit][i][j] += 1
            elseif howrecents[j] == 0 && howrecents[i] < windowSize
                counts[unit][i][j] += 1
            end
        end
    end

    ## Normalize and overwrite the unit model's placeholders
    if sphereNormalize
        for unitRow in eachrow(unitModel)
            unit = unitRow[unitVar]
            vector = [counts[unit][i][j] for (i,j) in values(relationships)]
            s = sqrt(sum(vector .^ 2))
            if s != 0
                vector /= s
            end

            for (k, r) in enumerate(keys(relationships))
                unitRow[r] = vector[k]
            end
        end
    end

    # Filtering
    ## User-defined unit subsetting
    filter!(subsetFilter, unitModel)

    ## Drop empty values
    if dropEmpty
        filter!(unitModel) do unitRow
            return !all(iszero.(values(unitRow[networkModel[!, :relationship]])))
        end
    end

    ## Copy the refit model here since now we have all the data in the original unit model
    refitUnitModel = copy(unitModel)

    # Compute the position of the codes in the approximated high-dimensional space
    ## Compute the density of each network row line
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        networkRow[:density] = sum(unitModel[!, r])
    end

    ## Normalize the network densities
    s = maximum(networkModel[!, :density])
    networkModel[!, :density] /= s

    ## Compute the density of each code row dot, by "splitting" the density of each line between its two codes
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        i, j = relationships[r]
        codeModel[i, :density] += networkRow[:density] # TODO check, is this the right way to do this?
        codeModel[j, :density] += networkRow[:density]
    end

    ## Normalize the code densities
    s = maximum(codeModel[!, :density])
    codeModel[!, :density] /= s

    ## Regression model for placing the code dots into the approximated high-dimensional space
    X = Matrix{Float64}(zeros(nrow(unitModel), 1 + nrow(codeModel)))
    X[:, 1] .= 1 # Intercept term
    for (i, unitRow) in enumerate(eachrow(unitModel))
        for r in keys(relationships)
            a, b = relationships[r]
            X[i, a+1] += unitRow[r] / 2 # +1 because of intercept
            X[i, b+1] += unitRow[r] / 2 # +1 because of intercept
        end
    end

    X = (transpose(X) * X)^-1 * transpose(X)

    ## Fit each dimension of the original high-dimensional space
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        y = Vector{Float64}(unitModel[:, r])
        r_coefs = X * y
        codeModel[:, r] = r_coefs[2:end] # 2:end because of intercept
    end

    ## Refit the units: in high-d space, the refit units are as close as possible to their center of mass wrt the network
    ## This is by-definition the refit space
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        for (i, refitRow) in enumerate(eachrow(refitUnitModel))
            refitRow[r] =
                sum(
                    codeModel[relationships[k][1], r] * unitModel[i, k] +
                    codeModel[relationships[k][2], r] * unitModel[i, k]
                    for k in keys(relationships)
                ) / 2
        end
    end

    # Rotation step
    ## Use the given lambda, probably one of the out-of-the-box ENARotations, but could be anything user defined
    rotate!(rotateBy, networkModel, refitUnitModel)

    # Layout step
    ## Project the pos_x and pos_y for the original units onto the plane, now that we have the rotation
    ## This is really only for testing the goodness of fit
    for unitRow in eachrow(unitModel)
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
    for refitRow in eachrow(refitUnitModel)
        refitRow[:pos_x] = sum(
            networkRow[:weight_x] * refitRow[networkRow[:relationship]]
            for networkRow in eachrow(networkModel)
        )

        refitRow[:pos_y] = sum(
            networkRow[:weight_y] * refitRow[networkRow[:relationship]]
            for networkRow in eachrow(networkModel)
        )
    end

    ## Same for the codes
    ## These aren't used to compute what's really drawn, they are labels floating around
    ## in the same high-d space as what is really drawn, that we interpret in terms of the
    ## center of mass. If we *were* to use these to compute what's really drawn,
    ## it should actually give us the same result as the projection we used for refitRow's above,
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
    refitUnitModel[!, :pos_x] = refitUnitModel[!, :pos_x] .- mean(refitUnitModel[!, :pos_x])
    refitUnitModel[!, :pos_y] = refitUnitModel[!, :pos_y] .- mean(refitUnitModel[!, :pos_y])
    codeModel[!, :pos_x] = codeModel[!, :pos_x] .- mean(refitUnitModel[!, :pos_x])
    codeModel[!, :pos_y] = codeModel[!, :pos_y] .- mean(refitUnitModel[!, :pos_y])
    unitModel[!, :pos_x] = unitModel[!, :pos_x] .- mean(unitModel[!, :pos_x])
    unitModel[!, :pos_y] = unitModel[!, :pos_y] .- mean(unitModel[!, :pos_y])

    # Testing step
    ## Run tests for how well the refit space coregisters with the original space
    fitDiffs = Real[]
    dimDiffs = Real[]
    for (i, unitRowA) in enumerate(eachrow(unitModel))
        for (j, unitRowB) in enumerate(eachrow(unitModel))
            if i < j
                push!(fitDiffs, refitUnitModel[i, :pos_x] - refitUnitModel[j, :pos_x])
                push!(fitDiffs, refitUnitModel[i, :pos_y] - refitUnitModel[j, :pos_y])
                push!(dimDiffs, unitRowA[:pos_x] - unitRowB[:pos_x])
                push!(dimDiffs, unitRowA[:pos_y] - unitRowB[:pos_y])
            end
        end
    end

    p = pvalue(OneSampleTTest(fitDiffs, dimDiffs))
    pearson = cor(fitDiffs, dimDiffs)
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

    ## Compute the variance explained by the axes
    total_variance = sum(var.(eachcol(refitUnitModel[!, networkModel[!, :relationship]])))
    variance_x = var(refitUnitModel[!, :pos_x]) / total_variance
    variance_y = var(refitUnitModel[!, :pos_y]) / total_variance

    # Done!
    return ENAModel(unitJoinedData, codes, conversations, units,
                    windowSize, rotateBy, unitModel, networkModel, codeModel, refitUnitModel,
                    relationships, p, pearson, variance_x, variance_y)
end