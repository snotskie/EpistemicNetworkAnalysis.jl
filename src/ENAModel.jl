# TODO make a few different ENAModel overlap wrappers

"""
TODO: document
"""
struct ENAModel{T} <: AbstractENAModel{T}
    # Inherits:
    codes::Array{Symbol,1}
    conversations::Array{Symbol,1}
    units::Array{Symbol,1}
    rotation::T
    accumModel::DataFrame # all the unit-level data we compute
    centroidModel::DataFrame # accumModel with re-approximated relationship columns
    metadata::DataFrame
    codeModel::DataFrame # all the code-level data we compute
    networkModel::DataFrame # all the connections-level data we compute
    relationshipMap::Any

    # Adds:
    windowSize::Int
end

function ENAModel(data::DataFrame, codes::Array{Symbol,1}, conversations::Array{Symbol,1}, units::Array{Symbol,1};
    windowSize::Int=4, rotateBy::T=SVDRotation(), sphereNormalize::Bool=true, dropEmpty::Bool=true,
    subsetFilter::Function=x->true, accumulationMask::Function=x->1, blockingContext::Function=x->true) where {T<:AbstractENARotation}

    # Preparing model structures
    ## Relationships between codes
    relationshipMap = Dict(Symbol(string(code1, "_", code2)) => (i, j)
                         for (i, code1) in enumerate(codes)
                         for (j, code2) in enumerate(codes)
                         if i < j)

    ## Adding a new column to the raw data that labels each unit properly
    data = hcat(data, DataFrame(:ENA_UNIT => map(eachrow(data)) do dataRow
        return join(dataRow[units], ".")
    end))

    ## Unit model placeholders
    accumModel = combine(first, groupby(data, :ENA_UNIT))
    accumModel = hcat(accumModel, DataFrame(Dict(r => Real[0 for i in 1:nrow(accumModel)] # my relative strength for each relationship
                                               for r in keys(relationshipMap))))

    accumModel = hcat(accumModel, DataFrame(pos_x=Real[0 for i in 1:nrow(accumModel)], # my position on the x and y axes
                                          pos_y=Real[0 for i in 1:nrow(accumModel)])) 

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
                  for unit in accumModel[!, :ENA_UNIT])

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

        if blockingContext(howrecents)
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
    end

    ## Normalize and overwrite the unit model's placeholders
    for unitRow in eachrow(accumModel)
        unit = unitRow[:ENA_UNIT]
        vector = [counts[unit][i][j] for (i,j) in values(relationshipMap)]
        if sphereNormalize
            s = sqrt(sum(vector .^ 2))
            if s != 0
                vector /= s
            end
        end

        for (k, r) in enumerate(keys(relationshipMap))
            unitRow[r] = vector[k] * accumulationMask(r)
        end
    end

    # Filtering
    ## User-defined unit subsetting - we kept all the data in the count model up to this point so the user can define filters as they please
    filter!(subsetFilter, accumModel)

    ## Drop empty values
    if dropEmpty
        filter!(accumModel) do unitRow
            return !all(iszero.(values(unitRow[networkModel[!, :relationship]])))
        end
    end

    ## Now that we have all the data counted, divvy and copy it between accumModel, centroidModel, and metadata
    metadata = accumModel[!, setdiff(Symbol.(names(accumModel)), [:pos_x, :pos_y, networkModel[!, :relationship]...])]
    centroidModel = copy(accumModel[!, [:ENA_UNIT, :pos_x, :pos_y, networkModel[!, :relationship]...]])
    accumModel = accumModel[!, [:ENA_UNIT, :pos_x, :pos_y, networkModel[!, :relationship]...]]

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
    # X = Matrix{Float64}(zeros(nrow(accumModel), 1 + nrow(codeModel))) # +1 because of intercept
    X = Matrix{Float64}(zeros(nrow(accumModel), nrow(codeModel)))
    #X[:, 1] .= 1 # Intercept term
    for (i, unitRow) in enumerate(eachrow(accumModel))
        for r in keys(relationshipMap)
            a, b = relationshipMap[r]
            # X[i, a+1] += unitRow[r] / 2 # +1 because of intercept
            # X[i, b+1] += unitRow[r] / 2 # +1 because of intercept
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
        # codeModel[:, r] = r_coefs[2:end] # 2:end because of intercept
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
                    for k in keys(relationshipMap)
                ) / 2
        end
    end

    # Rotation step
    ## Use the given rotation method, probably one of the out-of-the-box ENARotations, but could be anything user defined
    rotate!(rotateBy, networkModel, centroidModel, metadata)

    # Layout step
    ## Project the pos_x and pos_y for the original units onto the plane, now that we have the rotation
    ## This is really only for testing the goodness of fit
    accumModel[!, :pos_x] = Matrix{Float64}(accumModel[!, networkModel[!, :relationship]]) * Vector{Float64}(networkModel[!, :weight_x])
    accumModel[!, :pos_y] = Matrix{Float64}(accumModel[!, networkModel[!, :relationship]]) * Vector{Float64}(networkModel[!, :weight_y])

    ## Same for the refit units
    ## These are what are really drawn
    centroidModel[!, :pos_x] = Matrix{Float64}(centroidModel[!, networkModel[!, :relationship]]) * Vector{Float64}(networkModel[!, :weight_x])
    centroidModel[!, :pos_y] = Matrix{Float64}(centroidModel[!, networkModel[!, :relationship]]) * Vector{Float64}(networkModel[!, :weight_y])

    ## Same for the codes
    ## These aren't used to compute what's really drawn, they are labels floating around
    ## in the same high-d space as what is really drawn, that we interpret in terms of the
    ## center of mass. If we *were* to use these to compute what's really drawn,
    ## it should actually give us the same result as the projection we used for centroidRow's above,
    ## since that's the property we defined the refit space to have. (ignoring the intercept)
    codeModel[!, :pos_x] = Matrix{Float64}(codeModel[!, networkModel[!, :relationship]]) * Vector{Float64}(networkModel[!, :weight_x])
    codeModel[!, :pos_y] = Matrix{Float64}(codeModel[!, networkModel[!, :relationship]]) * Vector{Float64}(networkModel[!, :weight_y])

    ## Translate everything so overall mean lies at the origin
    centroidModel[!, :pos_x] = centroidModel[!, :pos_x] .- mean(centroidModel[!, :pos_x])
    centroidModel[!, :pos_y] = centroidModel[!, :pos_y] .- mean(centroidModel[!, :pos_y])
    codeModel[!, :pos_x] = codeModel[!, :pos_x] .- mean(centroidModel[!, :pos_x])
    codeModel[!, :pos_y] = codeModel[!, :pos_y] .- mean(centroidModel[!, :pos_y])
    accumModel[!, :pos_x] = accumModel[!, :pos_x] .- mean(accumModel[!, :pos_x])
    accumModel[!, :pos_y] = accumModel[!, :pos_y] .- mean(accumModel[!, :pos_y])

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
        codes, conversations, units, rotateBy,
        accumModel, centroidModel, metadata, codeModel, networkModel,
        relationshipMap,
        windowSize
    )
end