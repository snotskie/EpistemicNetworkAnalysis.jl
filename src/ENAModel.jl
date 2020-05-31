struct ENAModel
    data::DataFrame # the raw data
    codes::Array{Symbol,1}
    conversations::Array{Symbol,1}
    units::Array{Symbol,1}
    windowSize::Int
    rotation::ENARotation
    # groupVar::Union{Nothing,Symbol}
    # treatmentGroup::Any # if groupVar == treatmentGroup, then 1
    # controlGroup::Any # if groupVar == treatmentGroup, then 0; if neither, leave it out
    # confounds::Union{Nothing,Array{Symbol,1}}
    # metadata::Array{Symbol,1}
    # rotateBy::Function
    # relationships::Any
    unitModel::DataFrame # all the unit-level data we compute
    networkModel::DataFrame # all the connections-level data we compute
    codeModel::DataFrame # all the code-level data we compute
    relationshipMap::Any
    pvalue::Real
end

function ENAModel(data::DataFrame, codes::Array{Symbol,1}, conversations::Array{Symbol,1}, units::Array{Symbol,1};
    windowSize::Int=4, rotateBy::T=SVDRotation(), sphereNormalize::Bool=true) where {T<:ENARotation}
# function ENAModel(data::DataFrame, codes::Array{Symbol,1}, conversations::Array{Symbol,1}, units::Array{Symbol,1};
#     metadata::Array{Symbol,1}=Symbol[], windowSize::Int=4, confounds::Union{Nothing,Array{Symbol,1}}=nothing,
#     groupVar::Union{Nothing,Symbol}=nothing, treatmentGroup::Any=nothing, controlGroup::Any=nothing,
#     rotateBy::Function=svd_rotation!)

    # Preparing model structures
    ## Relationships between codes
    relationships = Dict(Symbol(string(code1, "_", code2)) => (i, j)
                         for (i, code1) in enumerate(codes)
                         for (j, code2) in enumerate(codes)
                         if i < j)

    ## Unit model
    unitVar = :ENA_UNIT
    unitJoinedData = hcat(data, DataFrame(unitVar => map(eachrow(data)) do dataRow
        return join(dataRow[units], ".")
    end))
    unitModel = combine(first, groupby(unitJoinedData, unitVar)) #by(unitJoinedData, unitVar, first)
    # if !isempty(metadata)
    #     if !isnothing(groupVar) && !isnothing(confounds)
    #         unitModel = unitModel[:, [unitVar, metadata..., groupVar, confounds...]]
    #     elseif !isnothing(groupVar)
    #         unitModel = unitModel[:, [unitVar, metadata..., groupVar]]
    #     elseif !isnothing(confounds)
    #         unitModel = unitModel[:, [unitVar, metadata..., confounds...]]
    #     else
    #         unitModel = unitModel[:, [unitVar, metadata...]]
    #     end
    # else
    #     if !isnothing(groupVar) && !isnothing(confounds)
    #         unitModel = unitModel[:, [unitVar, groupVar, confounds...]]
    #     elseif !isnothing(groupVar)
    #         unitModel = unitModel[:, [unitVar, groupVar]]
    #     elseif !isnothing(confounds)
    #         unitModel = unitModel[:, [unitVar, confounds...]]
    #     else
    #         unitModel = unitModel[:, [unitVar]]
    #     end
    # end

    unitModel = hcat(unitModel, DataFrame(Dict(r => Real[0 for i in 1:nrow(unitModel)]
                                               for r in keys(relationships))))

    unitModel = hcat(unitModel, DataFrame(dim_x=Real[0 for i in 1:nrow(unitModel)], # real sum of my counts times relationship weights
                                          dim_y=Real[0 for i in 1:nrow(unitModel)], 
                                          fit_x=Real[0 for i in 1:nrow(unitModel)], # fitted sum of the above
                                          fit_y=Real[0 for i in 1:nrow(unitModel)])) 
    
    ## Network model
    networkModel = DataFrame(relationship=collect(keys(relationships)),
                             density=Real[0 for r in relationships], # how thick to make the line
                             weight_x=Real[0 for r in relationships], # the weight I contribute to dim_x's
                             weight_y=Real[0 for r in relationships]) # the weight I contribute to dim_y's
    
    ## Code model
    codeModel = DataFrame(code=codes,
                          density=Real[0 for c in codes], # how thick to make the dot
                          fit_x=Real[0 for c in codes], # where to plot this code on the fitted plot's x-axis
                          fit_y=Real[0 for c in codes]) # where to plot this code on the fitted plot's y-axis

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

    ## Normalize and overwrite the model's placeholders
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

    ## Simplify the unitModel down to just those in the treatment/control when groupVar is present
    # if !isnothing(groupVar)
    #     filter!(unitModel) do unitRow
    #         if unitRow[groupVar] in [treatmentGroup, controlGroup]
    #             return true
    #         else
    #             return false
    #         end
    #     end
    # end

    ## TODO drop relationships and units that have all zeros

    # Clean up zeros
    ## Remove unit rows with all zeros
    filter!(unitModel) do unitRow
        if all(unitRow[networkRow[:relationship]] == 0 for networkRow in eachrow(networkModel))
            return false
        else
            return true
        end
    end

    ## Remove unit cols with all zeros
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        if all(unitRow[r] == 0 for unitRow in eachrow(unitModel))
            delete!(unitModel, r)
        end
    end

    ## Remove network rows that no longer have a column in the unit model
    filter!(networkModel) do networkRow
        r = networkRow[:relationship]
        if !(r in names(unitModel))
            return false
        else
            return true
        end
    end

    ## Mean center each relationship
    # for networkRow in eachrow(networkModel)
    #     r = networkRow[:relationship]
    #     mu = mean(unitModel[!, r])
    #     unitModel[!, r] = unitModel[!, r] .- mu
    # end

    # Rotation step
    ## Prepare the config
    # config = Dict{Symbol,Any}()
    # if !isnothing(groupVar)
    #     config[:groupVar] = groupVar
    #     config[:treatmentGroup] = treatmentGroup
    #     config[:controlGroup] = controlGroup
    # end

    # if !isnothing(confounds)
    #     config[:confounds] = confounds
    # end

    ## Use the given lambda, probably one of the ENARotations functions
    # rotateBy(networkModel, unitModel, config)
    rotateBy(networkModel, unitModel)

    ## Compute dim_x and dim_y for the units, now that we have the rotation
    for unitRow in eachrow(unitModel)
        unitRow[:dim_x] = sum(
            networkRow[:weight_x] * unitRow[networkRow[:relationship]]
            for networkRow in eachrow(networkModel)
        )

        unitRow[:dim_y] = sum(
            networkRow[:weight_y] * unitRow[networkRow[:relationship]]
            for networkRow in eachrow(networkModel)
        )
    end

    # Layout step
    ## Compute the density of each network row line
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        networkRow[:density] = sum(unitModel[!, r])
    end

    ## Normalize
    s = maximum(networkModel[!, :density])
    networkModel[!, :density] /= s

    ## Compute the density of each code row dot, by "splitting" the density of each line between its two codes
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        i, j = relationships[r]
        codeModel[i, :density] += networkRow[:density] # TODO check, is this the right way to do this?
        codeModel[j, :density] += networkRow[:density]
    end

    ## Normalize
    s = maximum(codeModel[!, :density])
    codeModel[!, :density] /= s

    ## Place the code dots into fit_x, fit_y, by predicting unit dim_x and dim_y values
    ## This is a projection of the relationship space into the code space
    X = Matrix{Float64}(zeros(nrow(unitModel), 1 + nrow(codeModel)))
    X[:, 1] .= 1 # intercept
    for (i, unitRow) in enumerate(eachrow(unitModel))
        # denom = 2
        denom = 2 * sum(unitRow[t] for t in keys(relationships))
        for r in keys(relationships)
            a, b = relationships[r]
            X[i, a+1] += unitRow[r] / denom # +1 because we started with the intercept
            X[i, b+1] += unitRow[r] / denom # +1 because we started with the intercept
        end
    end

    ## fit_x
    y = Vector{Float64}(unitModel[:, :dim_x])
    # x_axis_coefs = coef(lm(X, y))
    x_axis_coefs = (transpose(X) * X)^-1 * transpose(X) * y
    for i in 1:nrow(codeModel)
        codeModel[i, :fit_x] = x_axis_coefs[i+1] # +1 because we started with the intercept
    end

    ## fit_y
    y = Vector{Float64}(unitModel[:, :dim_y])
    # y_axis_coefs = coef(lm(X, y))
    y_axis_coefs = (transpose(X) * X)^-1 * transpose(X) * y
    for i in 1:nrow(codeModel)
        codeModel[i, :fit_y] = y_axis_coefs[i+1] # +1 because we started with the intercept
    end

    ## Fit the units
    # TODO verify
    for unitRow in eachrow(unitModel)
        denom = 2 * sum(unitRow[t] for t in keys(relationships))
        if denom != 0
            unitRow[:fit_x] = x_axis_coefs[1] +
                sum(
                    codeModel[relationships[r][1], :fit_x] * unitRow[r] +
                    codeModel[relationships[r][2], :fit_x] * unitRow[r]
                    for r in keys(relationships)
                ) / denom

            unitRow[:fit_y] = y_axis_coefs[1] +
                sum(
                    codeModel[relationships[r][1], :fit_y] * unitRow[r] +
                    codeModel[relationships[r][2], :fit_y] * unitRow[r]
                    for r in keys(relationships)
                ) / denom
        end
    end

    ## Translate code model fits to account for the intercepts
    codeModel[!, :fit_x] = codeModel[!, :fit_x] .+ x_axis_coefs[1]
    codeModel[!, :fit_y] = codeModel[!, :fit_y] .+ y_axis_coefs[1]

    ## Translate everything to account for overall mean?
    mu_x = mean(unitModel[!, :fit_x])
    mu_y = mean(unitModel[!, :fit_y])
    unitModel[!, :dim_x] = unitModel[!, :dim_x] .- mu_x
    unitModel[!, :dim_y] = unitModel[!, :dim_y] .- mu_y
    unitModel[!, :fit_x] = unitModel[!, :fit_x] .- mu_x
    unitModel[!, :fit_y] = unitModel[!, :fit_y] .- mu_y
    codeModel[!, :fit_x] = codeModel[!, :fit_x] .- mu_x
    codeModel[!, :fit_y] = codeModel[!, :fit_y] .- mu_y

    # Run tests for how well the fit models the dims
    fitDiffs = Real[]
    dimDiffs = Real[]
    for (i, unitRowA) in enumerate(eachrow(unitModel))
        for (j, unitRowB) in enumerate(eachrow(unitModel))
            if i < j
                push!(fitDiffs, unitRowA[:fit_x] - unitRowB[:fit_x])
                push!(fitDiffs, unitRowA[:fit_y] - unitRowB[:fit_y])
                push!(dimDiffs, unitRowA[:dim_x] - unitRowB[:dim_x])
                push!(dimDiffs, unitRowA[:dim_y] - unitRowB[:dim_y])
            end
        end
    end

    p = pvalue(OneSampleTTest(fitDiffs, dimDiffs))
    # pvalue(EqualVarianceTTest(x, y))
    # pvalue(UnequalVarianceTTest(x, y))
    # pvalue(MannWhitneyUTest(x, y))
    # pvalue(SignedRankTest(x, y))

    # Done!
    # return ENAModel(unitJoinedData, codes, conversations, units, windowSize, groupVar, treatmentGroup, controlGroup,
    #                 confounds, metadata, rotateBy, relationships, unitModel, networkModel, codeModel, p)
    return ENAModel(unitJoinedData, codes, conversations, units, windowSize, rotateBy, unitModel, networkModel, codeModel, relationships, p)
end