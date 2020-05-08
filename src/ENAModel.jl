struct ENAModel
    data::DataFrame # the raw data
    codes::Array{Symbol,1}
    conversations::Array{Symbol,1}
    units::Array{Symbol,1}
    windowSize::Int
    groupVar::Union{Nothing,Symbol}
    treatmentGroup::Any # if groupVar == treatmentGroup, then 1
    controlGroup::Any # if groupVar == treatmentGroup, then 0; if neither, leave it out
    confounds::Union{Nothing,Array{Symbol,1}}
    metadata::Array{Symbol,1}
    # plots::Dict{String,Scene} # dropping this in favor of doing it with display/plot methods, for now maybe
    rotateBy::Function
    relationships::Array{Symbol,1}
    unitModel::DataFrame # all the unit-level data we compute
    networkModel::DataFrame # all the connections-level data we compute
    codeModel::DataFrame # all the code-level data we compute
end

function ENAModel(data::DataFrame, codes::Array{Symbol,1}, conversations::Array{Symbol,1}, units::Array{Symbol,1};
    metadata::Array{Symbol,1}=Symbol[], windowSize::Int=4, confounds::Union{Nothing,Array{Symbol,1}}=nothing,
    groupVar::Union{Nothing,Symbol}=nothing, treatmentGroup::Any=nothing, controlGroup::Any=nothing,
    rotateBy::Function=svd_rotation!)

    # Preparing model structures
    ## Relationships between codes
    relationships = Dict(Symbol(string(code1, "_", code2)) => (i, j)
                         for (i, code1) in enumerate(codes)
                         for (j, code2) in enumerate(codes)
                         if i < j)

    ## Unit model
    unitVar = :ENA_UNIT
    unitJoinedData = hcat(data, DataFrame(unitVar => map(eachrow(data)) do dataRow # TODO do this automatically when a group var, simplify the units param down to just a unitVar scalar again?
        return join(dataRow[units], ".")
    end))
    unitModel = by(unitJoinedData, unitVar, first)
    if !isempty(metadata)
        if !isnothing(groupVar) && !isnothing(confounds)
            unitModel = unitModel[:, [unitVar, metadata..., groupVar, confounds...]]
        elseif !isnothing(groupVar)
            unitModel = unitModel[:, [unitVar, metadata..., groupVar]]
        elseif !isnothing(confounds)
            unitModel = unitModel[:, [unitVar, metadata..., confounds...]]
        else
            unitModel = unitModel[:, [unitVar, metadata...]]
        end
    else
        if !isnothing(groupVar) && !isnothing(confounds)
            unitModel = unitModel[:, [unitVar, groupVar, confounds...]]
        elseif !isnothing(groupVar)
            unitModel = unitModel[:, [unitVar, groupVar]]
        elseif !isnothing(confounds)
            unitModel = unitModel[:, [unitVar, confounds...]]
        else
            unitModel = unitModel[:, [unitVar]]
        end
    end

    unitModel = hcat(unitModel, DataFrame(Dict(r => Real[0 for i in 1:nrow(unitModel)]
                                               for r in keys(relationships))))

    unitModel = hcat(unitModel, DataFrame(dim_x=Real[0 for i in 1:nrow(unitModel)], # real sum of my counts times relationship weights
                                          dim_y=Real[0 for i in 1:nrow(unitModel)], 
                                          fit_x=Real[0 for i in 1:nrow(unitModel)], # fitted sum of the above
                                          fit_y=Real[0 for i in 1:nrow(unitModel)])) 
    
    ## Network model
    networkModel = DataFrame(relationship=collect(keys(relationships)),
                             thickness=Real[0 for r in relationships], # how thick to make the line
                             weight_x=Real[0 for r in relationships], # the weight I contribute to dim_x's
                             weight_y=Real[0 for r in relationships]) # the weight I contribute to dim_y's
    
    ## Code model
    codeModel = DataFrame(code=codes,
                          thickness=Real[0 for c in codes], # how thick to make the dot
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
    for unitRow in eachrow(unitModel)
        unit = unitRow[unitVar]
        vector = [counts[unit][i][j] for (i,j) in values(relationships)]
        s = sqrt(sum(vector .^ 2))
        vector /= s
        for (k, r) in enumerate(keys(relationships))
            unitRow[r] = vector[k]
        end
    end

    ## Simplify the unitModel down to just those in the treatment/control when groupVar is present
    if !isnothing(groupVar)
        filter!(unitModel) do unitRow
            if unitRow[groupVar] in [treatmentGroup, controlGroup]
                return true
            else
                return false
            end
        end
    end

    # Rotation step
    ## Prepare the config
    config = Dict{Symbol,Any}()
    if !isnothing(groupVar)
        config[:groupVar] = groupVar
        config[:treatmentGroup] = treatmentGroup
        config[:controlGroup] = controlGroup
    end

    if !isnothing(confounds)
        config[:confounds] = confounds
    end

    ## Use the given lambda, probably one of the ENARotations functions
    rotateBy(networkModel, unitModel, config)

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
    ## Compute the thickness of each network row line
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        networkRow[:thickness] = sum(unitModel[!, r])
    end

    ## Normalize
    s = sqrt(sum(networkModel[!, :thickness] .^ 2))
    networkModel[!, :thickness] /= s

    ## Compute the thickness of each code row dot, by "splitting" the thickness of each line between its two codes
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        i, j = relationships[r]
        codeModel[i, :thickness] += networkRow[:thickness] # TODO check, is this the right way to do this?
        codeModel[j, :thickness] += networkRow[:thickness]
    end

    ## Normalize
    s = sqrt(sum(codeModel[!, :thickness] .^ 2))
    codeModel[!, :thickness] /= s

    ## Place the code dots into fit_x, fit_y, by predicting unit dim_x and dim_y values
    ## This is a projection of the relationship space into the code space
    # TODO verify, if there's a part that's wrong when doing the fitting, I believe it's this and/or the coef's below
    X = Matrix{Float64}(zeros(nrow(unitModel), 1 + nrow(codeModel)))
    X[:, 1] .= 1 # intercept
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        i, j = relationships[r]
        for (k, unitRow) in enumerate(eachrow(unitModel))
            X[k, i+1] += unitRow[r] / 2 # +1 because we started with the intercept
            X[k, j+1] += unitRow[r] / 2 # +1 because we started with the intercept
        end
    end

    ## fit_x
    y = Vector{Float64}(unitModel[:, :dim_x])
    ols = lm(X, y)
    for i in 1:nrow(codeModel)
        codeModel[i, :fit_x] = coef(ols)[i+1] # +1 because we started with the intercept
    end

    ## fit_y
    y = Vector{Float64}(unitModel[:, :dim_y])
    ols = lm(X, y)
    for i in 1:nrow(codeModel)
        codeModel[i, :fit_y] = coef(ols)[i+1] # +1 because we started with the intercept
    end

    ## Fit the units
    # TODO verify, but I'm pretty confident that *this* part is right
    for unitRow in eachrow(unitModel)
        unitRow[:fit_x] = sum(
            codeModel[relationships[networkRow[:relationship]][1], :fit_x] * unitRow[networkRow[:relationship]] / 2 +
            codeModel[relationships[networkRow[:relationship]][2], :fit_x] * unitRow[networkRow[:relationship]] / 2
            for networkRow in eachrow(networkModel) # TODO simplify, here and below
        )

        unitRow[:fit_y] = sum(
            codeModel[relationships[networkRow[:relationship]][1], :fit_y] * unitRow[networkRow[:relationship]] / 2 +
            codeModel[relationships[networkRow[:relationship]][2], :fit_y] * unitRow[networkRow[:relationship]] / 2
            for networkRow in eachrow(networkModel)
        )
    end

    # Done!
    return ENAModel(unitJoinedData, codes, conversations, units, windowSize, groupVar, treatmentGroup, controlGroup,
                    confounds, metadata, rotateBy, collect(keys(relationships)), unitModel, networkModel, codeModel)
end