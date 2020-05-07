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
        # s = std(vector, corrected=false)
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
    config = Dict{Symbol,Any}()
    if !isnothing(groupVar)
        config[:groupVar] = groupVar
        config[:treatmentGroup] = treatmentGroup
        config[:controlGroup] = controlGroup
    end

    if !isnothing(confounds)
        config[:confounds] = confounds
    end

    rotateBy(networkModel, unitModel, config)

    # Layout step
    # TODO compute dim_x and dim_y of the unit-level model
    # TODO fit the x and y positions of the unit-level model and code-level model
    # TODO compute the dot sizes for the code-level model

    return ENAModel(unitJoinedData, codes, conversations, units, windowSize, groupVar, treatmentGroup, controlGroup,
                    confounds, metadata, rotateBy, collect(keys(relationships)), unitModel, networkModel, codeModel)
end

# TODO verify this
function svd_rotation!(networkModel, unitModel, config)
    if haskey(config, :confounds)
        goodRows = completecases(unitModel[!, config[:confounds]])
        filteredUnitModel = unitModel[goodRows, :]
        controlModel = filteredUnitModel[!, config[:confounds]]
        pcaModel = help_ac_svd(networkModel, filteredUnitModel, controlModel)
        networkModel[!, :weight_x] = pcaModel[:, 1]
        networkModel[!, :weight_y] = pcaModel[:, 2]
    else
        pcaModel = help_ac_svd(networkModel, unitModel)
        networkModel[!, :weight_x] = pcaModel[:, 1]
        networkModel[!, :weight_y] = pcaModel[:, 2]
    end

    # TODO compute the thickness of the relationships
end

function means_rotation!(networkModel, unitModel, config)
    ## Must have group variable to use (m)mr1
    if haskey(config, :groupVar)

        ## Tidy up the unit model to just those in the control/treatment group
        filteredUnitModel = filter(unitModel) do unitRow
            if unitRow[config[:groupVar]] == config[:controlGroup]
                return true
            elseif unitRow[config[:groupVar]] == config[:treatmentGroup]
                return true
            else
                return false
            end
        end

        ## When confounds present, drop rows with missing confound data
        if haskey(config, :confounds)
            goodRows = completecases(filteredUnitModel[!, config[:confounds]])
            filteredUnitModel = filteredUnitModel[goodRows, :]
        end

        ## Convert the control/treatment label to just 0/1
        factors = map(eachrow(filteredUnitModel)) do filteredRow
            if filteredRow[config[:groupVar]] == config[:controlGroup]
                return 0.0
            else
                return 1.0
            end
        end

        factoredUnitModel = hcat(filteredUnitModel, DataFrame(:factoredGroupVar => factors))

        ## Predictors for the linear model. When confounds present, moderate the model.
        X = DataFrame(:Intercept => [1 for i in 1:nrow(factoredUnitModel)])
        X = hcat(X, factoredUnitModel[!, :factoredGroupVar])
        if haskey(config, :confounds)
            # TODO think of a more consistent way to specify columns. right now, confounds must be simple numeric, but groups can be whatever data type, factored into a 0, 1, or drop
            interactions = factoredUnitModel[!, [config[:confounds]...]] .* factoredUnitModel[!, :factoredGroupVar]
            X = hcat(X, factoredUnitModel[!, [config[:confounds]...]])
            X = hcat(X, interactions, makeunique=true)
        end

        X = Matrix{Float64}(X)

        ## For each relationship, find the (moderated) difference of means
        for networkRow in eachrow(networkModel)
            r = networkRow[:relationship]
            y = Vector{Float64}(factoredUnitModel[!, r])
            ols = lm(X, y)
            # print(ols) # TEMP
            slope = coef(ols)[2]
            networkRow[:thickness] = slope # TODO check that this is the right way to do this; negative means blue, positive means red
            networkRow[:weight_x] = slope
        end

        ## Normalize the differences of the means
        s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
        networkModel[!, :weight_x] /= s

        ## Find the first svd dim of the data orthogonal to the x weights, use these as the y weights
        controlModel = DataFrame(:x_axis => [
            sum(
                networkRow[:weight_x] * unitRow[networkRow[:relationship]]
                for networkRow in eachrow(networkModel)
            ) for unitRow in eachrow(factoredUnitModel)
        ])
        pcaModel = help_ac_svd(networkModel, factoredUnitModel, controlModel)

        if haskey(config, :confounds) # TODO: why does this work to match R??
            networkModel[!, :weight_y] = pcaModel[:, 1]
        else
            networkModel[!, :weight_y] = pcaModel[:, 2]
        end
    else
        error("means_rotation requires a groupVar")
    end
end

# weights = Vector{Float64}(networkModel[!, :weight_x]) # CHECKED same as R
# rawCounts = Matrix{Float64}(factoredUnitModel[!, [networkRow[:relationship] for networkRow in eachrow(networkModel)]]) # CHECKED same as R
# meanCenteredCounts = rawCounts .- transpose(collect(mean(rawCounts[:, i]) for i in 1:size(rawCounts)[2])) # TODO say this simpler # CHECKED same as R
# XBar = (meanCenteredCounts - meanCenteredCounts*weights*transpose(weights)) * qr(weights).Q[:, 2:end] # 2:end to remove the x-axis (the 1 col) from the deflation
# display(XBar) # TODO HERE
# pcaModel = fit(PCA, Matrix{Float64}(transpose(XBar)), # TODO check what config of Julia's PCA runs the same algorithm as R's prcomp(X, scale=FALSE)
#     pratio=1.0, mean=0, method=:svd)
# display(projection(pcaModel))
# orthosvd = qr(weights).Q[:, 2:end] * projection(pcaModel) # 2:end to remove the x-axis (the 1 col) from the reinflation
# display(orthosvd)
# networkModel[!, :weight_y] = orthosvd[:, 2]