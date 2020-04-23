struct ENAModel
    data::DataFrame # the raw data
    codes::Array{Symbol,1}
    conversations::Array{Symbol,1}
    unitVar::Symbol
    windowSize::Int
    groupVar::Union{Nothing,Symbol}
    treatmentGroup::Any # if groupVar == treatmentGroup, then 1
    controlGroup::Any # if groupVar == treatmentGroup, then 0; if neither, leave it out
    confounds::Array{Symbol,1}
    metadata::Array{Symbol,1}
    # plots::Dict{String,Scene} # dropping this in favor of doing it with display/plot methods, for now maybe
    rotateBy::Function
    relationships::Array{Symbol,1}
    unitModel::DataFrame # all the unit-level data we compute
    networkModel::DataFrame # all the connections-level data we compute
    codeModel::DataFrame # all the code-level data we compute
end

function ENAModel(data::DataFrame, codes::Array{Symbol,1}, conversations::Array{Symbol,1}, unitVar::Symbol,
    metadata::Array{Symbol,1}=Symbol[], windowSize::Int=4, confounds::Union{Array{Symbol,1},Nothing}=nothing,
    groupVar::Union{Symbol,Nothing}=nothing, treatmentGroup::Any=nothing, controlGroup::Any=nothing,
    rotateBy::Function=svd_rotation!)

    # Preparing model structures
    ## Relationships between codes
    relationships = Dict(Symbol(string(code1, "_", code2)) => (i, j)
                         for (i, code1) in enumerate(codes)
                         for (j, code2) in enumerate(codes)
                         if i < j)

    ## Unit model
    unitModel = DataFrame()
    if !isnothing(groupVar) && !isnothing(confounds)
        unitModel = by(data, unitVar,
                       [m=>first for m in metadata]...,
                       groupVar=>first,
                       [c=>first for c in confounds]...)
    elseif !isnothing(groupVar)
        unitModel = by(data, unitVar,
                       [m=>first for m in metadata]...,
                       groupVar=>first)
    elseif !isnothing(confounds)
        unitModel = by(data, unitVar,
                       [m=>first for m in metadata]...,
                       [c=>first for c in confounds]...)
    else
        unitModel = by(data, unitVar,
                       [m=>first for m in metadata]...)
    end

    unitModel = hcat(unitModel, DataFrame(Dict(r => Real[0 for i in 1:nrow(unitModel)]
                                               for r in keys(relationships))))

    unitModel = hcat(unitModel, DataFrame(dim_x=Real[0 for i in 1:nrow(unitModel)], # real sum of my counts times relationship weights
                                          dim_y=Real[0 for i in 1:nrow(unitModel)], 
                                          fit_x=Real[0 for i in 1:nrow(unitModel)], # fitted sum of the above
                                          fit_y=Real[0 for i in 1:nrow(unitModel)])) 
    
    ## Network model
    networkModel = DataFrame(relationship=keys(relationships),
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

        unit = line[unitVar]
        for r in keys(relationships)
            i, j = relationships[r]
            if howrecents[i] == 0 && howrecents[j] <= windowSize
                counts[unit][i][j] += 1
            elseif howrecents[j] == 0 && howrecents[i] <= windowSize
                counts[unit][i][j] += 1
            end
        end
    end

    ## Normalize and overwrite the model's placeholders
    for unitRow in eachrow(unitModel)
        unit = unitRow[unitVar]
        vector = [counts[unit][i][j] for (i,j) in values(relationships)]
        s = std(vector)
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
end

function svd_rotation!(networkModel, unitModel, config)
    # TODO compute the thickness and weights of the network-level model
    if haskey(config, :confounds)
        # TODO use AC-PCA when confounds present
    else
        # TODO use PCA otherwise
    end
end

function means_rotation!(networkModel, unitModel, config)
    # TODO compute the thickness and weights of the network-level model
    if haskey(config, :confounds) && haskey(config, :groupVar)
        # TODO use moderated MR1 when confounds present
    elseif haskey(config, :groupVar)
        # TODO use default MR1 otherwise (this can probably be generalized into the above)
    else
        error("means_rotation requires a groupVar")
    end
end