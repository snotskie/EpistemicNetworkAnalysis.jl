struct ENAModel
    data::DataFrame # the raw data
    codes::Array{Symbol,1}
    conversations::Array{Symbol,1}
    units::Array{Symbol,1}
    windowSize::Int
    groupVar::Union{Nothing,Symbol}
    treatmentGroup::Any # if groupVar == treatmentGroup, then 1, else 0
    confounds::Array{Symbol,1}
    metadata::Array{Symbol,1}
    # plots::Dict{String,Scene} # dropping this in favor of doing it with display/plot methods, for now maybe
    rotateBy::Function
    relationships::Array{Symbol,1}
    unitModel::DataFrame # all the unit-level data we compute
    networkModel::DataFrame # all the connections-level data we compute
    codeModel::DataFrame # all the code-level data we compute
end

function ENAModel(data::DataFrame, codes::Array{Symbol,1}, conversations::Array{Symbol,1}, units::Array{Symbol,1},
    metadata::Array{Symbol,1}=Symbol[], windowSize::Int=4, confounds::Union{Array{Symbol,1},Nothing}=nothing,
    groupVar::Union{Symbol,Nothing}=nothing, treatmentGroup::Any=nothing, rotateBy::Function=svd_rotation)

    # Preparing model structures
    relationships = Dict(Symbol(string(code1, "_", code2)) => (i, j)
                         for (i, code1) in enumerate(codes)
                         for (j, code2) in enumerate(codes)
                         if i < j)

    unitModel = data[[], [units..., conversations..., metadata...]]
    if !isnothing(groupVar)
        unitModel = hcat(unitModel, data[[], groupVar])
    end

    if !isnothing(confounds)
        unitModel = hcat(unitModel, data[[], [confounds...]])
    end

    unitModel = hcat(unitModel, data[[], [codes...]])
    unitModel = hcat(unitModel, DataFrame(Dict(r => Real[]
                                               for r in keys(relationships))))

    unitModel = hcat(unitModel, DataFrame(dim_x=Real[], dim_y=Real[], # real sum of my counts times relationship weights
                                          fit_x=Real[], fit_y=Real[])) # fitted sum of the above
    
    networkModel = DataFrame(relationship=keys(relationships),
                             thickness=Real[0 for r in relationships], # how thick to make the line
                             weight_x=Real[0 for r in relationships], # the weight I contribute to dim_x's
                             weight_y=Real[0 for r in relationships]) # the weight I contribute to dim_y's
    
    codeModel = DataFrame(code=codes,
                          thickness=Real[0 for c in codes], # how thick to make the dot
                          fit_x=Real[0 for c in codes], # where to plot this code on the fitted plot's x-axis
                          fit_y=Real[0 for c in codes]) # where to plot this code on the fitted plot's y-axis

    # Accumulation step
    for convo in groupby(data, [conversations..., units...])
        counts = [[0 for j in codes] for i in codes]
        howrecents = [Inf for i in codes]

        ## Raw counts
        for line in eachrow(convo)
            for (i, code) in enumerate(codes)
                if line[code] > 0
                    howrecents[i] = 0
                else
                    howrecents[i] += 1
                end
            end

            for r in keys(relationships)
                i, j = relationships[r]
                if howrecents[i] == 0 && howrecents[j] <= windowSize
                    counts[i][j] += 1
                elseif howrecents[j] == 0 && howrecents[i] <= windowSize
                    counts[i][j] += 1
                end
            end
        end

        ## Put raw counts into the data frame
        row = convo[1, [units..., conversations..., metadata..., codes...]]
        for r in keys(relationships)
            i, j = relationships[r]
            row = hcat(row, DataFrame(Dict(r => Real[counts[i][j]])))
        end

        ## Sphere normalize the raw counts
        s = std(row[keys(relationships)])
        row[keys(relationships)] = map(row[keys(relationships)]) do raw
            return raw/s
        end

        ## Add the row to the unitModel data frame
        row = hcat(row, DataFrame(dim_x=Real[0], dim_y=Real[0], fit_x=Real[0], fit_y=Real[0]))
        push!(unitModel, row)
    end

    # Rotation step
    config = Dict{Symbol,Any}()
    if !isnothing(groupVar)
        config[:groupVar] = groupVar
        config[:treatmentGroup] = treatmentGroup
    end

    if !isnothing(confounds)
        config[:confounds] = confounds
    end

    networkModel = rotateBy(unitModel, networkModel, config)

    # Layout step
    # TODO fit the x and y positions of the unit-level model and code-level model
    # TODO compute the dot sizes for the code-level model
end

function svd_rotation(unitModel, networkModel, config)
    # TODO compute the thickness and weights of the network-level model
    # TODO use LASSO when confounds present
end

function means_rotation(unitModel, networkModel, config)
    # TODO compute the thickness and weights of the network-level model
    # TODO use moderated MR1 when confounds present
    # HERE
end