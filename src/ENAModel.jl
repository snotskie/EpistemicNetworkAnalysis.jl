struct ENAModel

    data::DataFrame # the raw data
    codes::Array{Symbol,1}
    conversations::Array{Symbol,1}
    units::Array{Symbol,1}
    windowSize::Int
    groupVar::Union{Nothing,Symbol}
    groupVals::Array{Any,1}
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
    metadata::Array{Symbol,1}=Symbol[], windowSize::Int=4)

    # Preparing model structures
    relationships = Dict(Symbol(string(code1, "_", code2)) => (i, j)
                         for (i, code1) in enumerate(codes)
                         for (j, code2) in enumerate(codes)
                         if i < j)

    # TODO handle groupvar when it's present
    unitModel = data[[], [units..., conversations..., metadata..., codes...]]
    unitModel = hcat(unitModel, DataFrame(Dict(relationship => Int[]
                                               for relationship in keys(relationships))))

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
    # TODO compute the unit-level model
    for convo in groupby(data, [conversations..., units...])
        counts = [[0 for j in codes] for i in codes]
        howrecents = [Inf for i in codes]
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

        row = convo[1, [units..., conversations..., metadata..., codes...]]
        for r in keys(relationships)
            i, j = relationships[r]
            row = hcat(row, DataFrame(Dict(relationship => Int[counts[i][j]])))
        end

        # TODO sphere normalize

        row = hcat(row, DataFrame(dim_x=Real[0], dim_y=Real[0], fit_x=Real[0], fit_y=Real[0]))
        push!(unitModel, row)
    end

    # Rotation step
    # TODO compute the thickness and weights of the network-level model
    # TODO set default rotateBy in different cases
    networkModel = rotateBy(unitModel, networkModel) # TODO pass in group var, group vals, and confounds

    # Layout step
    # TODO fit the x and y positions of the unit-level model and code-level model
    # TODO compute the dot sizes for the code-level model
end