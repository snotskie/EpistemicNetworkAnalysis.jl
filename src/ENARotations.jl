# TODO LASSO

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
        networkModel[!, :coefs] = [Real[0, 0] for networkRow in eachrow(networkModel)]
        for networkRow in eachrow(networkModel)
            r = networkRow[:relationship]
            y = Vector{Float64}(factoredUnitModel[!, r])
            # ols = lm(X, y)
            # slope = coef(ols)[2]
            coefs = (transpose(X) * X)^-1 * transpose(X) * y
            slope = coefs[2]
            networkRow[:weight_x] = slope
            networkRow[:coefs] = coefs[2:end]
        end

        ## Normalize the differences of the means
        s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
        networkModel[!, :weight_x] /= s

        ## Find the first svd dim of the data orthogonal to the x weights, use these as the y weights
        unitValues = Matrix{Float64}(factoredUnitModel[!, networkModel[!, :relationship]])
        coefValues = transpose(Matrix{Float64}(DataFrame(hcat(networkModel[!, :coefs]...))))
        controlModel = DataFrame(unitValues * coefValues)
        pcaModel = help_ac_svd(networkModel, factoredUnitModel, controlModel) # TODO optimize the lambda, don't just guess at 10 all the time
        networkModel[!, :weight_y] = pcaModel[:, 1]

        ## Warnings for colinear axes
        for unitRow in eachrow(factoredUnitModel)
            unitRow[:dim_x] = sum(
                networkRow[:weight_x] * unitRow[networkRow[:relationship]]
                for networkRow in eachrow(networkModel)
            )
    
            unitRow[:dim_y] = sum(
                networkRow[:weight_y] * unitRow[networkRow[:relationship]]
                for networkRow in eachrow(networkModel)
            )
        end

        norm_y = sqrt(sum(factoredUnitModel[!, :dim_y] .* factoredUnitModel[!, :dim_y]))
        for i in 1:ncol(controlModel)
            norm_x = sqrt(sum(controlModel[!, i] .* controlModel[!, i]))
            cos_angle = sum(factoredUnitModel[!, :dim_y] .* controlModel[!, i]) / norm_y / norm_x
            if abs(cos_angle) >= 0.15
                @warn """The y-axis appears colinear with dimension #$(i) of the x-axis.
Dimension #1 is the regular x-axis, aka, the (moderated) effect of the grouping variable on the connection weights.
Dimension #2 is the effect from the first confound, #3 is the effect from the second confound, and so on.
The next dimension is the effect from the interacton between the grouping variable and the first confound, and so on.
"""
            end
        end
    else
        error("means_rotation requires a groupVar")
    end
end