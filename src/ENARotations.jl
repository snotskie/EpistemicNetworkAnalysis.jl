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
        for networkRow in eachrow(networkModel)
            r = networkRow[:relationship]
            y = Vector{Float64}(factoredUnitModel[!, r])
            ols = lm(X, y)
            # print(ols) # TEMP
            slope = coef(ols)[2]
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