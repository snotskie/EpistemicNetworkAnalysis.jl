# TODO LASSO
# LASSO rotation finds the x-axis of the network model that
# best predicts the first listed confound, holding the other
# confounds constant, and with penalties for weak correlation.
# This should result in a best attempt at a model that would pass
# k-fold cross validation.

# TODO xy LASSO

# TODO verify this
# SVD rotation finds the x- and y-axes that explain the most
# variance of the network model while being orthogonal to one
# another. When confounds are given, these axes are penalized
# for colinearity with those confounds, attempting to optimize
# for orthogonality with the first listed confound. No orthogonality
# with the confounds is guaranteed, just approximated. This approximation
# should be good enough to interpret the quadrants.
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

# TODO test
# Regression rotation finds the x-axis of the network space
# such that it captures the variance explained by the effect
# of the first listed confound, holding the other confounds constant.
# The y-axis is then set to capture the most variance while being
# orthogonal to the x-axis. This y-axis could very well be colinear
# with one or more of the other confounds given.
function regression_rotation!(networkModel, unitModel, config)
    ## Must have confounds to use rr1
    if haskey(config, :confounds)

        ## Drop rows with missing confound data
        goodRows = completecases(unitModel[!, config[:confounds]])
        filteredUnitModel = unitModel[goodRows, :]

        ## Predictors for the linear model.
        X = DataFrame(:Intercept => [1 for i in 1:nrow(filteredUnitModel)])
        X = hcat(X, filteredUnitModel[!, [config[:confounds]...]])
        X = Matrix{Float64}(X)

        ## For each relationship, find the effect of the first listed confound
        networkModel[!, :coefs] = [Real[0, 0] for networkRow in eachrow(networkModel)]
        for networkRow in eachrow(networkModel)
            r = networkRow[:relationship]
            y = Vector{Float64}(filteredUnitModel[!, r])
            # ols = lm(X, y)
            # slope = coef(ols)[2]
            coefs = (transpose(X) * X)^-1 * transpose(X) * y
            slope = coefs[2]
            networkRow[:weight_x] = slope
            # networkRow[:coefs] = coefs[2:end]
        end

        ## Normalize the weights
        s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
        networkModel[!, :weight_x] /= s

        ## Find the first svd dim of the data orthogonal to the x weights, use these as the y weights
        unitValues = Matrix{Float64}(filteredUnitModel[!, networkModel[!, :relationship]])
        # coefValues = transpose(Matrix{Float64}(DataFrame(hcat(networkModel[!, :coefs]...))))
        # controlModel = DataFrame(unitValues * coefValues)
        axisValues = Matrix{Float64}(DataFrame(:weight_x => networkModel[!, :weight_x]))
        controlModel = DataFrame(unitValues * axisValues)
        pcaModel = help_ac_svd(networkModel, filteredUnitModel, controlModel) # TODO replace this with an ortho svd
        networkModel[!, :weight_y] = pcaModel[:, 1]
    else
        error("regression_rotation requires confounds")
    end
end

# TODO call the above twice, then combine them into an x and y axis; must have at least two confounds
# function xy_regression_rotation!(networkModel, unitModel, config)
# end

# Means rotation is a special case of regression rotation. First,
# means rotation factors the grouping variable into a dummy 0/1,
# then it runs a regression rotation with that dummy as the first
# listed confound, followed by the given confounds (if any), followed
# by the interactions between the group dummy and each of those confounds (if any).
# This finds the x-axis of the network space such that it captures
# the variance explained by the effect of the grouping variable, moderated
# by the confounds. The y-axis is then set to capture the most variance while being
# orthogonal to the x-axis. This y-axis could very well be colinear
# with one or more of the other confounds given.
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

        ## Reconfigure the confounds to pass to the regression_rotation
        ## When confounds were passed here, moderate the regression too.
        reconfounds = [:factoredGroupVar]
        if haskey(config, :confounds)
            reconfounds = [:factoredGroupVar, config[:confounds]...]
            for confound in config[:confounds]
                interaction = Symbol(string(:factoredGroupVar, "_", confound, "_interaction"))
                factoredUnitModel[!, interaction] = factoredUnitModel[!, confound] .* factoredUnitModel[!, :factoredGroupVar]
                push!(reconfounds, interaction)
            end
        end

        config[:confounds] = reconfounds
        return regression_rotation!(networkModel, factoredUnitModel, config)
    else
        error("means_rotation requires a groupVar")
    end
end

# TODO xy_means_rotation!, call xy_regression_rotation; must have group and confounds