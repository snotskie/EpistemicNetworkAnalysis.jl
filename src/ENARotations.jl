abstract type ENARotation
end

struct SVDRotation <: ENARotation
end

struct MeansRotation <: ENARotation
    groupVar::Symbol
    controlGroup::Any
    treatmentGroup::Any
end

struct ModeratedMeansRotation <: ENARotation
    groupVar::Symbol
    controlGroup::Any
    treatmentGroup::Any
    confounds::Array{Symbol,1}
end

struct RegressionRotation <: ENARotation
    predictor::Symbol
    confounds::Array{Symbol,1}
end

struct TwoGroupRotation <: ENARotation
    groupVar1::Symbol
    controlGroup1::Any
    treatmentGroup1::Any
    groupVar2::Symbol
    controlGroup2::Any
    treatmentGroup2::Any
    confounds::Array{Symbol,1}
end

# SVD Rotation: TODO document
function (rotation::SVDRotation)(networkModel, unitModel)
    pcaModel = help_deflating_svd(networkModel, unitModel)
    networkModel[!, :weight_x] = pcaModel[:, 1]
    networkModel[!, :weight_y] = pcaModel[:, 2]
end

# Means Rotation: TODO document
function (rotation::MeansRotation)(networkModel, unitModel)
    ## Filter the unit model to just those in the control/treatment group
    filteredUnitModel = filter(unitModel) do unitRow
        if unitRow[rotation.groupVar] == rotation.controlGroup
            return true
        elseif unitRow[rotation.groupVar] == rotation.treatmentGroup
            return true
        else
            return false
        end
    end

    ## Factor the control/treatment label to just 0/1
    factors = map(eachrow(filteredUnitModel)) do filteredRow
        if filteredRow[rotation.groupVar] == rotation.controlGroup
            return 0.0
        else
            return 1.0
        end
    end

    factoredUnitModel = hcat(filteredUnitModel, DataFrame(:factoredGroupVar => factors))

    ## Use a RegressionRotation to get the rest of the job done
    rr = RegressionRotation(:factoredGroupVar, [])
    rr(networkModel, factoredUnitModel)
end

# Moderated Means Rotation: TODO document
function (rotation::ModeratedMeansRotation)(networkModel, unitModel)
    ## TODO
end

# Regression Rotation: TODO document
function (rotation::RegressionRotation)(networkModel, unitModel)
    ## Collect the columns we need
    columns = [rotation.predictor, rotation.confounds...]

    ## Filter rows to exclude those with missing data
    goodRows = completecases(unitModel[!, columns])
    filteredUnitModel = unitModel[goodRows, :]

    ## Predictors for the linear model.
    X = DataFrame(:Intercept => [1 for i in 1:nrow(filteredUnitModel)])
    X = hcat(X, filteredUnitModel[!, columns])
    X = Matrix{Float64}(X)

    ## For each relationship, find the effect of the first listed column
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        y = Vector{Float64}(filteredUnitModel[!, r])
        # ols = lm(X, y)
        # slope = coef(ols)[2]
        coefs = (transpose(X) * X)^-1 * transpose(X) * y
        slope = coefs[2]
        networkRow[:weight_x] = slope
    end

    ## Normalize the weights
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    ## Find the first svd dim of the data orthogonal to the x weights, use these as the y weights
    unitValues = Matrix{Float64}(filteredUnitModel[!, networkModel[!, :relationship]])
    # coefValues = transpose(Matrix{Float64}(DataFrame(hcat(networkModel[!, :coefs]...))))
    # controlModel = DataFrame(unitValues * coefValues)
    axisValues = Matrix{Float64}(DataFrame(:weight_x => networkModel[!, :weight_x]))
    controlModel = DataFrame(unitValues * axisValues)
    pcaModel = help_deflating_svd(networkModel, filteredUnitModel, controlModel)
    networkModel[!, :weight_y] = pcaModel[:, 1]
end

# Two Group Rotation: TODO document
function (rotation::TwoGroupRotation)(networkModel, unitModel)
    ## Filter rows to exclude those with missing data
    filteredUnitModel = filter(unitModel) do unitRow
        if unitRow[rotation.groupVar1] == rotation.controlGroup1
            return true
        elseif unitRow[rotation.groupVar1] == rotation.treatmentGroup1
            return true
        else
            return false
        end
    end

    filteredUnitModel = filter(filteredUnitModel) do unitRow
        if unitRow[rotation.groupVar2] == rotation.controlGroup2
            return true
        elseif unitRow[rotation.groupVar2] == rotation.treatmentGroup2
            return true
        else
            return false
        end
    end

    ## When we have confounds, drop rows with missing confound data
    if !isempty(rotation.confounds)
        goodRows = completecases(filteredUnitModel[!, rotation.confounds])
        filteredUnitModel = filteredUnitModel[goodRows, :]
    end

    ## Factor the control/treatment label to just 0/1
    factors1 = map(eachrow(filteredUnitModel)) do filteredRow
        if filteredRow[rotation.groupVar1] == rotation.controlGroup1
            return 0.0
        else
            return 1.0
        end
    end

    factors2 = map(eachrow(filteredUnitModel)) do filteredRow
        if filteredRow[rotation.groupVar2] == rotation.controlGroup2
            return 0.0
        else
            return 1.0
        end
    end

    factoredUnitModel = hcat(filteredUnitModel,
        DataFrame(:factoredGroupVar1 => factors1, :factoredGroupVar2 => factors2))
    
    ## Interact the two groups, mean centered
    factoredUnitModel[!, :factoredGroupInteraction] =
        (factoredUnitModel[!, :factoredGroupVar1] .- mean(factoredUnitModel[!, :factoredGroupVar1])) .*
        (factoredUnitModel[!, :factoredGroupVar2] .- mean(factoredUnitModel[!, :factoredGroupVar2]))

    ## Collect the columns we need
    columns = [:factoredGroupVar1, :factoredGroupVar2, :factoredGroupInteraction, rotation.confounds...]

    ## Predictors for the linear model.
    X = DataFrame(:Intercept => [1 for i in 1:nrow(factoredUnitModel)])
    X = hcat(X, factoredUnitModel[!, columns])
    X = Matrix{Float64}(X)
    X = (transpose(X) * X)^-1 * transpose(X)
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        y = Vector{Float64}(factoredUnitModel[!, r])
        coefs = X * y
        networkRow[:weight_x] = coefs[2]
        networkRow[:weight_y] = coefs[3]
    end

    ## Normalize the weights
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
    if s != 0
        networkModel[!, :weight_y] /= s
    end
end

#### OLD STUFF BELOW ####


## TODO refactor so that rotations take arguments and produce a rotationg function

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
        pcaModel = help_deflating_svd(networkModel, filteredUnitModel, controlModel)
        networkModel[!, :weight_x] = pcaModel[:, 1]
        networkModel[!, :weight_y] = pcaModel[:, 2]
    else
        pcaModel = help_deflating_svd(networkModel, unitModel)
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
        # networkModel[!, :coefs] = [Real[0, 0] for networkRow in eachrow(networkModel)]
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
        if s != 0
            networkModel[!, :weight_x] /= s
        end

        ## Find the first svd dim of the data orthogonal to the x weights, use these as the y weights
        unitValues = Matrix{Float64}(filteredUnitModel[!, networkModel[!, :relationship]])
        # coefValues = transpose(Matrix{Float64}(DataFrame(hcat(networkModel[!, :coefs]...))))
        # controlModel = DataFrame(unitValues * coefValues)
        axisValues = Matrix{Float64}(DataFrame(:weight_x => networkModel[!, :weight_x]))
        controlModel = DataFrame(unitValues * axisValues)
        pcaModel = help_deflating_svd(networkModel, filteredUnitModel, controlModel)
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

# IN DEV: put G_1 effect on x and G_2 effect on y, make sure they are orthogonal
function two_group_rotation!(networkModel, unitModel, config)
    ## Must have group variable to use (m)mr1
    if haskey(config, :groupVar)

        ## Must have confounds, using the first confound as the G_2
        if haskey(config, :confounds)

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

            ## Drop rows with missing confound data
            goodRows = completecases(filteredUnitModel[!, config[:confounds]])
            filteredUnitModel = filteredUnitModel[goodRows, :]

            ## Convert the control/treatment label to just 0/1
            factors = map(eachrow(filteredUnitModel)) do filteredRow
                if filteredRow[config[:groupVar]] == config[:controlGroup]
                    return 0.0
                else
                    return 1.0
                end
            end

            factoredUnitModel = hcat(filteredUnitModel, DataFrame(:factoredGroupVar1 => factors))

            ## Convert the first confound to just 0/1
            controlConfound = unique(sort(factoredUnitModel[!, config[:confounds][1]]))[1]
            factors2 = map(eachrow(filteredUnitModel)) do filteredRow
                if filteredRow[config[:confounds][1]] == controlConfound
                    return 0.0
                else
                    return 1.0
                end
            end

            # println(controlConfound)

            factoredUnitModel = hcat(factoredUnitModel, DataFrame(:factoredGroupVar2 => factors2))

            # ## Interact the two groups
            # factoredUnitModel[!, :factoredGroupInteraction] =
            #     factoredUnitModel[!, :factoredGroupVar1] .* factoredUnitModel[!, :factoredGroupVar2]
            
            ## Interact the two groups, mean centered
            factoredUnitModel[!, :factoredGroupInteraction] =
                (factoredUnitModel[!, :factoredGroupVar1] .- mean(factoredUnitModel[!, :factoredGroupVar1])) .*
                (factoredUnitModel[!, :factoredGroupVar2] .- mean(factoredUnitModel[!, :factoredGroupVar2]))

            # println(factoredUnitModel[!, [:factoredGroupVar1, :factoredGroupVar2, :factoredGroupInteraction]])

            ## Reconfigure the confounds
            reconfounds = [:factoredGroupVar1, :factoredGroupVar2, :factoredGroupInteraction]
            config[:confounds] = reconfounds

            ## Predictors for the linear model.
            X = DataFrame(:Intercept => [1 for i in 1:nrow(factoredUnitModel)])
            X = hcat(X, factoredUnitModel[!, [config[:confounds]...]])
            X = Matrix{Float64}(X)

            # println(X)

            # println((transpose(X) * X)^-1)

            # println((transpose(X) * X)^-1 * transpose(X))

            X = (transpose(X) * X)^-1 * transpose(X)

            # println(any(isinf.(X[:])))
            # println(any(isnan.(X[:])))
            # println(any(ismissing.(X[:])))

            # println(X * Vector{Float64}([rand() for i in eachrow(factoredUnitModel)]))

            ## For each relationship, find the effect of the first listed confound
            # networkModel[!, :coefs] = [Real[0, 0] for networkRow in eachrow(networkModel)]
            for networkRow in eachrow(networkModel)
                r = networkRow[:relationship]
                y = Vector{Float64}(factoredUnitModel[!, r])

                # println(any(isinf.(y)))
                # println(any(isnan.(y)))
                # println(any(ismissing.(y)))
                
                # println(r)
                # println(y)
                # ols = lm(X, y)
                # slope = coef(ols)[2]
                # coefs = (transpose(X) * X)^-1 * transpose(X) * y
                coefs = X * y
                # println(coefs)
                # slope = coefs[2]
                # networkRow[:weight_x] = slope
                # networkRow[:coefs] = coefs[2:end]
                networkRow[:weight_x] = coefs[2]
                networkRow[:weight_y] = coefs[3]
            end

            ## Normalize the weights
            s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
            if s != 0
                networkModel[!, :weight_x] /= s
            end

            s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
            if s != 0
                networkModel[!, :weight_y] /= s
            end
        else
            error("this requires a confound")
        end
    else
        error("means_rotation requires a groupVar")
    end
end