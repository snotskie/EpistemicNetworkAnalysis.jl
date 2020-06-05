# TODO borrow a rotation from another ena object

"""
TODO: document
"""
abstract type ENARotation
end

"""
SVD rotation finds the x- and y-axes that explain the most
variance of the network model while being orthogonal to one
another. When confounds are given, these axes are penalized
for colinearity with those confounds, attempting to optimize
for orthogonality with the first listed confound. No orthogonality
with the confounds is guaranteed, just approximated. This approximation
should be good enough to interpret the quadrants.
"""
struct SVDRotation <: ENARotation
end

"""
Means rotation is a special case of regression rotation. First,
means rotation factors the grouping variable into a dummy 0/1,
then it runs a regression rotation with that dummy as the first
listed confound, followed by the given confounds (if any), followed
by the interactions between the group dummy and each of those confounds (if any).
This finds the x-axis of the network space such that it captures
the variance explained by the effect of the grouping variable, moderated
by the confounds. The y-axis is then set to capture the most variance while being
orthogonal to the x-axis. This y-axis could very well be colinear
with one or more of the other confounds given.
"""
struct MeansRotation <: ENARotation
    groupVar::Symbol
    controlGroup::Any
    treatmentGroup::Any
end

"""
TODO: code and document
"""
struct ModeratedMeansRotation <: ENARotation
    groupVar::Symbol
    controlGroup::Any
    treatmentGroup::Any
    confounds::Array{Symbol,1}
end

"""
Regression rotation finds the x-axis of the network space
such that it captures the variance explained by the effect
of predictor, holding the given confounds constant.
The y-axis is then set to capture the most variance while being
orthogonal to the x-axis. This y-axis could very well be colinear
with one or more of the other confounds given.
"""
struct RegressionRotation <: ENARotation
    predictor::Symbol
    confounds::Array{Symbol,1}
end

"""
TODO: document
"""
struct TwoGroupRotation <: ENARotation
    groupVar1::Symbol
    controlGroup1::Any
    treatmentGroup1::Any
    groupVar2::Symbol
    controlGroup2::Any
    treatmentGroup2::Any
    confounds::Array{Symbol,1}
end

# SVD Rotation
function (rotation::SVDRotation)(networkModel, unitModel)
    pcaModel = help_deflating_svd(networkModel, unitModel)
    networkModel[!, :weight_x] = pcaModel[:, 1]
    networkModel[!, :weight_y] = pcaModel[:, 2]
end

# Means Rotation
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

# Moderated Means Rotation
function (rotation::ModeratedMeansRotation)(networkModel, unitModel)
    ## TODO
end

# Regression Rotation
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
    X = (transpose(X) * X)^-1 * transpose(X)

    ## For each relationship, find the effect of the first listed column
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        y = Vector{Float64}(filteredUnitModel[!, r])
        coefs = X * y
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

# Two Group Rotation
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

    ## Find the line of the first group's effect through the high dimensional space
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        y = Vector{Float64}(factoredUnitModel[!, r])
        coefs = X * y
        networkRow[:weight_x] = coefs[2]
    end

    ## Normalize the weights
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    ## Find the values for each unit projected onto the x-axis (for the whole unit model?)
    xAxisColumn =
        Matrix{Float64}(factoredUnitModel[!, networkModel[!, :relationship]]) *
        Vector{Float64}(networkModel[!, :weight_x])

    ## Find the line of the second group's effect through the high dimensional space,
    ## orthogonal to (independent from) the first
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        y = Vector{Float64}(factoredUnitModel[!, r])
        scalar = dot(y, xAxisColumn) / dot(xAxisColumn, xAxisColumn)
        y -= scalar * xAxisColumn
        coefs = X * y
        networkRow[:weight_y] = coefs[3]
    end

    ## Normalize the weights
    s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
    if s != 0
        networkModel[!, :weight_y] /= s
    end
end

#### OLD STUFF BELOW ####


# ## TODO refactor so that rotations take arguments and produce a rotationg function

# # TODO LASSO
# # LASSO rotation finds the x-axis of the network model that
# # best predicts the first listed confound, holding the other
# # confounds constant, and with penalties for weak correlation.
# # This should result in a best attempt at a model that would pass
# # k-fold cross validation.

# # TODO xy LASSO

# # TODO verify this


# # TODO test


