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
Formula rotation finds the x-axis of the network space
such that it captures the variance explained by the effect
of the first predictor after the intercept, holding the given
confounds constant. The y-axis is then set to capture the most
variance while being orthogonal to the x-axis. This y-axis could
very well be colinear with one or more of the other confounds given.
"""
struct FormulaRotation{T <: RegressionModel} <: ENARotation
    regression_model::Type{T}
    f1::FormulaTerm
end

## Aliasing
const Formula1Rotation = FormulaRotation

"""
TODO document
"""
struct Formula2Rotation{T <: RegressionModel} <: ENARotation
    regression_model1::Type{T}
    f1::FormulaTerm
    regression_model2::Type{T}
    f2::FormulaTerm
end

# SVD Rotation
function rotate!(::SVDRotation, networkModel::DataFrame, unitModel::DataFrame)
    pcaModel = help_deflating_svd(networkModel, unitModel)
    networkModel[!, :weight_x] = pcaModel[:, 1]
    networkModel[!, :weight_y] = pcaModel[:, 2]
end

# Means Rotation
function rotate!(rotation::MeansRotation, networkModel::DataFrame, unitModel::DataFrame)
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

    ## Use a FormulaRotation to do the rest of the work
    fr = FormulaRotation(LinearModel, @formula(y ~ 1 + factoredGroupVar))
    rotate!(fr, networkModel, factoredUnitModel)
end

# Formula Rotation
function rotate!(rotation::FormulaRotation, networkModel::DataFrame, unitModel::DataFrame)
    ## TODO check assumptions about f1

    ## For each relationship, find the effect of the first predictor after the intercept
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        f1 = FormulaTerm(term(r), rotation.f1.rhs)
        m1 = fit(rotation.regression_model, f1, unitModel)
        slope = coef(m1)[2]
        networkRow[:weight_x] = slope
        println(r)
        display(m1)
        display(coef(m1))
        println()
    end

    ## Normalize the weights
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    ## Find the first svd dim of the data orthogonal to the x weights, use these as the y weights
    unitValues = Matrix{Float64}(unitModel[!, networkModel[!, :relationship]])
    axisValues = Matrix{Float64}(DataFrame(:weight_x => networkModel[!, :weight_x]))
    controlModel = DataFrame(unitValues * axisValues)
    pcaModel = help_deflating_svd(networkModel, unitModel, controlModel)
    networkModel[!, :weight_y] = pcaModel[:, 1]
end

# Formula2 Rotation
function rotate!(rotation::Formula2Rotation, networkModel::DataFrame, unitModel::DataFrame)
    ## TODO check assumptions about f1 and f2

    ## For each relationship, find the effect of the first predictor after the intercept
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        f1 = FormulaTerm(term(r), rotation.f1.rhs)
        m1 = fit(rotation.regression_model1, f1, unitModel)
        slope = coef(m1)[2]
        networkRow[:weight_x] = slope
        println(r)
        display(m1)
        display(coef(m1))
        println()
    end

    ## Normalize the weights
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    ## Orthogonalization
    unitValues = Matrix{Float64}(unitModel[!, networkModel[!, :relationship]])
    axisValues = Matrix{Float64}(DataFrame(:weight_x => networkModel[!, :weight_x]))
    controlModel = DataFrame(unitValues * axisValues)
    x = controlModel[!, 1]
    orthoUnitModel = copy(unitModel)
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        y = orthoUnitModel[:, r]
        scalar = dot(y, x) / dot(x, x)
        y -= scalar * x
        orthoUnitModel[!, r] = y
    end

    ## For each relationship in the ortho model, find the effect of the first predictor after the intercept
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        f2 = FormulaTerm(term(r), rotation.f2.rhs)
        m2 = fit(rotation.regression_model2, f2, orthoUnitModel)
        slope = coef(m2)[2]
        networkRow[:weight_y] = slope
        println(r)
        display(m2)
        display(coef(m2))
        println()
    end

    ## Normalize the weights
    s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
    if s != 0
        networkModel[!, :weight_y] /= s
    end
end


# # Hierarchical Rotation
# function rotate!(rotation::HierarchicalRotation, networkModel::DataFrame, unitModel::DataFrame)
#     ## Collect our columns
#     columns = unique([
#         rotation.xAxisVar,
#         rotation.yAxisYar,
#         rotation.confounds...,
#         (values(rotation.nestings)...)...
#     ])

#     ## Filter out rows with missing data
#     goodRows = completecases(unitModel[!, columns])
#     filteredUnitModel = unitModel[goodRows, :]

#     ## Interact the two axis variables, mean centered, and add to our columns list
#     filteredUnitModel[!, :axisInteraction] =
#         (filteredUnitModel[!, rotation.xAxisVar] .- mean(filteredUnitModel[!, rotation.xAxisVar])) .*
#         (filteredUnitModel[!, rotation.yAxisYar] .- mean(filteredUnitModel[!, rotation.yAxisYar]))

#     push!(columns, :axisInteraction)

#     ## Independent variables for the main linear model (more efficient to compute this just once out here)
#     X = DataFrame(:Intercept => [1 for i in 1:nrow(filteredUnitModel)])
#     X = hcat(X, filteredUnitModel[!, columns])
#     X = Matrix{Float64}(X)
#     X = (transpose(X) * X)^-1 * transpose(X)

#     ## Find the line of the first group's effect through the high dimensional space, controlling for such and such
#     ## We do this once for the x-axis, then once for the y-axis, with small differences
#     xAxisColumn = Vector{Float64}([0]) # Placeholder
#     for T in 1:2
#         for networkRow in eachrow(networkModel)
#             r = networkRow[:relationship]

#             ## Dependent variable
#             y = Vector{Float64}(filteredUnitModel[!, r])

#             ## For the y-axis iteration, make the columns independent of the x-axis we've already found
#             if T == 2
#                 scalar = dot(y, xAxisColumn) / dot(xAxisColumn, xAxisColumn)
#                 y -= scalar * xAxisColumn
#             end

#             ### Adjust for nestings
#             #### TODO FIXME this part seems not to work at all, give it a look
#             toExclude = []
#             for nest in keys(rotation.nestings)
#                 groups = unique(filteredUnitModel[!, nest])
#                 append!(toExclude, rotation.nestings[nest])
#                 if rotation.xAxisVar in toExclude || rotation.yAxisYar in toExclude
#                     push!(toExclude, :axisInteraction)
#                 end

#                 groupColumns = setdiff(columns, toExclude)
#                 for group in groups

#                     #### Group's units
#                     unitsInGroup = filter(filteredUnitModel) do unitRow
#                         return unitRow[nest] == group
#                     end

#                     ### (re-)Predict the column
#                     groupX = DataFrame(:Intercept => [1 for i in 1:nrow(unitsInGroup)])
#                     if !isempty(groupColumns)
#                         groupX = hcat(groupX, unitsInGroup[!, groupColumns])
#                     end

#                     groupX = Matrix{Float64}(groupX)
#                     groupy = Vector{Float64}(map(filter(x->x[2][nest]==group, collect(enumerate(eachrow(filteredUnitModel))))) do (i, unitRow)
#                         return y[i]
#                     end)

#                     println(group)
#                     println(groupColumns)
#                     coefs = (transpose(groupX) * groupX)^-1 * transpose(groupX) * groupy
#                     yhat = groupX * transpose(coefs)

#                     ### Update the column by copying yhat into it at the right locations
#                     j = 1
#                     for (i, unitRow) in enumerate(eachrow(filteredUnitModel))
#                         if unitRow[nest] == group
#                             y[i] = yhat[j]
#                             j += 1
#                         end
#                     end
#                 end
#             end
            
#             ## Run the regression
#             coefs = X * y
#             if T == 1
#                 networkRow[:weight_x] = coefs[2]
#             else
#                 networkRow[:weight_y] = coefs[3]
#             end
#         end

#         ## Normalize the weights
#         if T == 1
#             s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
#             if s != 0
#                 networkModel[!, :weight_x] /= s
#             end
#         else
#             s = sqrt(sum(networkModel[!, :weight_y] .^ 2))
#             if s != 0
#                 networkModel[!, :weight_y] /= s
#             end
#         end

#         ## Find the values for each unit projected onto the x-axis
#         if T == 1
#             xAxisColumn =
#                 Matrix{Float64}(filteredUnitModel[!, networkModel[!, :relationship]]) *
#                 Vector{Float64}(networkModel[!, :weight_x])
#         end
#     end
# end