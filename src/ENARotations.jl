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
struct Formula2Rotation{T <: RegressionModel, U <: RegressionModel} <: ENARotation
    regression_model1::Type{T}
    f1::FormulaTerm
    regression_model2::Type{U}
    f2::FormulaTerm
end

# SVD Rotation
function rotate!(::SVDRotation, networkModel::DataFrame, unitModel::DataFrame)
    pcaModel = projection(help_deflating_svd(networkModel, unitModel))
    networkModel[!, :weight_x] = pcaModel[:, 1]
    networkModel[!, :weight_y] = pcaModel[:, 2]
end

# Means Rotation
function rotate!(rotation::MeansRotation, networkModel::DataFrame, unitModel::DataFrame)
    ## Filter the unit model to just those in the control/treatment group
    filteredUnitModel = filter(unitModel) do unitRow
        return unitRow[rotation.groupVar] in [rotation.controlGroup, rotation.treatmentGroup]
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

    ## Filter missing data
    filteredUnitModel = unitModel
    for t in rotation.f1.rhs
        if isa(t, Term)
            col = Symbol(t)
            goodRows = completecases(filteredUnitModel[!, [col]])
            filteredUnitModel = filteredUnitModel[goodRows, :]
        end
    end

    ## Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        filteredUnitModel[!, r] = map(Float64, filteredUnitModel[!, r])
    end

    ## For each relationship, find the effect of the first predictor after the intercept
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        f1 = FormulaTerm(term(r), rotation.f1.rhs)
        try
            m1 = fit(rotation.regression_model, f1, filteredUnitModel)
            slope = coef(m1)[2]
            networkRow[:weight_x] = slope
        catch e
            error("""
            An error occured running a regression during the rotation step of this ENA model.
            Usually, this occurs because the data, the regression model, and regression formula are not in agreement.
            If you are using a MeansRotation, then this usually means that your accidentally grouped your
            units on a different variable than the variable you passed to your MeansRotation.
            """)
        end
    end

    ## Normalize the weights
    s = sqrt(sum(networkModel[!, :weight_x] .^ 2))
    if s != 0
        networkModel[!, :weight_x] /= s
    end

    ## Find the first svd dim of the data orthogonal to the x weights, use these as the y weights
    xAxis = Matrix{Float64}(unitModel[!, networkModel[!, :relationship]]) *
            Matrix{Float64}(networkModel[!, [:weight_x]])
    xAxis = xAxis .- mean(xAxis)
    controlModel = DataFrame(xAxis)
    pcaModel = projection(help_deflating_svd(networkModel, unitModel, controlModel))
    networkModel[!, :weight_y] = pcaModel[:, 1]
end

# Formula2 Rotation
function rotate!(rotation::Formula2Rotation, networkModel::DataFrame, unitModel::DataFrame)
    ## TODO check assumptions about f1 and f2

    ## Filter missing data
    filteredUnitModel = unitModel
    for t in rotation.f1.rhs
        if isa(t, Term)
            col = Symbol(t)
            goodRows = completecases(filteredUnitModel[!, [col]])
            filteredUnitModel = filteredUnitModel[goodRows, :]
        end
    end

    for t in rotation.f2.rhs
        if isa(t, Term)
            col = Symbol(t)
            goodRows = completecases(filteredUnitModel[!, [col]])
            filteredUnitModel = filteredUnitModel[goodRows, :]
        end
    end

    ## Bugfix: https://github.com/JuliaStats/GLM.jl/issues/239
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        filteredUnitModel[!, r] = map(Float64, filteredUnitModel[!, r])
    end

    ## For each relationship, find the effect of the first predictor after the intercept
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        f1 = FormulaTerm(term(r), rotation.f1.rhs)
        try
            m1 = fit(rotation.regression_model1, f1, filteredUnitModel)
            slope = coef(m1)[2]
            networkRow[:weight_x] = slope
        catch e
            error("""
            An error occured running a regression during the rotation step of this ENA model.
            Usually, this occurs because the data, the regression model, and regression formula are not in agreement.
            """)
        end
    end

    ## Again, but for the y-axis and using the second formula
    for networkRow in eachrow(networkModel)
        r = networkRow[:relationship]
        f2 = FormulaTerm(term(r), rotation.f2.rhs)
        try
            m2 = fit(rotation.regression_model2, f2, filteredUnitModel)
            slope = coef(m2)[2]
            networkRow[:weight_y] = slope
        catch e
            error("""
            An error occured running a regression during the rotation step of this ENA model.
            Usually, this occurs because the data, the regression model, and regression formula are not in agreement.
            """)
        end
    end

    ## Orthogonalization: replace y weights with their rejection from the x weights
    before = copy(networkModel[!, :weight_y])
    scalar = dot(networkModel[!, :weight_y], networkModel[!, :weight_x]) / dot(networkModel[!, :weight_x], networkModel[!, :weight_x])
    networkModel[!, :weight_y] -= scalar * networkModel[!, :weight_x]    
    after = copy(networkModel[!, :weight_y])

    ## Raise a warning about interpreting the y-axis when before and after have a large angle between them
    theta = dot(before, after)
    theta /= sqrt(dot(before, before))
    theta /= sqrt(dot(after, after))
    angle = acos(theta) * 180 / pi
    if abs(angle) > 5
        @warn """The angle between the y-axis and the direction of the requested effect is larger than 5 degrees ($angle degrees).
This can undermine interpreting the y-axis in terms of the requested effect."""
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