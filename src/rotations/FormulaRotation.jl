abstract type AbstractFormulaRotation <: AbstractLinearENARotation end
struct FormulaRotation <: AbstractFormulaRotation
    regression_models::Array{Type{<:RegressionModel}}
    formulas::Array{<:FormulaTerm}
    coef_indexes::Array{<:Int}
    contrasts::Array{<:Union{Nothing,Dict}}
end

"""
    FormulaRotation(
        regression_model1::Type{T},
        formula1::FormulaTerm,
        coef_index1::Int,
        contrast1::Union{Nothing,Dict},
        args...
    ) where {T <: RegressionModel}

Define a rotation that uses regression models to determine axes most closely associated with some linear trend

Note: `RegressionModel`s must be imported from other stats packages

Note: [contrasts](https://juliastats.org/StatsModels.jl/stable/contrasts/) are used to model categorical data

## Example

```julia
using GLM

rotation = EpistemicNetworkAnalysis.FormulaRotation(
    LinearModel, @formula(edge ~ 1 + FinalGrade), 2, nothing
)
```

This will fit the x-axis to the `FinalGrade` metadata, because:

1. We use `LinearModel` from the `GLM` package
2. We use the formula `edge ~ 1 + FinalGrade`
3. And we use the `2`nd coefficient of the `LinearModel` (in this case `FinalGrade`) to determine the values of the embedding
4. We have no categorical data in the `LinearModel`, so we leave the contrasts as `nothing`

Additional formulae may be used to define subsequent axes:

```julia
rotation = EpistemicNetworkAnalysis.FormulaRotation(
    LinearModel, @formula(edge ~ 1 + PretestGrade + PosttestGrade), 2, nothing,
    LinearModel, @formula(edge ~ 1 + PretestGrade + PosttestGrade), 3, nothing
)
```

Note: When multiple formulae are given, `FormulaRotation` finds the *plane* of the effects in the `accum` space, rotating it such that the first formula aligns with the x-axis, and the second formula aligns *approximately* with the y-axis. Unless this approximation is strong, a warning will be raised describing possible issues.

## Statistical Tests

Models using a `MeansRotation` will run the following statistical tests:

- R^2 and adjusted-R^2 for each dimension with a formula
"""
FormulaRotation

function FormulaRotation(
        regression_model1::Type{T},
        formula1::FormulaTerm,
        coef_index1::Int,
        contrast1::Union{Nothing,Dict},
        args...
    ) where {T <: RegressionModel}

    @assert length(args) % 4 == 0 "FormulaRotation expects a multiple of 4 arguments"

    regression_models = Type[regression_model1]
    coef_indexes = Int[coef_index1]
    formulas = FormulaTerm[formula1]
    contrasts = Union{Nothing,Dict}[contrast1]

    for i in 1:4:length(args)
        push!(regression_models, args[i+0])
        push!(coef_indexes, args[i+1])
        if args[i+2] isa Term # BUGFIX
            push!(formulas, [args[i+2]])
        else
            push!(formulas, args[i+2])
        end
        push!(contrasts, args[i+3])
        i += 4
    end

    return FormulaRotation(regression_models, formulas, coef_indexes, contrasts)
end

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractFormulaRotation, M<:AbstractLinearENAModel{R}}

    regressionData = leftjoin(model.accum, model.metadata, on=:unitID)

    ## Filter out rows with missing data - not needed?
    # for formula in model.rotation.formulas
    #     for t in formula.rhs
    #         if t isa Term
    #             col = Symbol(t)
    #             goodRows = completecases(regressionData[!, [col]])
    #             regressionData = regressionData[goodRows, :]
    #         end
    #     end
    # end

    ## BUGFIX: https://github.com/JuliaStats/GLM.jl/issues/239
    edgeIDs = model.edges.edgeID
    for edge in edgeIDs
        regressionData[!, edge] = map(Float64, regressionData[!, edge])
    end

    ## For each formula, for each relationship, find the effect of the chosen predictor, use that as the axis weights
    embedding = similar(model.embedding, length(model.rotation.formulas))
    for (i, formula) in enumerate(model.rotation.formulas)
        embedding[i, :label] = string(formula.rhs[model.rotation.coef_indexes[i]])
        for edge in edgeIDs
            f1 = FormulaTerm(term(edge), formula.rhs)
            try
                if model.rotation.contrasts[i] isa Nothing
                    m1 = fit(model.rotation.regression_models[i], f1, regressionData)
                    slope = coef(m1)[model.rotation.coef_indexes[i]]
                    embedding[i, edge] = slope
                else
                    m1 = fit(model.rotation.regression_models[i], f1, regressionData, contrasts=model.rotation.contrasts[i])
                    slope = coef(m1)[model.rotation.coef_indexes[i]]
                    embedding[i, edge] = slope
                end
            catch e
                println(e)
                error("""
                An error occured running a regression during the model.rotation step of this ENA model.
                Usually, this occurs because the data, the regression model, and regression formula are not in agreement.
                Double check that:
                (1) TODO
                """)
            end
        end
    end

    append!(model.embedding, embedding)

    # let parent handle the rest
    super = rotationsupertype(M, AbstractFormulaRotation)
    rotate!(super, model)
end

function test!(
        ::Type{M}, trainmodel::AbstractLinearENAModel, testmodel::AbstractLinearENAModel
    ) where {R<:AbstractFormulaRotation, M<:AbstractLinearENAModel{R}}

    super = rotationsupertype(M, AbstractFormulaRotation)
    test!(super, trainmodel, testmodel)
    for (i, formula) in enumerate(trainmodel.rotation.formulas)
        test!(M, trainmodel, testmodel, trainmodel.rotation.regression_models[i], dim=i,
            formula=formula,
            contrasts=trainmodel.rotation.contrasts[i]
        )
    end
end