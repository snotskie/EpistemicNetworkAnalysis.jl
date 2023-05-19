abstract type AbstractFormulaRotation <: AbstractLinearENARotation end
struct FormulaRotation <: AbstractFormulaRotation
    regression_models::Array{Type{<:RegressionModel}}
    formulas::Array{<:FormulaTerm}
    coef_indexes::Array{<:Int}
    contrasts::Array{<:Union{Nothing,Dict}}
end

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