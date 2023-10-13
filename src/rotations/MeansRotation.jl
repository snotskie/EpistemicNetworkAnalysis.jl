abstract type AbstractMeansRotation <: AbstractFormulaRotation end
struct MeansRotation <: AbstractMeansRotation
    regression_models::Array{Type{<:RegressionModel}}
    formulas::Array{<:FormulaTerm}
    coef_indexes::Array{<:Int}
    contrasts::Array{<:Union{Nothing,Dict}}
    groupVars::Array{Symbol}
    controlGroups::Array{<:Any}
    treatmentGroups::Array{<:Any}
    moderated::Bool
end

"""
    MeansRotation(
        # Required, groups to compare along the x-axis
        groupVar1::Symbol,
        controlGroup1::Any,
        treatmentGroup1::Any,

        # Optional, groups to compare along subsequent axes
        args...;

        # Optional, whether to moderate the interactions between group dimensions
        moderated=false
    )

Define a rotation for comparing pairs of groups, by maximizing the variance between pairs

See also: `LDARotation` and `MulticlassRotation`

## Example

```julia
rotation = MeansRotation(
    :Play, "Romeo and Juliet", "Hamlet",
    :Act, 1, 5
)
```

## Statistical Tests

Models using a `MeansRotation` will run the following statistical tests:

- `KruskalWallisTest` for each dimension with a group pair
"""
MeansRotation

function MeansRotation(
        groupVar1::Symbol,
        controlGroup1::Any,
        treatmentGroup1::Any,
        args...;
        moderated=false
    )

    @assert length(args) % 3 == 0 "MeansRotation expects a multiple of 3 arguments"

    groupVars = Symbol[groupVar1]
    controlGroups = Any[controlGroup1]
    treatmentGroups = Any[treatmentGroup1]
    f1 = @formula(y ~ 1)
    f1 = FormulaTerm(f1.lhs, f1.rhs + Term(Symbol(string("MCFactored_", groupVar1))))

    for i in 1:4:length(args)
        push!(groupVars, args[i+0])
        push!(controlGroups, args[i+1])
        push!(treatmentGroups, args[i+2])
        fterm = Term(Symbol(string("MCFactored_", args[i+0])))
        f1 = FormulaTerm(f1.lhs, f1.rhs + fterm)
        i += 4
    end

    if moderated
        for i in 1:length(controlGroups)
            for j in (i+1):length(controlGroups)
                iterm = Term(Symbol(string("MCInteraction_", groupVars[i], "_", groupVars[j])))
                f1 = FormulaTerm(f1.lhs, f1.rhs + iterm)
            end
        end
    end

    N = length(groupVars)
    regression_models = repeat([LinearModel], N)
    formulas = repeat([f1], N)
    coef_indexes = [1:N...] .+ 1
    contrasts = repeat([nothing], N)

    return MeansRotation(
        regression_models,
        formulas,
        coef_indexes,
        contrasts,
        groupVars,
        controlGroups,
        treatmentGroups,
        moderated
    )
end

function rotate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractMeansRotation, M<:AbstractLinearENAModel{R}}

    ## Manually factor the grouping variables to 0/1/missing, mean centered
    for (i, coef_index) in enumerate(model.rotation.coef_indexes)
        col_name = Symbol(model.rotation.formulas[i].rhs[coef_index])
        groupVar = model.rotation.groupVars[i]
        controlGroup = model.rotation.controlGroups[i]
        treatmentGroup = model.rotation.treatmentGroups[i]
        col_data = map(eachrow(model.metadata)) do unitRow
            if unitRow[groupVar] == controlGroup
                return 0.0
            elseif unitRow[groupVar] == treatmentGroup
                return 1.0
            else
                return missing
            end
        end

        col_data = Vector(col_data) .- mean(skipmissing(col_data))
        model.metadata[!, col_name] = col_data
    end

    ## optionally compute interactions between those terms
    if model.rotation.moderated
        coef_index_interaction = maximum(model.rotation.coef_indexes) + 1
        for (i, coef_index_i) in enumerate(model.rotation.coef_indexes)
            col_name_i = Symbol(model.rotation.formulas[i].rhs[coef_index_i])
            for (j, coef_index_j) in enumerate(model.rotation.coef_indexes[i+1:end])
                col_name_j = Symbol(model.rotation.formulas[j].rhs[coef_index_j])
                col_name_interaction = Symbol(model.rotation.formulas[j].rhs[coef_index_interaction])
                col_data = Vector(model.metadata[!, col_name_i]) .* Vector(model.metadata[!, col_name_j])
                model.metadata[!, col_name_interaction] = col_data
                coef_index_interaction += 1
            end
        end
    end

    # let parent handle the rest
    super = rotationsupertype(M, AbstractMeansRotation)
    rotate!(super, model)

    # update embedding labels
    for (i, label) in enumerate(model.rotation.groupVars)
        model.embedding[i, :label] = string(label)
    end
end

function test!(
        ::Type{M}, trainmodel::AbstractLinearENAModel, testmodel::AbstractLinearENAModel
    ) where {R<:AbstractMeansRotation, M<:AbstractLinearENAModel{R}}

    # BUGFIX https://github.com/snotskie/EpistemicNetworkAnalysis.jl/issues/6
    # when we have a true train/test split, we have to "redo" the work of
    # factoring the model so the testmodel.metadata has columns expected by
    # FormulaRotation's tests
    if trainmodel != testmodel
        ## factor the grouping variables to 0/1/missing, mean centered
        for (i, coef_index) in enumerate(trainmodel.rotation.coef_indexes)
            col_name = Symbol(trainmodel.rotation.formulas[i].rhs[coef_index])
            groupVar = trainmodel.rotation.groupVars[i]
            controlGroup = trainmodel.rotation.controlGroups[i]
            treatmentGroup = trainmodel.rotation.treatmentGroups[i]
            col_data = map(eachrow(testmodel.metadata)) do unitRow
                if unitRow[groupVar] == controlGroup
                    return 0.0
                elseif unitRow[groupVar] == treatmentGroup
                    return 1.0
                else
                    return missing
                end
            end

            col_data = Vector(col_data) .- mean(skipmissing(col_data))
            testmodel.metadata[!, col_name] = col_data
        end

        ## compute interactions between terms
        if trainmodel.rotation.moderated
            coef_index_interaction = maximum(trainmodel.rotation.coef_indexes) + 1
            for (i, coef_index_i) in enumerate(trainmodel.rotation.coef_indexes)
                col_name_i = Symbol(trainmodel.rotation.formulas[i].rhs[coef_index_i])
                for (j, coef_index_j) in enumerate(trainmodel.rotation.coef_indexes[i+1:end])
                    col_name_j = Symbol(trainmodel.rotation.formulas[j].rhs[coef_index_j])
                    col_name_interaction = Symbol(trainmodel.rotation.formulas[j].rhs[coef_index_interaction])
                    col_data = Vector(testmodel.metadata[!, col_name_i]) .* Vector(testmodel.metadata[!, col_name_j])
                    testmodel.metadata[!, col_name_interaction] = col_data
                    coef_index_interaction += 1
                end
            end
        end
    end

    super = rotationsupertype(M, AbstractMeansRotation)
    test!(super, trainmodel, testmodel)
    for (i, label) in enumerate(trainmodel.rotation.groupVars)
        test!(M, trainmodel, testmodel, GroupwiseCoregistrationTest, dim=i, groupVar=label, groups=[
            trainmodel.rotation.controlGroups[i],
            trainmodel.rotation.treatmentGroups[i]
        ])

        test!(M, trainmodel, testmodel, KruskalWallisTest, dim=i, groupVar=label, groups=[
            trainmodel.rotation.controlGroups[i],
            trainmodel.rotation.treatmentGroups[i]
        ])
    end
end

# insert colors and groups
function defaultplotkwargs(
        ::Type{M},
        model::AbstractLinearENAModel;
        x::Int=1,
        y::Int=2,
        # negColor::Colorant=DEFAULT_NEG_COLOR,
        # posColor::Colorant=DEFAULT_POS_COLOR,
        # extraColors::Array{<:Colorant,1}=(
        #     # ensure the "left" group is the "red" group, using alphabetical order
        #     model.rotation.controlGroups[x] <= model.rotation.treatmentGroups[x] ?
        #     [DEFAULT_NEG_COLOR, DEFAULT_POS_COLOR, DEFAULT_EXTRA_COLORS...] :
        #     [DEFAULT_POS_COLOR, DEFAULT_NEG_COLOR, DEFAULT_EXTRA_COLORS...]
        # ),
        # if x has a group var, group by that
        # elif y has a group var, group by that
        # else, nothing
        groupBy::Union{Symbol,Nothing}=(
            x <= length(model.rotation.groupVars) ?
            model.rotation.groupVars[x] :
            (
                y <= length(model.rotation.groupVars) ?
                model.rotation.groupVars[y] :
                nothing
            )
        ),
        # if x and y have group vars, inner group by y
        innerGroupBy::Union{Symbol,Nothing}=(
            x <= length(model.rotation.groupVars) && y <= length(model.rotation.groupVars) ?
            model.rotation.groupVars[y] :
            nothing
        ),
        kwargs...
    ) where {R<:AbstractMeansRotation, M<:AbstractLinearENAModel{R}}

    kwargs = NamedTuple(kwargs)
    defaults = (
        x=x,
        y=y,
        # negColor=negColor,
        # posColor=posColor,
        # extraColors=extraColors,
        groupBy=groupBy,
        innerGroupBy=innerGroupBy,
        kwargs...
    )

    super = rotationsupertype(M, AbstractMeansRotation)
    return defaultplotkwargs(super, model, merge(defaults, kwargs))
end
