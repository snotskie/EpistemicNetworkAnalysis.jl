abstract type AbstractMeansRotation <: AbstractModeratedRotation end
struct MeansRotation <: AbstractMeansRotation
    regression_models::Array{Type{<:RegressionModel}}
    formulas::Array{<:FormulaTerm}
    coef_indexes::Array{<:Int}
    contrasts::Array{<:Union{Nothing,Dict}}
    groupVars::Array{Symbol}
    controlGroups::Array{<:Any}
    treatmentGroups::Array{<:Any}
end

function MeansRotation(
        groupVar1::Symbol,
        controlGroup1::Any,
        treatmentGroup1::Any,
        args...
    )

    @assert length(args) % 3 == 0 "MeansRotation expects a multiple of 3 arguments"

    groupVars = Symbol[groupVar1]
    controlGroups = Any[controlGroup1]
    treatmentGroups = Any[treatmentGroup1]
    f1 = @formula(y ~ 1)
    f1 = FormulaTerm(f1.lhs, f1.rhs + Term(Symbol(string("MCFactored_", groupVar1))))

    i = 0
    while i < length(args)
        push!(groupVars, args[i+0])
        push!(controlGroups, args[i+1])
        push!(treatmentGroups, args[i+2])
        f1 = FormulaTerm(f1.lhs, f1.rhs + Term(Symbol(string("MCFactored_", args[i+0]))))
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
        treatmentGroups
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

        # col_data = Vector{Float64}(col_data) .- mean(col_data)
        model.metadata[!, col_name] = col_data
    end

    # let parent handle the rest
    super = rotationsupertype(M, AbstractMeansRotation)
    rotate!(super, model)
end

# # Simplified constructor
# function MeansRotation(groupVar::Symbol, controlGroup::Any, treatmentGroup::Any)
#     ## Always use a univariate model for the formula rotation
#     regression_model = LinearModel
#     coefindex = 2
#     f1 = @formula(y ~ 1 + FactoredGroupVar)
#     contrasts = nothing
#     return MeansRotation(regression_model, coefindex, f1, contrasts, groupVar, controlGroup, treatmentGroup)
# end

# # Implement rotate
# function rotate!(rotation::AbstractMeansRotation, networkModel::DataFrame, codeModel::DataFrame, metadata::DataFrame, subspaceModel::DataFrame)

#     ## Manually factor the grouping variable to 0/1/missing
#     metadata[!, :FactoredGroupVar] = map(eachrow(metadata)) do unitRow
#         if unitRow[rotation.groupVar] == rotation.controlGroup
#             return 0.0
#         elseif unitRow[rotation.groupVar] == rotation.treatmentGroup
#             return 1.0
#         else
#             return missing
#         end
#     end

#     ## Use a FormulaRotation to do the rest of the work
#     invoke(rotate!, Tuple{AbstractFormulaRotation, DataFrame, DataFrame, DataFrame, DataFrame}, rotation, networkModel, codeModel, metadata, subspaceModel)
# end

# # Override plotting pieces
# ## Base - Inject a groupBy and some labels when none are given
# function plot(ena::AbstractENAModel{<:AbstractMeansRotation};
#     negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
#     extraColors::Array{<:Colorant,1}=DEFAULT_EXTRA_COLORS,
#     groupBy=nothing,
#     xlabel=nothing, ylabel=nothing,
#     kwargs...)

#     if isnothing(groupBy) || groupBy == ena.rotation.groupVar
#         groupBy = ena.rotation.groupVar
#         if ena.rotation.controlGroup < ena.rotation.treatmentGroup
#             extraColors = [negColor, posColor, extraColors...]
#         else
#             extraColors = [posColor, negColor, extraColors...]
#         end
#     end

#     if isnothing(xlabel)
#         xlabel = ena.rotation.groupVar
#     end

#     if isnothing(ylabel)
#         ylabel = "SVD"
#     end

#     return invoke(plot, Tuple{AbstractENAModel{<:AbstractFormulaRotation}}, ena;
#                   negColor=negColor, posColor=posColor, extraColors=extraColors,
#                   groupBy=groupBy, xlabel=xlabel, ylabel=ylabel, kwargs...)
# end

# ## Units - different default for labels
# function plot_units!(p::Plot, ena::AbstractENAModel{<:AbstractMeansRotation}, displayRows::Array{Bool,1};
#     flipX::Bool=false, flipY::Bool=false, minLabel::Union{Nothing,String}=nothing, maxLabel::Union{Nothing,String}=nothing,
#     kwargs...)

#     ### Use meaningful legend labels for the units
#     if isnothing(minLabel)
#         minLabel = "$(ena.rotation.controlGroup) Units"
#     end

#     if isnothing(maxLabel)
#         maxLabel = "$(ena.rotation.treatmentGroup) Units"
#     end

#     ### Let formula rotation do the rest of the work
#     invoke(plot_units!, Tuple{Plot, AbstractENAModel{<:AbstractFormulaRotation}, Array{Bool,1}},
#         p, ena, displayRows; flipX=flipX, flipY=flipY, minLabel=minLabel, maxLabel=maxLabel, kwargs...)
# end

# function test(ena::AbstractENAModel{<:AbstractMeansRotation})
#     results = invoke(test, Tuple{AbstractENAModel{<:AbstractFormulaRotation}}, ena)
#     controlRows = ena.metadata[!, ena.rotation.groupVar] .== ena.rotation.controlGroup
#     treatmentRows = ena.metadata[!, ena.rotation.groupVar] .== ena.rotation.treatmentGroup
#     controlVals = ena.accumModel[controlRows, :pos_x]
#     treatmentVals = ena.accumModel[treatmentRows, :pos_x]
#     mwut = MannWhitneyUTest(controlVals, treatmentVals)
#     results[:mann_whitney_u] = mwut.U
#     results[:mann_whitney_pvalue] = pvalue(mwut)
#     results[:n_control] = length(controlVals)
#     results[:n_treatment] = length(treatmentVals)
#     results[:mdn_control] = median(controlVals)
#     results[:mdn_treatment] = median(treatmentVals)
#     return results
# end
