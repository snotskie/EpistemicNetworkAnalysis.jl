"""
TODO document
"""
struct MeansRotation <: AbstractMeansRotation
    regression_model::Type{LinearModel}
    coefindex::Int
    f1::FormulaTerm
    contrasts::Union{Nothing,Dict}
    groupVar::Symbol
    controlGroup::Any
    treatmentGroup::Any
end

# Simplified constructor
function MeansRotation(groupVar::Symbol, controlGroup::Any, treatmentGroup::Any)
    ## Always use a univariate model for the formula rotation
    regression_model = LinearModel
    coefindex = 2
    f1 = @formula(y ~ 1 + FactoredGroupVar)
    contrasts = nothing
    return MeansRotation(regression_model, coefindex, f1, contrasts, groupVar, controlGroup, treatmentGroup)
end

# Implement rotate
function rotate!(rotation::AbstractMeansRotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)

    ## Manually factor the grouping variable to 0/1/missing
    metadata[!, :FactoredGroupVar] = map(eachrow(metadata)) do unitRow
        if unitRow[rotation.groupVar] == rotation.controlGroup
            return 0.0
        elseif unitRow[rotation.groupVar] == rotation.treatmentGroup
            return 1.0
        else
            return missing
        end
    end

    ## Use a FormulaRotation to do the rest of the work
    invoke(rotate!, Tuple{AbstractFormulaRotation, DataFrame, DataFrame, DataFrame}, rotation, networkModel, unitModel, metadata)
end

# Override plotting pieces
## Base - Inject a groupBy and some labels when none are given
function plot(ena::AbstractENAModel{<:AbstractMeansRotation};
    negColor::Colorant=DEFAULT_NEG_COLOR, posColor::Colorant=DEFAULT_POS_COLOR,
    extraColors::Array{<:Colorant,1}=DEFAULT_EXTRA_COLORS,
    groupBy=nothing,
    xlabel=nothing, ylabel=nothing,
    kwargs...)

    if isnothing(groupBy) || groupBy == ena.rotation.groupVar
        groupBy = ena.rotation.groupVar
        extraColors = [negColor, posColor, extraColors...]
    end

    if isnothing(xlabel)
        xlabel = ena.rotation.groupVar
    end

    if isnothing(ylabel)
        ylabel = "SVD"
    end

    return invoke(plot, Tuple{AbstractENAModel{<:AbstractFormulaRotation}}, ena;
                  negColor=negColor, posColor=posColor, extraColors=extraColors,
                  groupBy=groupBy, xlabel=xlabel, ylabel=ylabel, kwargs...)
end

## Units - different default for labels
function plot_units!(p::Plot, ena::AbstractENAModel{<:AbstractMeansRotation}, displayRows::Array{Bool,1};
    flipX::Bool=false, flipY::Bool=false, minLabel::Union{Nothing,String}=nothing, maxLabel::Union{Nothing,String}=nothing,
    kwargs...)

    ### Use meaningful legend labels for the units
    if isnothing(minLabel)
        minLabel = "$(ena.rotation.controlGroup) Units"
    end

    if isnothing(maxLabel)
        maxLabel = "$(ena.rotation.treatmentGroup) Units"
    end

    ### Let formula rotation do the rest of the work
    invoke(plot_units!, Tuple{Plot, AbstractENAModel{<:AbstractFormulaRotation}, Array{Bool,1}},
        p, ena, displayRows; flipX=flipX, flipY=flipY, minLabel=minLabel, maxLabel=maxLabel, kwargs...)
end