"""
TODO document
"""
struct Means2Rotation <: AbstractMeans2Rotation
    regression_model::Type{LinearModel}
    coefindex::Int
    f1::FormulaTerm
    contrasts::Union{Nothing,Dict}

    regression_model2::Type{LinearModel}
    coefindex2::Int
    f2::FormulaTerm
    contrasts2::Union{Nothing,Dict}

    groupVar::Symbol
    controlGroup::Any
    treatmentGroup::Any

    groupVar2::Symbol
    controlGroup2::Any
    treatmentGroup2::Any
end

# Simplified constructor
function Means2Rotation(
    groupVar::Symbol, controlGroup::Any, treatmentGroup::Any,
    groupVar2::Symbol, controlGroup2::Any, treatmentGroup2::Any
)
    ## Always use a univariate model for the formula rotation
    regression_model = LinearModel
    coefindex = 2
    f1 = @formula(y ~ 1 + MCFactoredGroupVar + MCFactoredGroupVar2 + InteractionOfMCFactors)
    contrasts = nothing
    regression_model2 = LinearModel
    coefindex2 = 3
    f2 = @formula(y ~ 1 + MCFactoredGroupVar + MCFactoredGroupVar2 + InteractionOfMCFactors)
    contrasts2 = nothing
    return Means2Rotation(
        regression_model, coefindex, f1, contrasts,
        regression_model2, coefindex2, f2, contrasts2,
        groupVar, controlGroup, treatmentGroup,
        groupVar2, controlGroup2, treatmentGroup2
    )
end

# Implement rotate
function rotate!(rotation::AbstractMeans2Rotation, networkModel::DataFrame, unitModel::DataFrame, metadata::DataFrame)

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

    metadata[!, :FactoredGroupVar2] = map(eachrow(metadata)) do unitRow
        if unitRow[rotation.groupVar2] == rotation.controlGroup2
            return 0.0
        elseif unitRow[rotation.groupVar2] == rotation.treatmentGroup2
            return 1.0
        else
            return missing
        end
    end

    ## Mean center the factors
    metadata[!, :MCFactoredGroupVar] = metadata[!, :FactoredGroupVar] .- mean(metadata[!, :FactoredGroupVar])
    metadata[!, :MCFactoredGroupVar2] = metadata[!, :FactoredGroupVar2] .- mean(metadata[!, :FactoredGroupVar2])

    ## Interact the MCs
    metadata[!, :InteractionOfMCFactors] = metadata[!, :MCFactoredGroupVar] .* metadata[!, :MCFactoredGroupVar2]

    ## Use a Formula2Rotation to do the rest of the work
    invoke(rotate!, Tuple{AbstractFormula2Rotation, DataFrame, DataFrame, DataFrame}, rotation, networkModel, unitModel, metadata)
end

# Override plotting pieces
## Base - Inject a groupBy and labels when none are given
function plot(ena::AbstractENAModel{<:AbstractMeans2Rotation};
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
        ylabel = ena.rotation.groupVar2
    end

    return invoke(plot, Tuple{AbstractENAModel{<:AbstractFormula2Rotation}}, ena;
                  negColor=negColor, posColor=posColor, extraColors=extraColors,
                  groupBy=groupBy, xlabel=xlabel, ylabel=ylabel, kwargs...)
end

## Units - different default for labels (identical to the MeansRotation implementation)
function plot_units!(p::Plot, ena::AbstractENAModel{<:AbstractMeans2Rotation}, displayRows::Array{Bool,1};
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

### CIs - also draw triangles for the sub groups
function plot_cis!(p::Plot, ena::AbstractENAModel{<:AbstractMeans2Rotation}, displayRows::Array{Bool,1}, groupName::String;
    color::Colorant=colorant"black",
    flipX::Bool=false, flipY::Bool=false,
    kwargs...)

    #### Whole group
    # xs = ena.centroidModel[displayRows, :pos_x] * (flipX ? -1 : 1)
    # ys = ena.centroidModel[displayRows, :pos_y] * (flipY ? -1 : 1)
    # help_plot_ci(p, xs, ys, color, :square, "$(groupName) Mean")

    #### Partition by the second group var
    controlRows2 = map(ena.metadata[!, ena.rotation.groupVar2]) do group
        if group == ena.rotation.controlGroup2
            return true
        else
            return false
        end
    end

    treatmentRows2 = map(ena.metadata[!, ena.rotation.groupVar2]) do group
        if group == ena.rotation.treatmentGroup2
            return true
        else
            return false
        end
    end

    controlDisplayRows2 = displayRows .& controlRows2
    treatmentDisplayRows2 = displayRows .& treatmentRows2

    #### Show them as up/down triangles
    xs = ena.centroidModel[controlDisplayRows2, :pos_x] * (flipX ? -1 : 1)
    ys = ena.centroidModel[controlDisplayRows2, :pos_y] * (flipY ? -1 : 1)
    help_plot_ci(p, xs, ys, color, :dtriangle, "$(groupName) Mean where $(ena.rotation.groupVar2) = $(ena.rotation.controlGroup2)")

    xs = ena.centroidModel[treatmentDisplayRows2, :pos_x] * (flipX ? -1 : 1)
    ys = ena.centroidModel[treatmentDisplayRows2, :pos_y] * (flipY ? -1 : 1)
    help_plot_ci(p, xs, ys, color, :utriangle, "$(groupName) Mean where $(ena.rotation.groupVar2) = $(ena.rotation.treatmentGroup2)")
end