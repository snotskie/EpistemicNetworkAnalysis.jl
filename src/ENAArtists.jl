# TODO show one ENAUnit's model
# TODO refactor: reduce multiple artists to calling a single, generic LinesAboveAverageArtist or similar
# TODO refactor: don't make callable objects, overload a paint! or similar function

"""
TODO: document
"""
abstract type ENAArtist
end

"""
TODO: document
"""
struct DefaultArtist <: ENAArtist
end

"""
TODO: document
"""
struct MeansArtist <: ENAArtist
    groupVar::Symbol
    controlGroup::Any
    treatmentGroup::Any
end

"""
TODO: document
"""
struct WindowsArtist <: ENAArtist
    groupVar1::Symbol
    controlGroup1::Any
    treatmentGroup1::Any
    groupVar2::Symbol
    controlGroup2::Any
    treatmentGroup2::Any
end

"""
TODO: document
"""
struct TVRemoteArtist <: ENAArtist
    groupVar1::Symbol
    controlGroup1::Any
    treatmentGroup1::Any
    groupVar2::Symbol
    controlGroup2::Any
    treatmentGroup2::Any
end

# DefaultArtist
function (artist::DefaultArtist)(cb, ena, scene)
    ## Colors
    networkColors = map(eachrow(ena.networkModel)) do networkRow
        return :black
    end

    unitColors = map(eachrow(ena.refitUnitModel)) do refitRow
        return :black
    end

    codeColors = map(eachrow(ena.codeModel)) do codeRow
        return :black
    end

    ## Shapes
    unitShapes = map(eachrow(ena.refitUnitModel)) do refitRow
        return :circle
    end

    codeShapes = map(eachrow(ena.codeModel)) do codeRow
        return :circle
    end

    ## Sizes
    networkLineWidths = map(eachrow(ena.networkModel)) do networkRow
        return networkRow[:density]
    end

    unitMarkerSizes = map(eachrow(ena.refitUnitModel)) do refitRow
        return 1
    end

    codeMarkerSizes = map(eachrow(ena.codeModel)) do codeRow
        return codeRow[:density]
    end

    ## Normalize
    s = maximum(networkLineWidths)
    networkLineWidths /= s
    networkLineWidths *= 2

    s = maximum(unitMarkerSizes)
    unitMarkerSizes /= s
    unitMarkerSizes *= 1.5

    s = maximum(codeMarkerSizes)
    codeMarkerSizes /= s
    codeMarkerSizes *= 8

    ## Confidence Intervals
    confidenceIntervals = []
    color = :black
    shape = :square
    size = 4
    if nrow(ena.refitUnitModel) > 1
        mu_x = mean(ena.refitUnitModel[!, :pos_x])
        mu_y = mean(ena.refitUnitModel[!, :pos_y])
        ci_x = collect(confint(OneSampleTTest(ena.refitUnitModel[!, :pos_x])))
        ci_y = collect(confint(OneSampleTTest(ena.refitUnitModel[!, :pos_y])))
        CI = (mu_x, mu_y, ci_x, ci_y, color, shape, size)
        push!(confidenceIntervals, CI)
    elseif nrow(ena.refitUnitModel) == 1
        mu_x = mean(ena.refitUnitModel[!, :pos_x])
        mu_y = mean(ena.refitUnitModel[!, :pos_y])
        CI = (mu_x, mu_y, [mu_x, mu_x], [mu_y, mu_y], color, shape, size)
        push!(confidenceIntervals, CI)
    else
        # do nothing
    end

    ## Do Callback so ENADisplay can do the bulk of the work
    cb(networkColors,
       unitColors,
       codeColors,
       unitShapes,
       codeShapes,
       networkLineWidths,
       unitMarkerSizes,
       codeMarkerSizes,
       confidenceIntervals)
end

# MeansArtist
function (artist::MeansArtist)(cb, ena, scene)
    ##Pre-splitting the groups
    controlUnits = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar] == artist.controlGroup
            return true
        else
            return false
        end
    end

    treatmentUnits = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar] == artist.treatmentGroup
            return true
        else
            return false
        end
    end

    # ## Pre-processing the size of the network lines
    # ### Find the "direction", "strength", and "angle" for the line size
    # ### Thicker lines are those whose rotation weights are more towards
    # ### one side of the difference of the means
    # mu_x_control = mean(controlUnits[!, :pos_x])
    # mu_y_control = mean(controlUnits[!, :pos_y])
    # mu_x_treatment = mean(treatmentUnits[!, :pos_x])
    # mu_y_treatment = mean(treatmentUnits[!, :pos_y])
    # mu_x_all = (mu_x_treatment + mu_x_control) / 2
    # mu_y_all = (mu_y_treatment + mu_y_control) / 2
    # vt = Vector{Float64}([
    #     mu_x_treatment - mu_x_all,
    #     mu_y_treatment - mu_y_all
    # ])

    # norm_vt = sqrt(dot(vt, vt))
    # lineStrengths = Dict{Symbol,Float64}()
    # for networkRow in eachrow(ena.networkModel)
    #     r = networkRow[:relationship]
    #     vl = Vector{Float64}([
    #         networkRow[:weight_x],
    #         networkRow[:weight_y]
    #     ])

    #     lineStrengths[r] = dot(vl, vt) / norm_vt
    # end

    ## Pre-processing the size of the network lines
    ### Find out which group, on average, has the biggest
    ### "lead" over other groups, on average, for each line.
    #### Use that "lead" as the size for the line
    lineStrengthsControl = Dict{Symbol,Float64}()
    lineStrengthsTreatment = Dict{Symbol,Float64}()
    for networkRow in eachrow(ena.networkModel)
        r = networkRow[:relationship]
        lineStrengthsControl[r] = sum(controlUnits[!, r]) / nrow(controlUnits)
        lineStrengthsTreatment[r] = sum(treatmentUnits[!, r]) / nrow(treatmentUnits)
        if isnan(lineStrengthsControl[r])
            lineStrengthsControl[r] = 0
        end

        if isnan(lineStrengthsTreatment[r])
            lineStrengthsTreatment[r] = 0
        end
    end

    ## Colors
    networkColors = map(eachrow(ena.networkModel)) do networkRow
        r = networkRow[:relationship]
        if lineStrengthsControl[r] - lineStrengthsTreatment[r] > 0
            return :purple
        elseif lineStrengthsControl[r] - lineStrengthsTreatment[r] < 0
            return :orange
        else
            return :white
        end
    end

    unitColors = map(eachrow(ena.refitUnitModel)) do refitRow
        if refitRow[artist.groupVar] == artist.controlGroup
            return :purple
        elseif refitRow[artist.groupVar] == artist.treatmentGroup
            return :orange
        else
            return :black
        end
    end

    codeColors = map(eachrow(ena.codeModel)) do codeRow
        return :black
    end

    ## Shapes
    unitShapes = map(eachrow(ena.refitUnitModel)) do refitRow
        return :circle
    end

    codeShapes = map(eachrow(ena.codeModel)) do codeRow
        return :circle
    end

    ## Sizes
    networkLineWidths = map(eachrow(ena.networkModel)) do networkRow
        # return lineStrengths[networkRow[:relationship]] ^ 2
        r = networkRow[:relationship]
        return abs(lineStrengthsControl[r] - lineStrengthsTreatment[r])
    end

    unitMarkerSizes = map(eachrow(ena.refitUnitModel)) do refitRow
        return 1
    end

    codeMarkerSizes = map(eachrow(ena.codeModel)) do codeRow
        return codeRow[:density]
    end

    ## Normalize
    s = maximum(networkLineWidths)
    networkLineWidths /= s
    networkLineWidths *= 2

    s = maximum(unitMarkerSizes)
    unitMarkerSizes /= s
    unitMarkerSizes *= 1.5

    s = maximum(codeMarkerSizes)
    codeMarkerSizes /= s
    codeMarkerSizes *= 8

    ## Confidence Intervals
    confidenceIntervals = []

    ### Control
    color = :purple
    shape = :square
    size = 4
    if nrow(controlUnits) > 1
        mu_x = mean(controlUnits[!, :pos_x])
        mu_y = mean(controlUnits[!, :pos_y])
        ci_x = collect(confint(OneSampleTTest(controlUnits[!, :pos_x])))
        ci_y = collect(confint(OneSampleTTest(controlUnits[!, :pos_y])))
        CI = (mu_x, mu_y, ci_x, ci_y, color, shape, size)
        push!(confidenceIntervals, CI)
    elseif nrow(controlUnits) == 1
        mu_x = mean(controlUnits[!, :pos_x])
        mu_y = mean(controlUnits[!, :pos_y])
        CI = (mu_x, mu_y, [mu_x, mu_x], [mu_y, mu_y], color, shape, size)
        push!(confidenceIntervals, CI)
    else
        # do nothing
    end

    ### Treatment
    color = :orange
    shape = :square
    size = 4
    if nrow(treatmentUnits) > 1
        mu_x = mean(treatmentUnits[!, :pos_x])
        mu_y = mean(treatmentUnits[!, :pos_y])
        ci_x = collect(confint(OneSampleTTest(treatmentUnits[!, :pos_x])))
        ci_y = collect(confint(OneSampleTTest(treatmentUnits[!, :pos_y])))
        CI = (mu_x, mu_y, ci_x, ci_y, color, shape, size)
        push!(confidenceIntervals, CI)
    elseif nrow(treatmentUnits) == 1
        mu_x = mean(treatmentUnits[!, :pos_x])
        mu_y = mean(treatmentUnits[!, :pos_y])
        CI = (mu_x, mu_y, [mu_x, mu_x], [mu_y, mu_y], color, shape, size)
        push!(confidenceIntervals, CI)
    else
        # do nothing
    end

    ## Do Callback so ENADisplay can do the bulk of the work
    cb(networkColors,
       unitColors,
       codeColors,
       unitShapes,
       codeShapes,
       networkLineWidths,
       unitMarkerSizes,
       codeMarkerSizes,
       confidenceIntervals)
end

# WindowsArtist
function (artist::WindowsArtist)(cb, ena, scene)
    # ##Pre-splitting the groups
    group00Units = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar1] == artist.controlGroup1 &&
           refitRow[artist.groupVar2] == artist.controlGroup2
            return true
        else
            return false
        end
    end

    group01Units = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar1] == artist.controlGroup1 &&
           refitRow[artist.groupVar2] == artist.treatmentGroup2
            return true
        else
            return false
        end
    end

    group10Units = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar1] == artist.treatmentGroup1 &&
           refitRow[artist.groupVar2] == artist.controlGroup2
            return true
        else
            return false
        end
    end

    group11Units = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar1] == artist.treatmentGroup1 &&
           refitRow[artist.groupVar2] == artist.treatmentGroup2
            return true
        else
            return false
        end
    end

    allGroupedUnits = [group00Units,
        group01Units,
        group10Units,
        group11Units]
    
    ## Pre-processing the line weights and colors
    lineStrengths = Dict{Symbol,Float64}(
        r => -1
        for r in ena.networkModel[!, :relationship]
    )

    lineColors = Dict{Symbol,Symbol}(
        r => :black
        for r in ena.networkModel[!, :relationship]
    )

    for networkRow in eachrow(ena.networkModel)
        r = networkRow[:relationship]
        colors = [:blue, :purple, :green, :orange]
        for (i, groupedUnits) in enumerate(allGroupedUnits)
            color = colors[i]
            antiGroupedUnits = vcat(setdiff(allGroupedUnits, [groupedUnits])...)
            groupStrength = sum(groupedUnits[!, r]) / nrow(groupedUnits)
            antiGroupStrength = sum(antiGroupedUnits[!, r]) / nrow(antiGroupedUnits)
            if isnan(groupStrength)
                groupStrength = 0
            end
            
            if isnan(antiGroupStrength)
                antiGroupStrength = 0
            end
            
            netStrength = groupStrength - antiGroupStrength
            if netStrength > lineStrengths[r]
                lineStrengths[r] = netStrength
                lineColors[r] = color
            end
        end
    end

    ## Colors
    networkColors = map(eachrow(ena.networkModel)) do networkRow
        r = networkRow[:relationship]
        return lineColors[r]
    end

    unitColors = map(eachrow(ena.refitUnitModel)) do refitRow
        if refitRow[artist.groupVar1] == artist.controlGroup1 &&
           refitRow[artist.groupVar2] == artist.controlGroup2
            return :blue
        elseif refitRow[artist.groupVar1] == artist.controlGroup1 &&
               refitRow[artist.groupVar2] == artist.treatmentGroup2
            return :purple
        elseif refitRow[artist.groupVar1] == artist.treatmentGroup1 &&
               refitRow[artist.groupVar2] == artist.controlGroup2
            return :green
        elseif refitRow[artist.groupVar1] == artist.treatmentGroup1 &&
               refitRow[artist.groupVar2] == artist.treatmentGroup2
            return :orange
        else
            return :black
        end
    end

    codeColors = map(eachrow(ena.codeModel)) do codeRow
        return :black
    end

    ## Shapes
    unitShapes = map(eachrow(ena.refitUnitModel)) do refitRow
        return :circle
    end

    codeShapes = map(eachrow(ena.codeModel)) do codeRow
        return :circle
    end

    ## Sizes
    networkLineWidths = map(eachrow(ena.networkModel)) do networkRow
        r = networkRow[:relationship]
        return lineStrengths[r]
    end

    unitMarkerSizes = map(eachrow(ena.refitUnitModel)) do refitRow
        return 1
    end

    codeMarkerSizes = map(eachrow(ena.codeModel)) do codeRow
        return codeRow[:density]
    end

    ## Normalize
    s = maximum(networkLineWidths)
    networkLineWidths /= s
    networkLineWidths *= 2

    s = maximum(unitMarkerSizes)
    unitMarkerSizes /= s
    unitMarkerSizes *= 1.5

    s = maximum(codeMarkerSizes)
    codeMarkerSizes /= s
    codeMarkerSizes *= 8

    ## Confidence Intervals
    confidenceIntervals = []
    colors = [:blue, :purple, :green, :orange]
    for (i, groupedUnits) in enumerate(allGroupedUnits)
        color = colors[i]
        shape = :square
        size = 4
        if nrow(groupedUnits) > 1
            mu_x = mean(groupedUnits[!, :pos_x])
            mu_y = mean(groupedUnits[!, :pos_y])
            ci_x = collect(confint(OneSampleTTest(groupedUnits[!, :pos_x])))
            ci_y = collect(confint(OneSampleTTest(groupedUnits[!, :pos_y])))
            CI = (mu_x, mu_y, ci_x, ci_y, color, shape, size)
            push!(confidenceIntervals, CI)
        elseif nrow(groupedUnits) == 1
            mu_x = mean(groupedUnits[!, :pos_x])
            mu_y = mean(groupedUnits[!, :pos_y])
            CI = (mu_x, mu_y, [mu_x, mu_x], [mu_y, mu_y], color, shape, size)
            push!(confidenceIntervals, CI)
        else
            # do nothing
        end
    end

    ## Do Callback so ENADisplay can do the bulk of the work
    cb(networkColors,
       unitColors,
       codeColors,
       unitShapes,
       codeShapes,
       networkLineWidths,
       unitMarkerSizes,
       codeMarkerSizes,
       confidenceIntervals)
end

# TVRemoteArtist
function (artist::TVRemoteArtist)(cb, ena, scene)
    ##Pre-splitting the groups
    group0xUnits = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar1] == artist.controlGroup1
            return true
        else
            return false
        end
    end

    group1xUnits = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar1] == artist.treatmentGroup1
            return true
        else
            return false
        end
    end

    group00Units = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar1] == artist.controlGroup1 &&
           refitRow[artist.groupVar2] == artist.controlGroup2
            return true
        else
            return false
        end
    end

    group01Units = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar1] == artist.controlGroup1 &&
           refitRow[artist.groupVar2] == artist.treatmentGroup2
            return true
        else
            return false
        end
    end

    group10Units = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar1] == artist.treatmentGroup1 &&
           refitRow[artist.groupVar2] == artist.controlGroup2
            return true
        else
            return false
        end
    end

    group11Units = filter(ena.refitUnitModel) do refitRow
        if refitRow[artist.groupVar1] == artist.treatmentGroup1 &&
           refitRow[artist.groupVar2] == artist.treatmentGroup2
            return true
        else
            return false
        end
    end

    allGroupedUnits = [group00Units,
        group01Units,
        group10Units,
        group11Units,
        group0xUnits,
        group1xUnits]

    # ## Pre-processing the size of the network lines
    # ### Find the "direction", "strength", and "angle" for the line size
    # ### Thicker lines are those whose rotation weights are more towards
    # ### one side of the difference of the means
    # mu_x_control = mean(group0xUnits[!, :pos_x])
    # mu_y_control = mean(group0xUnits[!, :pos_y])
    # mu_x_treatment = mean(group1xUnits[!, :pos_x])
    # mu_y_treatment = mean(group1xUnits[!, :pos_y])
    # mu_x_all = (mu_x_treatment + mu_x_control) / 2
    # mu_y_all = (mu_y_treatment + mu_y_control) / 2
    # vt = Vector{Float64}([
    #     mu_x_treatment - mu_x_all,
    #     mu_y_treatment - mu_y_all
    # ])

    # norm_vt = sqrt(dot(vt, vt))
    # lineStrengths = Dict{Symbol,Float64}()
    # for networkRow in eachrow(ena.networkModel)
    #     r = networkRow[:relationship]
    #     vl = Vector{Float64}([
    #         networkRow[:weight_x],
    #         networkRow[:weight_y]
    #     ])

    #     lineStrengths[r] = dot(vl, vt) / norm_vt
    # end

    # ## Colors
    # networkColors = map(eachrow(ena.networkModel)) do networkRow
    #     if lineStrengths[networkRow[:relationship]] < 0
    #         return :purple
    #     else
    #         return :orange
    #     end
    # end

    lineStrengths0x = Dict{Symbol,Float64}()
    lineStrengths1x = Dict{Symbol,Float64}()
    for networkRow in eachrow(ena.networkModel)
        r = networkRow[:relationship]
        lineStrengths0x[r] = sum(group0xUnits[!, r] .- minimum(ena.refitUnitModel[!, r])) / nrow(group0xUnits)
        lineStrengths1x[r] = sum(group1xUnits[!, r] .- minimum(ena.refitUnitModel[!, r])) / nrow(group1xUnits)
        if isnan(lineStrengths0x[r])
            lineStrengths0x[r] = 0
        end

        if isnan(lineStrengths1x[r])
            lineStrengths1x[r] = 0
        end
    end

    ## Colors
    networkColors = map(eachrow(ena.networkModel)) do networkRow
        r = networkRow[:relationship]
        if lineStrengths0x[r] - lineStrengths1x[r] > 0
            return :purple
        elseif lineStrengths0x[r] - lineStrengths1x[r] < 0
            return :orange
        else
            return :white
        end
    end

    unitColors = map(eachrow(ena.refitUnitModel)) do refitRow
        if refitRow[artist.groupVar1] == artist.controlGroup1 &&
           refitRow[artist.groupVar2] == artist.controlGroup2
            return :purple
        elseif refitRow[artist.groupVar1] == artist.controlGroup1 &&
               refitRow[artist.groupVar2] == artist.treatmentGroup2
            return :purple
        elseif refitRow[artist.groupVar1] == artist.treatmentGroup1 &&
               refitRow[artist.groupVar2] == artist.controlGroup2
            return :orange
        elseif refitRow[artist.groupVar1] == artist.treatmentGroup1 &&
               refitRow[artist.groupVar2] == artist.treatmentGroup2
            return :orange
        else
            return :black
        end
    end

    codeColors = map(eachrow(ena.codeModel)) do codeRow
        return :black
    end

    codeColors = map(eachrow(ena.codeModel)) do codeRow
        return :black
    end

    ## Shapes
    unitShapes = map(eachrow(ena.refitUnitModel)) do refitRow
        return :circle
    end

    codeShapes = map(eachrow(ena.codeModel)) do codeRow
        return :circle
    end

    ## Sizes
    networkLineWidths = map(eachrow(ena.networkModel)) do networkRow
        # return lineStrengths[networkRow[:relationship]] ^ 2
        r = networkRow[:relationship]
        return abs(lineStrengths0x[r] - lineStrengths1x[r])
    end

    unitMarkerSizes = map(eachrow(ena.refitUnitModel)) do refitRow
        return 1
    end

    codeMarkerSizes = map(eachrow(ena.codeModel)) do codeRow
        return codeRow[:density]
    end

    ## Normalize
    s = maximum(networkLineWidths)
    networkLineWidths /= s
    networkLineWidths *= 2

    s = maximum(unitMarkerSizes)
    unitMarkerSizes /= s
    unitMarkerSizes *= 1.5

    s = maximum(codeMarkerSizes)
    codeMarkerSizes /= s
    codeMarkerSizes *= 8

    ## Confidence Intervals
    confidenceIntervals = []
    colors = [:purple, :purple, :orange, :orange, :purple, :orange]
    shapes = [:dtriangle, :utriangle, :dtriangle, :utriangle, :square, :square]
    for (i, groupedUnits) in enumerate(allGroupedUnits)
        color = colors[i]
        shape = shapes[i]
        size = 4
        if (i < 5 && nrow(groupedUnits) > 1) || # draw a quadrant CI if we have enough
            (i == 5 && nrow(group00Units) > 0 && nrow(group01Units) > 0) || # draw an omnibus CI only if each of its quadrants have something
            (i == 6 && nrow(group10Units) > 0 && nrow(group11Units) > 0)
            
            mu_x = mean(groupedUnits[!, :pos_x])
            mu_y = mean(groupedUnits[!, :pos_y])
            ci_x = collect(confint(OneSampleTTest(groupedUnits[!, :pos_x])))
            ci_y = collect(confint(OneSampleTTest(groupedUnits[!, :pos_y])))
            CI = (mu_x, mu_y, ci_x, ci_y, color, shape, size)
            push!(confidenceIntervals, CI)
        elseif (i < 5 && nrow(groupedUnits) == 1) # only for quadrant CIs
            mu_x = mean(groupedUnits[!, :pos_x])
            mu_y = mean(groupedUnits[!, :pos_y])
            CI = (mu_x, mu_y, [mu_x, mu_x], [mu_y, mu_y], color, shape, size)
            push!(confidenceIntervals, CI)
        else
            # do nothing
        end
    end

    ## Do Callback so ENADisplay can do the bulk of the work
    cb(networkColors,
       unitColors,
       codeColors,
       unitShapes,
       codeShapes,
       networkLineWidths,
       unitMarkerSizes,
       codeMarkerSizes,
       confidenceIntervals)
end