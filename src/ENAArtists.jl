# TODO show one ENAUnit's model

abstract type ENAArtist
end

struct DefaultArtist <: ENAArtist
end

struct MeansArtist <: ENAArtist
    groupVar::Symbol
    controlGroup::Any
    treatmentGroup::Any
end

struct WindowsArtist <: ENAArtist
    groupVar1::Symbol
    controlGroup1::Any
    treatmentGroup1::Any
    groupVar2::Symbol
    controlGroup2::Any
    treatmentGroup2::Any
end

struct TVRemoteArtist <: ENAArtist
    groupVar1::Symbol
    controlGroup1::Any
    treatmentGroup1::Any
    groupVar2::Symbol
    controlGroup2::Any
    treatmentGroup2::Any
end

# DefaultArtist: TODO document
function (artist::DefaultArtist)(cb, ena, scene)
    ## Colors
    networkColors = map(eachrow(ena.networkModel)) do networkRow
        return :black
    end

    unitColors = map(eachrow(ena.unitModel)) do unitRow
        return :black
    end

    codeColors = map(eachrow(ena.codeModel)) do codeRow
        return :black
    end

    ## Shapes
    unitShapes = map(eachrow(ena.unitModel)) do unitRow
        return "o"
    end

    codeShapes = map(eachrow(ena.codeModel)) do codeRow
        return "o"
    end

    ## Sizes
    networkLineWidths = map(eachrow(ena.networkModel)) do networkRow
        return networkRow[:density]
    end

    unitMarkerSizes = map(eachrow(ena.unitModel)) do unitRow
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
    unitMarkerSizes *= 0.05

    s = maximum(codeMarkerSizes)
    codeMarkerSizes /= s
    codeMarkerSizes *= 0.1

    ## Confidence Intervals
    confidenceIntervals = []
    mu_x = mean(ena.unitModel[!, :fit_x])
    mu_y = mean(ena.unitModel[!, :fit_y])
    ci_x = collect(confint(OneSampleTTest(ena.unitModel[!, :fit_x])))
    ci_y = collect(confint(OneSampleTTest(ena.unitModel[!, :fit_y])))
    color = :black
    shape = :square
    size = 0.05
    CI = (mu_x, mu_y, ci_x, ci_y, color, shape, size)
    push!(confidenceIntervals, CI)

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

# MeansArtist: TODO document
function (artist::MeansArtist)(cb, ena, scene)
    ##Pre-splitting the groups
    controlUnits = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar] == artist.controlGroup
            return true
        else
            return false
        end
    end

    treatmentUnits = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar] == artist.treatmentGroup
            return true
        else
            return false
        end
    end

    ## Pre-processing the size of the network lines
    ### Find the "direction", "strength", and "angle" for the line size
    ### Thicker lines are those whose rotation weights are more towards
    ### one side of the difference of the means
    mu_x_control = mean(controlUnits[!, :fit_x])
    mu_y_control = mean(controlUnits[!, :fit_y])
    mu_x_treatment = mean(treatmentUnits[!, :fit_x])
    mu_y_treatment = mean(treatmentUnits[!, :fit_y])
    mu_x_all = (mu_x_treatment + mu_x_control) / 2
    mu_y_all = (mu_y_treatment + mu_y_control) / 2
    vt = Vector{Float64}([
        mu_x_treatment - mu_x_all,
        mu_y_treatment - mu_y_all
    ])

    norm_vt = sqrt(dot(vt, vt))
    lineStrengths = Dict{Symbol,Float64}()
    for networkRow in eachrow(ena.networkModel)
        r = networkRow[:relationship]
        vl = Vector{Float64}([
            networkRow[:weight_x],
            networkRow[:weight_y]
        ])

        lineStrengths[r] = dot(vl, vt) / norm_vt
    end

    ## Colors
    networkColors = map(eachrow(ena.networkModel)) do networkRow
        if lineStrengths[networkRow[:relationship]] < 0
            return :purple
        else
            return :orange
        end
    end

    unitColors = map(eachrow(ena.unitModel)) do unitRow
        if unitRow[artist.groupVar] == artist.controlGroup
            return :purple
        elseif unitRow[artist.groupVar] == artist.treatmentGroup
            return :orange
        else
            return :black
        end
    end

    codeColors = map(eachrow(ena.codeModel)) do codeRow
        return :black
    end

    ## Shapes
    unitShapes = map(eachrow(ena.unitModel)) do unitRow
        return "o"
    end

    codeShapes = map(eachrow(ena.codeModel)) do codeRow
        return "o"
    end

    ## Sizes
    networkLineWidths = map(eachrow(ena.networkModel)) do networkRow
        return lineStrengths[networkRow[:relationship]] ^ 2
    end

    unitMarkerSizes = map(eachrow(ena.unitModel)) do unitRow
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
    # unitMarkerSizes *= 0.05

    s = maximum(codeMarkerSizes)
    codeMarkerSizes /= s
    # codeMarkerSizes *= 0.1

    ## Confidence Intervals
    confidenceIntervals = []

    ### Control
    mu_x = mean(controlUnits[!, :fit_x])
    mu_y = mean(controlUnits[!, :fit_y])
    ci_x = collect(confint(OneSampleTTest(controlUnits[!, :fit_x])))
    ci_y = collect(confint(OneSampleTTest(controlUnits[!, :fit_y])))
    color = :purple
    shape = :square
    size = 0.05
    CI = (mu_x, mu_y, ci_x, ci_y, color, shape, size)
    push!(confidenceIntervals, CI)

    ### Treatment
    mu_x = mean(treatmentUnits[!, :fit_x])
    mu_y = mean(treatmentUnits[!, :fit_y])
    ci_x = collect(confint(OneSampleTTest(treatmentUnits[!, :fit_x])))
    ci_y = collect(confint(OneSampleTTest(treatmentUnits[!, :fit_y])))
    color = :orange
    shape = :square
    size = 0.05
    CI = (mu_x, mu_y, ci_x, ci_y, color, shape, size)
    push!(confidenceIntervals, CI)

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

# WindowsArtist: TODO document
function (artist::WindowsArtist)(cb, ena, scene)
    # ##Pre-splitting the groups
    group00Units = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar1] == artist.controlGroup1 &&
           unitRow[artist.groupVar2] == artist.controlGroup2
            return true
        else
            return false
        end
    end

    group01Units = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar1] == artist.controlGroup1 &&
           unitRow[artist.groupVar2] == artist.treatmentGroup2
            return true
        else
            return false
        end
    end

    group10Units = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar1] == artist.treatmentGroup1 &&
           unitRow[artist.groupVar2] == artist.controlGroup2
            return true
        else
            return false
        end
    end

    group11Units = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar1] == artist.treatmentGroup1 &&
           unitRow[artist.groupVar2] == artist.treatmentGroup2
            return true
        else
            return false
        end
    end

    allGroupedUnits = [group00Units,
        group01Units,
        group10Units,
        group11Units]

    ## Colors
    networkColors = map(eachrow(ena.networkModel)) do networkRow
        return :black
    end

    unitColors = map(eachrow(ena.unitModel)) do unitRow
        if unitRow[artist.groupVar1] == artist.controlGroup1 &&
           unitRow[artist.groupVar2] == artist.controlGroup2
            return :blue
        elseif unitRow[artist.groupVar1] == artist.controlGroup1 &&
               unitRow[artist.groupVar2] == artist.treatmentGroup2
            return :purple
        elseif unitRow[artist.groupVar1] == artist.treatmentGroup1 &&
               unitRow[artist.groupVar2] == artist.controlGroup2
            return :green
        elseif unitRow[artist.groupVar1] == artist.treatmentGroup1 &&
               unitRow[artist.groupVar2] == artist.treatmentGroup2
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
    unitShapes = map(eachrow(ena.unitModel)) do unitRow
        return "o"
    end

    codeShapes = map(eachrow(ena.codeModel)) do codeRow
        return "o"
    end

    ## Sizes
    networkLineWidths = map(eachrow(ena.networkModel)) do networkRow
        return networkRow[:density]
    end

    unitMarkerSizes = map(eachrow(ena.unitModel)) do unitRow
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
    unitMarkerSizes *= 0.05

    s = maximum(codeMarkerSizes)
    codeMarkerSizes /= s
    codeMarkerSizes *= 0.1

    ## Confidence Intervals
    confidenceIntervals = []
    colors = [:blue, :purple, :green, :orange]
    
    for (i, groupedUnits) in enumerate(allGroupedUnits)
        color = colors[i]
        mu_x = mean(groupedUnits[!, :fit_x])
        mu_y = mean(groupedUnits[!, :fit_y])
        ci_x = collect(confint(OneSampleTTest(groupedUnits[!, :fit_x])))
        ci_y = collect(confint(OneSampleTTest(groupedUnits[!, :fit_y])))
        shape = :square
        size = 0.05
        CI = (mu_x, mu_y, ci_x, ci_y, color, shape, size)
        push!(confidenceIntervals, CI)
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

# TVRemoteArtist: TODO document
function (artist::TVRemoteArtist)(cb, ena, scene)
    ##Pre-splitting the groups
    group0xUnits = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar1] == artist.controlGroup1
            return true
        else
            return false
        end
    end

    group1xUnits = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar1] == artist.treatmentGroup1
            return true
        else
            return false
        end
    end

    group00Units = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar1] == artist.controlGroup1 &&
           unitRow[artist.groupVar2] == artist.controlGroup2
            return true
        else
            return false
        end
    end

    group01Units = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar1] == artist.controlGroup1 &&
           unitRow[artist.groupVar2] == artist.treatmentGroup2
            return true
        else
            return false
        end
    end

    group10Units = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar1] == artist.treatmentGroup1 &&
           unitRow[artist.groupVar2] == artist.controlGroup2
            return true
        else
            return false
        end
    end

    group11Units = filter(ena.unitModel) do unitRow
        if unitRow[artist.groupVar1] == artist.treatmentGroup1 &&
           unitRow[artist.groupVar2] == artist.treatmentGroup2
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

    ## Pre-processing the size of the network lines
    ### Find the "direction", "strength", and "angle" for the line size
    ### Thicker lines are those whose rotation weights are more towards
    ### one side of the difference of the means
    mu_x_control = mean(group0xUnits[!, :fit_x])
    mu_y_control = mean(group0xUnits[!, :fit_y])
    mu_x_treatment = mean(group1xUnits[!, :fit_x])
    mu_y_treatment = mean(group1xUnits[!, :fit_y])
    mu_x_all = (mu_x_treatment + mu_x_control) / 2
    mu_y_all = (mu_y_treatment + mu_y_control) / 2
    vt = Vector{Float64}([
        mu_x_treatment - mu_x_all,
        mu_y_treatment - mu_y_all
    ])

    norm_vt = sqrt(dot(vt, vt))
    lineStrengths = Dict{Symbol,Float64}()
    for networkRow in eachrow(ena.networkModel)
        r = networkRow[:relationship]
        vl = Vector{Float64}([
            networkRow[:weight_x],
            networkRow[:weight_y]
        ])

        lineStrengths[r] = dot(vl, vt) / norm_vt
    end

    ## Colors
    networkColors = map(eachrow(ena.networkModel)) do networkRow
        if lineStrengths[networkRow[:relationship]] < 0
            return :purple
        else
            return :orange
        end
    end

    unitColors = map(eachrow(ena.unitModel)) do unitRow
        if unitRow[artist.groupVar1] == artist.controlGroup1 &&
           unitRow[artist.groupVar2] == artist.controlGroup2
            return :purple
        elseif unitRow[artist.groupVar1] == artist.controlGroup1 &&
               unitRow[artist.groupVar2] == artist.treatmentGroup2
            return :purple
        elseif unitRow[artist.groupVar1] == artist.treatmentGroup1 &&
               unitRow[artist.groupVar2] == artist.controlGroup2
            return :orange
        elseif unitRow[artist.groupVar1] == artist.treatmentGroup1 &&
               unitRow[artist.groupVar2] == artist.treatmentGroup2
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
    unitShapes = map(eachrow(ena.unitModel)) do unitRow
        return :circle
    end

    codeShapes = map(eachrow(ena.codeModel)) do codeRow
        return :circle
    end

    ## Sizes
    networkLineWidths = map(eachrow(ena.networkModel)) do networkRow
        return lineStrengths[networkRow[:relationship]] ^ 2
    end

    unitMarkerSizes = map(eachrow(ena.unitModel)) do unitRow
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
    unitMarkerSizes *= 4

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
        mu_x = mean(groupedUnits[!, :fit_x])
        mu_y = mean(groupedUnits[!, :fit_y])
        ci_x = collect(confint(OneSampleTTest(groupedUnits[!, :fit_x])))
        ci_y = collect(confint(OneSampleTTest(groupedUnits[!, :fit_y])))
        size = 0.05
        CI = (mu_x, mu_y, ci_x, ci_y, color, shape, size)
        push!(confidenceIntervals, CI)
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