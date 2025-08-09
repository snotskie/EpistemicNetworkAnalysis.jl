# DONE

# use macro helper to define a standard ENA struct with all the bells
@enamodel CodewiseENAModel AbstractBiplotENAModel

# Documentation
"""
    CodewiseENAModel(
        # Required
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1};

        # Optional
        rotation::AbstractLinearENARotation=SVDRotation(),
        unitFilter::Function=unit->row[:CodewiseCode] != "__prefix__", # fixed, cannot change
        edgeFilter::Function=edge->edge.kind == :count, # fixed, cannot change
        windowSize::Real=1,
        sphereNormalize::Bool=true,
        lineNormalize::Bool=false,
        dropEmpty::Bool=false,
        recenterEmpty::Bool=false
    )

Construct a biplot model of code-wise counts of code co-occurences. Model will have perfect goodness of fit between `points` and `pointsHat`, will be much simpler than other model types, but unlike `BiplotENAModel` will not lose most information.

`CodewiseENAModel` follows the same argument and field structure as `ENAModel`, except `unitFilter` and `edgeFilter` are in effect ignored. Note also that the default `windowSize` is `1` instead of `Inf`, to prevent accidentally overloading RAM.

Two columns, `:CodewiseChorus` and `:CodewiseCode` will be added to `data`, and at least one of these should be included in the `units` parameter.
"""
CodewiseENAModel

# override default model constructor kwargs
function defaultmodelkwargs(
        ::Type{M};
        prev_config::NamedTuple=NamedTuple(),
        kwargs...
    ) where {R<:AbstractLinearENARotation, M<:AbstractCodewiseENAModel{R}}

    kwargs = NamedTuple(kwargs)
    super = modelsupertype(M, AbstractCodewiseENAModel)
    parentdefaults = defaultmodelkwargs(super)
    definitivedefaults = (
        edgeFilter=(row)->(
            row[:kind] == :count
        ),
        unitFilter=(row)->(
            row[:CodewiseCode] != "__prefix__"
        ),
        windowSize=get(prev_config, :windowSize, 1) # allow overriding what BiplotENAModel does
    )

    return merge(parentdefaults, prev_config, definitivedefaults, kwargs)
end

# override default model field populator
function populateENAfields(
        ::Type{M},
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1},
        rotation::AbstractLinearENARotation;
        config...
    ) where {R<:AbstractLinearENARotation, M<:AbstractCodewiseENAModel{R}}

    # Step 1. add columns that mark the presense of each chorus of each code
    data = copy(data)
    choruses = Symbol[] # in format of code + strike number
    current_convo = data[1, conversations]
    for code in codes
        strike = 0 # which strike number are we on
        sound = 0 # is the current strike still making sound
        chorus = nothing # temp
        for (i, val) in enumerate(data[!, code])
            if data[i, conversations] != current_convo # new convo, reset sound
                sound = 0
            end

            if sound == 0 && val == 1 # if good time to strike
                strike += 1
                chorus = Symbol(string(code, strike))
                push!(choruses, chorus)
                data[!, chorus] .= 0
            end

            if val == 1 # if mid-chorus, sustain the sound
                sound = config[:windowSize]
            elseif sound > 0 # else let it decrease
                sound -= 1
            end

            if sound > 0 # so long as still making sound, mark our presence
                data[i, chorus] = 1
            end
        end
    end

    # Step 2. Segment the data into each chorus and its prefix
    codewise_data = nothing # temp
    for chorus in choruses
        code = Symbol(replace(string(chorus), r"[0-9]+$" => ""))
        rows = data[!, chorus] .== 1
        subset = copy(data[rows, :])
        subset[!, :CodewiseChorus] .= chorus
        subset[!, :CodewiseCode] .= code
        start = findfirst(==(1), data[!, chorus])
        start -= config[:windowSize]
        if start >= 1
            prefix = copy(data[start:start+config[:windowSize], :])
            prefix[!, :CodewiseChorus] .= chorus
            prefix[!, :CodewiseCode] .= "__prefix__"
            subset = vcat(prefix, subset)
        end

        if isnothing(codewise_data)
            codewise_data = subset
        else
            codewise_data = vcat(codewise_data, subset)
        end
    end

    # Step 3. Let parent do the rest of the work
    codewise_convos = [:CodewiseChorus]
    super = modelsupertype(M, AbstractCodewiseENAModel)
    return populateENAfields(super, codewise_data, codes, codewise_convos, units, rotation; config...)
end
