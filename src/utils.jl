function loadExample(name::String)
    if name in ["shakespeare.data", "transitions.data"]
        filename = joinpath(dirname(@__FILE__), "..", "data", "$(name).csv")
        return DataFrame(CSV.File(filename, missingstring="#N/A"))
    else
        error("loadExample only works with loadExample(\"shakespeare.data\") and loadExample(\"transitions.data\") right now. If you have your own dataset, you can load it with DataFrame(CSV.File(file_name_here, other_settings_here)). Thanks!")
    end
end

function derivedAnyCode!(data, newCol, oldCols...)
    data[!, newCol] = ones(nrow(data))
    for col in oldCols
        data[!, newCol] = data[!, newCol] .* (1 .- data[!, col])
    end

    data[!, newCol] = 1 .- data[!, newCol]
end

function derivedAllCode!(data, newCol, oldCols...)
    data[!, newCol] = ones(nrow(data))
    for col in oldCols
        data[!, newCol] = data[!, newCol] .* data[!, col]
    end
end

# TODO to_xlsx and from_xlsx