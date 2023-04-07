function loadExample(name::String)
    if name in ["shakespeare.data", "transitions.data"]
        filename = joinpath(dirname(@__FILE__), "..", "data", "$(name).csv")
        return DataFrame(CSV.File(filename, missingstring="#N/A"))
    else
        error("loadExample only works with loadExample(\"shakespeare.data\") and loadExample(\"transitions.data\") right now. If you have your own dataset, you can load it with DataFrame(CSV.File(file_name_here, other_settings_here)). Thanks!")
    end
end
