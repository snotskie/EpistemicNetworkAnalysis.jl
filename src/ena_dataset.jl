function ena_dataset(name::String)
    if name in ["shakespeare.data", "transitions.data"]
        filename = joinpath(dirname(@__FILE__), "..", "data", "$(name).csv")
        return DataFrame(CSV.File(filename, missingstring="#N/A"))
    else
        error("ena_dataset only works with ena_dataset(\"shakespeare.data\") and ena_dataset(\"transitions.data\") right now. If you have your own dataset, you can load it with DataFrame(CSV.File(file_name_here, other_settings_here)). Thanks!")
    end
end
