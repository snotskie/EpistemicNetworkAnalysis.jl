"""
TODO: document
"""
function ena_dataset(name::String)
    if name in ["RS.data"]
        filename = joinpath(dirname(@__FILE__), "..", "data", "$(name).csv")
        return DataFrame(CSV.File(filename, missingstring="#N/A"))
    else
        error("ena_dataset only works with ena_dataset(\"RS.data\") right now. If you have your own dataset, you can load it with DataFrame(CSV.File(file_name_here, other_settings_here)). Thanks!")
    end
end