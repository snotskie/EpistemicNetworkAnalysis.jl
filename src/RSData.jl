function ena_dataset(name::String)
    filename = joinpath(dirname(@__FILE__), "..", "data", "$(name).csv")
    return DataFrame(CSV.File(filename, missingstring="#N/A"))
end