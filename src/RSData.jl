function ena_dataset(name::String)
    filename = joinpath(dirname(@__FILE__), "..", "data", "$(name).csv")
    return CSV.read(filename)
end