# TODO friendly error message for when they try to use an unknown name,
# tell them to just use CSV

function ena_dataset(name::String)
    filename = joinpath(dirname(@__FILE__), "..", "data", "$(name).csv")
    return DataFrame(CSV.File(filename, missingstring="#N/A"))
end