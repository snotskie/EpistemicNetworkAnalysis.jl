function loadExample(name::AbstractString)
    options = ["shakespeare", "transitions", "toy"]
    if name in options
        filename = joinpath(dirname(@__FILE__), "..", "data", "$(name).csv")
        return DataFrame(CSV.File(filename, missingstring="#N/A"))
    else
        error("Unrecognized dataset, $(name). Options are $(options)")
    end
end

"""
    loadExample(name::AbstractString)

Load an example dataset as a DataFrame

## Datasets

- `loadExample("shakespeare")`: Loads the [Shakespeare dataset](https://bookdown.org/tan78/intro_to_ena/examples.html), containing data on two plays, "Hamlet" and "Romeo and Juliet"
- `loadExample("transition")`: Loads the [Telling Stories of Transitions dataset](https://link.springer.com/chapter/10.1007/978-3-030-93859-8_8), containing metadata and codes only, due to the sensitive nature of the underlying text
- `loadExample("toy")`: Loads a minimal toy example, reproduced below

```csv
Group,Convo,Unit,Line,A,B,C
Red,1,X,1,0,0,1
Red,1,Y,2,1,0,0
Blue,1,Z,3,0,1,1
Blue,1,W,4,0,0,0
Red,1,X,5,0,0,1
Red,2,X,1,1,0,0
Red,2,Y,2,1,0,0
Blue,2,Z,3,0,1,1
Blue,2,W,4,0,0,0
Red,2,X,5,1,0,0
```

## Loading Your Own Data

To load your own datasets, use [`DataFrame`](https://dataframes.juliadata.org/stable/lib/types/#DataFrames.DataFrame) and [`CSV.File`](https://csv.juliadata.org/stable/reading.html#CSV.File), which requires the `DataFrames` and `CSV` packages

```julia
using Pkg
Pkg.add("DataFrames")
using DataFrames

Pkg.add("CSV")
using CSV

data = DataFrame(CSV.File("filename_here.csv"))
```
"""
loadExample

function deriveAnyCode!(data::DataFrame, newCol::Symbol, oldCols...)
    data[!, newCol] = ones(nrow(data))
    for col in oldCols
        data[!, newCol] = data[!, newCol] .* (1 .- data[!, col])
    end

    data[!, newCol] = 1 .- data[!, newCol]
end

"""
    deriveAnyCode!(
        data::DataFrame,
        newColumnName::Symbol,
        oldColumnNames...
    )

Add a new code column to `data`, derived from existing codes. The new code will be marked present where any of the old codes are present

## Example

```julia
deriveAnyCode!(data, :Food, :Hamburgers, :Salads, :Cereal)
```
"""
deriveAnyCode!

function deriveAllCode!(data::DataFrame, newCol::Symbol, oldCols...)
    data[!, newCol] = ones(nrow(data))
    for col in oldCols
        data[!, newCol] = data[!, newCol] .* data[!, col]
    end
end

"""
    deriveAllCode!(
        data::DataFrame,
        newColumnName::Symbol,
        oldColumnNames...
    )

Add a new code column to `data`, derived from existing codes. The new code will be marked present only where all of the old codes are present on the same line

## Example

```julia
deriveAllCode!(data, :ObservingStudentsLearning, :Observing, :Students, :Learning)
```
"""
deriveAllCode!

# Serialization
function to_xlsx(filename::AbstractString, model::AbstractENAModel)
    i = 1
    namesheet(xf, sheetname) = begin
        if i == 1
            XLSX.rename!(xf[1], sheetname)
        else
            XLSX.addsheet!(xf, sheetname)
        end
    end

    fixcelltype(x) = begin
        if typeof(x) <: Union{Missing, Bool, Float64, Int64, Dates.Date, Dates.DateTime, Dates.Time, String}
            return x
        elseif typeof(x) <: Function
            return join(string.(code_lowered(x)[1].code), "; ")
        else
            return string(x)
        end
    end

    fixcelltypes(col) = fixcelltype.(col)
    prepare(df) = fixcelltypes.(eachcol(df))
    putsheet(xf, df, sheetname) = begin
        namesheet(xf, sheetname)
        XLSX.writetable!(xf[i], prepare(df), names(df))
        i += 1
    end

    putvector(xf, pair, sheetname) = begin
        namesheet(xf, sheetname)
        XLSX.writetable!(xf[i], [fixcelltypes(last(pair))], [first(pair)])
        i += 1
    end

    puttuple(xf, tp, sheetname) = begin
        namesheet(xf, sheetname)
        XLSX.writetable!(xf[i], fixcelltypes.([[x] for x in values(tp)]), keys(tp))
        i += 1
    end

    putrotation(xf, rotation, sheetname) = begin
        namesheet(xf, sheetname)
        cols = ["RotationType", propertynames(rotation)...]
        vals = [string(nameof(typeof(rotation))), fixcelltypes.(fieldvalues(rotation))...]
        XLSX.writetable!(xf[i], [[x] for x in vals], cols)
        i += 1
    end

    XLSX.openxlsx(filename, mode="w") do xf
        putsheet(xf, DataFrame([
            "Data"          "The original data given to the ENA model";
            "Codes"         "A list of the column names defining the qualitative codes used in the model";
            "Conversations" "A list of the column names used to distinguish conversations";
            "Units"         "A list of the column names used to distinguish units of analysis";
            "Model"         "A row recording the configuration of the model used";
            "Rotation"      "A row recording the configuration of the rotation used";
            "Metadata"      "The metadata associated with each unit of analysis";
            "Accum"         "The normalized accumulated connection weights for each unit of analysis"
            "Points"        "The plotted unit points, one column per unit of analysis, one row per dimension. By default, the first row is the X-axis and the second row is the Y-axis";
            "Edges"         "The edges used in the model";
            "Embedding"     "The embedding arrived at by the model, one column per edge, one row per dimension. Rows are sorted in parallel to the Points sheet";
            "Nodes"         "The pseudo-accumulated weights for each node, found by fitting nodes to each column of Accum, and used to calculate PointsNodes";
            "PointsNodes"   "The plotted node points, one column per node, one row per dimension, found by Embedding * Nodes. Rows are sorted in parallel to the Points sheet";
            "AccumHat"      "The approximated unit weights, used to calculate PointsHat";
            "PointsHat"     "The approximated plotted unit points (aka 'centroids'), found by Embedding * AccumHat, used to calculate model goodness of fit"
        ], ["Sheet", "Description"]), "Documentation")
        putsheet(xf, model.data, "Data")
        putvector(xf, "Codes" => model.codes, "Codes")
        putvector(xf, "Conversations" => model.conversations, "Conversations")
        putvector(xf, "Units" => model.units, "Units")
        puttuple(xf, (
            ModelType=string(nameof(typeof(model))),
            model.config...
        ), "Model")
        putrotation(xf, model.rotation, "Rotation")
        putsheet(xf, model.metadata, "Metadata")
        putsheet(xf, model.accum, "Accum")
        putsheet(xf, model.points, "Points")
        putsheet(xf, model.edges, "Edges")
        putsheet(xf, model.embedding, "Embedding")
        putsheet(xf, model.nodes, "Nodes")
        putsheet(xf, model.pointsNodes, "PointsNodes")
        putsheet(xf, model.accumHat, "AccumHat")
        putsheet(xf, model.pointsHat, "PointsHat")
    end
end

"""
    to_xlsx(filename::AbstractString, model::AbstractENAModel)

Save a model to the disk as an Excel spreadsheet, useful for sharing results with others

See also [serialize](https://docs.julialang.org/en/v1/stdlib/Serialization/#Serialization.serialize) for saving models in a more efficient format that can be reloaded into Julia using [deserialize](https://docs.julialang.org/en/v1/stdlib/Serialization/#Serialization.deserialize)

Note: a `from_xlsx` function does not exist, but is planned. The difficulty is that Excel data is at root a human-readable string format, and some components of some models are difficult to represent reliably as human-readable strings
"""
to_xlsx

# function from_xlsx(filename::AbstractString)
#     xf = XLSX.readxlsx(filename)
#     sheetnames = XLSX.sheetnames(xf)
#     expected = [
#         "Data",
#         "Codes",
#         "Conversations",
#         "Units",
#         "Model",
#         "Rotation",
#         "Metadata",
#         "Accum",
#         "Points",
#         "Edges",
#         "Embedding",
#         "Nodes",
#         "PointsNodes",
#         "AccumHat",
#         "PointsHat"
#     ]

#     for sheetname in expected
#         @assert(sheetname in sheetnames,
#             "Expected sheet named '$sheetname' not found in the spreadsheet: $filename")
#     end

#     # Helper
#     df(sheetname) = DataFrame(XLSX.gettable(xf[sheetname]))

#     # Basic config
#     data = df("Data")
#     conversations = df("Conversations").Conversations
#     units = df("Units").Units
#     codes = df("Codes").Codes

#     # Rotation
#     rotationTypeName = df("Rotation").RotationType[1]
#     rotationType = try
#         getfield(@__MODULE__, Symbol(rotationTypeName))
#     catch UndefVarError
#         @error "Unknown rotation type: $rotationTypeName"
#     end

#     @assert(rotationType <: AbstractENARotation,
#         "Invalid rotation type: $rotationTypeName")
    
#     rotationArgs = df("Rotation")[1, 2:end]
#     rotation = rotationType(rotationArgs...)
#     # This is the tricky part, zero-knowledge type conversion
#     # from whatever type XLSX and DataFrame comes up with,
#     # to whatever type (or array type) the rotation requires
    
#     # Module components
#     # Should be mostly just loading in dataframes

#     # Module config
#     # Same trickiness as the rotation

#     # Model re-construction
#     # The idea is to just call the underlying struct directly

#     return xf # TODO return the reconstructed model
# end