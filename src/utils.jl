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

struct PointCloud
    X::Matrix{Real}
    feature_names::Vector{Symbol}
    unit_names::Vector{Symbol}
    mode::Symbol
    z_normed::Bool
    z_means::Vector{Union{Real,Missing}}
    z_stds::Vector{Union{Real,Missing}}
end

function pointcloud(
        model::AbstractENAModel;
        ndims::Int=nrow(model.points),
        mode::Symbol=:wide,
        z_norm::Bool=false,
        metadata::Vector{Symbol}=Symbol[]
    )

    unit_names = model.metadata.unitID
    feature_names = [Symbol.(model.embedding.label[1:ndims])..., metadata...]
    X = Matrix(model.points[1:ndims, unit_names]) # wide format
    if !isempty(metadata)
        X = vcat(X, transpose(Matrix(model.metadata[!, metadata]))) # transpose(long format)
        if !z_norm
            @warn "When adding additional metadata to pointcloud, it is recommended to normalize values with pointcloud(model, z_norm=true, ...). See the documentation for pointcloud for more information."
        end
    end

    means = repeat(Union{Real,Missing}[missing], size(X, 1))
    stds = repeat(Union{Real,Missing}[missing], size(X, 1))
    if z_norm
        for i in axes(X, 1)
            means[i] = mean(X[i, :])
            stds[i] = std(X[i, :])
            X[i, :] .-= means[i]
            X[i, :] ./= stds[i]
        end
    end

    if mode == :tall
        X = transpose(X)
    elseif mode != :wide
        @error "pointcloud mode must be :tall or :wide"
    end

    return PointCloud(X, feature_names, unit_names, mode, z_norm, means, stds)
end

"""
    pointcloud(
        model::AbstractENAModel;
        ndims::Int=nrow(model.points),
        mode::Symbol=:wide,
        z_norm::Bool=false,
        metadata::Vector{Symbol}=Symbol[]
    )

Produce a point cloud matrix from a model's plotted points and optional additional metadata columns,
for preparing data to pass to other packages, e.g., for machine learning.

## Arguments

Required:

- `model`: The ENA model to produce a point cloud from

Optional:

- `ndims`: The number of dimensions from the ENA model's embedding to include in the point cloud. The first `ndim` dimensions will be included. By default, all dimensions will be included
- `mode`: The orientation of the point cloud, either in `:wide` format (default) or `:tall` format. In wide format, the point cloud's `X` matrix's rows will correspond to features. In tall format, they will correspond to units.
- `z_norm`: Whether to normalize the point cloud's features (default: false)
- `metadata`: A list of additional names of metadata columns from the model to include in the point cloud. Note, when including additional metadata, it is advised to also set `z_norm` to true

## Fields

Once the point cloud is constructed, it will have the following fields:

- `X`: A matrix containing the point cloud data, in either wide or tall format
- `feature_names`: A vector of the names of the features included in the point cloud. When `mode` is `:wide`, `feature_names` corresponds to the rows of `X`. When `mode` is `:tall`, it corresponds to the columns of `X` instead.
- `unit_names`: A vector of the IDs of the units included in the point cloud. When `mode` is `:wide`, `unit_names` corresponds to the columns of `X`. When `mode` is `:tall`, it corresponds to the rows of `X` instead.
- `z_normed`: A boolean representing whether the point cloud was normalized
- `z_means` and `z_stds`: When `z_normed` is true, these are vectors of the original means and standard deviations of the features of the point cloud

## Example

```julia
# Wide format DataFrame
pc = pointcloud(model)
df = DataFrame(pc.X, pc.unit_names)

# Tall format DataFrame
pc = pointcloud(model, mode=:tall)
df = DataFrame(pc.X, pc.feature_names)

# ndims, metadata, and z_norm
pc = pointcloud(model, ndims=4, mode=:tall, metadata=[:Act], z_norm=true)
df = DataFrame(pc.X, pc.feature_names)
```
"""
pointcloud