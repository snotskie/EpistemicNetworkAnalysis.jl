# Type Tree

## Rotations
abstract type AbstractENARotation end
abstract type AbstractLinearENARotation <: AbstractENARotation end
# abstract type AbstractNonlinearENARotation <: AbstractENARotation end

## Models
abstract type AbstractENAModel{T<:AbstractENARotation} end
abstract type AbstractLinearENAModel{T<:AbstractLinearENARotation} <: AbstractENAModel{T} end
# abstract type AbstractNonlinearENAModel{T<:AbstractNonlinearENARotation} <: AbstractENAModel{T} end

## EdgePainters
abstract type AbstractEdgePainter end

## Type Helpers
function modelsupertype(::Type{M}, ::Type{N}) where {
        R<:AbstractENARotation,
        M<:AbstractENAModel{R},
        N<:AbstractENAModel
    }
    
    return supertype(N){R}
end

function rotationsupertype(::Type{M}, ::Type{S}) where {
        R<:AbstractENARotation,
        M<:AbstractENAModel{R},
        S<:AbstractENARotation
    }
    
    return M.name.wrapper{supertype(S)}
end

## Type Macros
macro enamodel(
        self, parent,
        defaultrotation=:SVDRotation,
        rotationtype=:AbstractLinearENARotation,
        abstractname=Symbol(string(:Abstract, self))
    )

    return quote
        # make abstract type
        abstract type $(esc(abstractname)){T<:$(esc(rotationtype))} <: $(esc(parent)){T} end

        # make struct
        mutable struct $(esc(self)){T<:$(esc(rotationtype))} <: $(esc(abstractname)){T}
            # required arguments
            data::DataFrame
            codes::Array{Symbol,1}
            conversations::Array{Symbol,1}
            units::Array{Symbol,1}
        
            # rotation
            rotation::T
        
            # model
            metadata::DataFrame
            points::DataFrame
            pointsHat::DataFrame
            pointsNodes::DataFrame
            accum::DataFrame
            accumHat::DataFrame
            edges::DataFrame
            nodes::DataFrame
            embedding::DataFrame
            config::NamedTuple
        end

        # make default constructor
        function $(esc(self))(
                data::DataFrame,
                codes::Array{<:Any,1},
                conversations::Array{<:Any,1},
                units::Array{<:Any,1};
                rotateBy::R=$(esc(defaultrotation))(),
                kwargs...
            ) where {R<:$(esc(rotationtype))}

            # call common ENA constructor
            return constructENA(
                $(esc(self)){R},
                data, codes, conversations, units,
                rotateBy; kwargs...
            )
        end

        # make re-modeler constructor
        function $(esc(self))(
                prev_model::AbstractENAModel;
                rotateBy::R=prev_model.rotation, #$(esc(defaultrotation))(),
                kwargs...
            ) where {R<:$(esc(rotationtype))}

            # call common ENA re-modeler
            return remodelENA(
                $(esc(self)){R},
                prev_model,
                rotateBy;
                kwargs...
            )
        end

        # make re-rotator constructor
        function $(esc(self))(
                prev_model::$(esc(self));
                rotateBy::R=$(esc(defaultrotation))()
            ) where {R<:$(esc(rotationtype))}

            # call common ENA re-rerotator
            return rerotateENA(
                $(esc(self)){R},
                prev_model,
                rotateBy
            )
        end
    end
end

# Wrapper constructor, accepts any type for codes etc.
function constructENA(
        ::Type{M},
        data::DataFrame,
        codes::Array{<:Any,1},
        conversations::Array{<:Any,1},
        units::Array{<:Any,1},
        rotation::AbstractENARotation;
        kwargs...
    ) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    return constructENA(
        M,
        data,
        Symbol.(codes),
        Symbol.(conversations),
        Symbol.(units),
        rotation;
        kwargs...
    )
end

# Base constructor, requires symbols
function constructENA(
        ::Type{M},
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1},
        rotation::AbstractENARotation;
        kwargs...
    ) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}

    kwargs = defaultmodelkwargs(M; kwargs...)
    model = M(
        populateENAfields(
            M, data, codes, conversations, units, rotation;
            kwargs...
        )...
    )

    accumulate!(M, model)
    approximate!(M, model)
    rotate!(M, model)
    test!(M, model)
    return model
end

# Re-model constructor
function remodelENA(
        ::Type{M},
        prev_model::AbstractENAModel,
        rotation::AbstractENARotation;
        kwargs...
    ) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}

    kwargs = defaultmodelkwargs(M; prev_config=prev_model.config, kwargs...)
    model = M(
        populateENAfields(
            M,
            copy(prev_model.data),
            copy(prev_model.codes),
            copy(prev_model.conversations),
            copy(prev_model.units),
            rotation;
            kwargs...
        )...
    )

    accumulate!(M, model)
    approximate!(M, model)
    rotate!(M, model)
    test!(M, model)
    return model
end

# Re-rotate constructor
function rerotateENA(
        ::Type{M},
        prev_model::AbstractENAModel,
        rotation::AbstractENARotation
    ) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}

    model = M(
        copy(prev_model.data),
        copy(prev_model.codes),
        copy(prev_model.conversations),
        copy(prev_model.units),
        rotation,
        copy(prev_model.metadata),
        similar(prev_model.points, 0),
        similar(prev_model.pointsHat, 0),
        similar(prev_model.pointsNodes, 0),
        copy(prev_model.accum),
        copy(prev_model.accumHat),
        copy(prev_model.edges),
        copy(prev_model.nodes),
        similar(prev_model.embedding, 0),
        prev_model.config
    )

    rotate!(M, model)
    test!(M, model)
    return model
end

## Unimplemented Functions
function populateENAfields(
        ::Type{M},
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1},
        rotation::AbstractENARotation;
        config...
    ) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    
    error("Unimplemented")
end

function defaultmodelkwargs(::Type{M}, model::AbstractENAModel; kwargs...) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function defaultplotkwargs(::Type{M}, model::AbstractENAModel; kwargs...) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function defaultplotkwargs(::Type{M}, model::AbstractENAModel, config::NamedTuple) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    # wrapper for the above, to simplify logic of children that override their parent plotkwargs
    defaultplotkwargs(M, model; Dict(zip(keys(config), values(config)))...)
end

# NOTE: when implementing these functions elsewhere, M should be the *most* specific
# type that the function applies to, while model should be the *least* specific.
# Also, accumulate! and approximate! should be generic to the rotation type when possible,
# while rotate! should be generic to the ENA type when possible.
function accumulate!(::Type{M}, model::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function approximate!(::Type{M}, model::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function rotate!(::Type{M}, model::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function test!(::Type{M}, model::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    test!(M, model, model)
end

function test!(::Type{M}, trainmodel::AbstractENAModel, testmodel::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function test!(::Type{M}, trainmodel::AbstractENAModel, testmodel::AbstractENAModel, test::Type{<:HypothesisTests.HypothesisTest}; kwargs...) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

function summary(::Type{M}, model::AbstractENAModel) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    columns = [
        :label,
        setdiff(Symbol.(names(model.embedding)), [:label, model.edges.edgeID...])...
    ]

    return copy(model.embedding[!, columns])
end

"""
    summary(model::AbstractENAModel)

Produce a dataframe containing summary statistics for each dimension of the model embedding
    
## Example

```julia
model = ENAModel(data, codes, conversations, units)
stats = summary(model)
```
"""
summary

function show(io::IO, model::AbstractENAModel)
    details = (
        ENATool="EpistemicNetworkAnalysis.jl",
        ToolVersion=PkgVersion.Version(EpistemicNetworkAnalysis),
        ToolAuthor=PkgVersion.Author(EpistemicNetworkAnalysis),
        ModelConfig=(
            codes=model.codes,
            conversations=model.conversations,
            units=model.units,
            model.config...
        ),
        NumberOfUnits=nrow(model.metadata),
        RotationType=string(nameof(typeof(model.rotation))),
        RotationConfig=
            length(propertynames(model.rotation)) > 0 ?
            namedtuple(propertynames(model.rotation), fieldvalues(model.rotation)) :
            (),
        Dimensions=Tables.rowtable(summary(model))
    )

    pprint(io, details)
end

"""
    show(model::AbstractENAModel)

Display text summarizing a model's configuration and summary statistics.
    
## Example

```julia
model = ENAModel(data, codes, conversations, units)
show(model)
```
"""
show(::AbstractENAModel)

# in linear, do plot like the consruct helper, override its components under there
function plot(::Type{M}, model::AbstractENAModel, plotconfig::NamedTuple) where {R<:AbstractENARotation, M<:AbstractENAModel{R}}
    error("Unimplemented")
end

## Wrapper Functions
function summary(model::AbstractENAModel)
    return summary(typeof(model), model)
end

function plot(model::AbstractENAModel; kwargs...)
    plotconfig = defaultplotkwargs(typeof(model), model; kwargs...)
    return plot(typeof(model), model, plotconfig)
end