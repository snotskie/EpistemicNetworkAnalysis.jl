# DONE

# use macro helper to define a standard ENA struct with all the bells
@enamodel DigraphENAModel AbstractLinearENAModel

# Documentation
"""
    DigraphENAModel(
        # Required
        data::DataFrame,
        codes::Array{Symbol,1},
        conversations::Array{Symbol,1},
        units::Array{Symbol,1};

        # Optional
        rotation::AbstractLinearENARotation=SVDRotation(),
        unitFilter::Function=unit->true,
        edgeFilter::Function=edge->edge.kind == :directed,
        windowSize::Real=Inf,
        sphereNormalize::Bool=true,
        lineNormalize::Bool=false,
        dropEmpty::Bool=false,
        recenterEmpty::Bool=false
    )

Construct a directed ENA model. Nodes are positioned to maximize goodness of fit between plotted points and units' weighted average of edge vectors.

`DigraphENAModel` follows the same argument and field structure as `ENAModel`.

Ensure that `edgeFilter` only includes `:directed` edges.
"""
DigraphENAModel

# override default model constructor kwargs
function defaultmodelkwargs(
        ::Type{M};
        prev_config::NamedTuple=NamedTuple(),
        kwargs...
    ) where {R<:AbstractLinearENARotation, M<:AbstractDigraphENAModel{R}}

    kwargs = NamedTuple(kwargs)
    super = modelsupertype(M, AbstractDigraphENAModel)
    parentdefaults = defaultmodelkwargs(super)
    definitivedefaults = (
        edgeFilter=(row)->(
            row[:kind] == :directed
        ),
        windowSize=1
    )

    return merge(parentdefaults, prev_config, definitivedefaults, kwargs)
end

function substantiate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLinearENARotation, M<:AbstractDigraphENAModel{R}}

    # Regression model for placing the code dots into the approximated high-dimensional space
    X = Matrix{Float64}(zeros(nrow(model.accum), nrow(model.nodes)))
    nodeIndexMap = Dict(
        nodeID => i
        for (i, nodeID) in enumerate(model.nodes.nodeID)
    )

    ## for each unit's edges, split the edge between nodes on either end
    for (k, unitRow) in enumerate(eachrow(model.accum))
        for edge in eachrow(model.edges)
            i, j = nodeIndexMap[edge.ground], nodeIndexMap[edge.response]
            X[k, i] -= unitRow[edge.edgeID]
            X[k, j] += unitRow[edge.edgeID]
        end
    end

    ## run partial regression
    X = (transpose(X) * X)^-1 * transpose(X)

    ## fit each dimension (edge) of the original high-dimensional space
    for edge in model.edges.edgeID
        y = Vector{Float64}(model.accum[:, edge])
        coefs = X * y # regress on this edge
        model.nodes[!, edge] = coefs[1:end] .- mean(coefs) # NOTE in dENA, these MUST be mean centered, else the plot will be illegible
    end
end

function approximate!(
        ::Type{M}, model::AbstractLinearENAModel
    ) where {R<:AbstractLinearENARotation, M<:AbstractDigraphENAModel{R}}

    # find accumHat
    ## Refit the units: in high-d space, the refit units are as close as possible to their
    ## center of mass wrt the network
    ## This is by-definition the refit space
    nodeIndexMap = Dict(
        nodeID => i
        for (i, nodeID) in enumerate(model.nodes.nodeID)
    )
    
    for unitEdgeID in model.edges.edgeID
        for (k, unitRow) in enumerate(eachrow(model.accumHat))
            unitRow[unitEdgeID] =
                sum(
                    model.nodes[nodeIndexMap[nodeEdge.response], unitEdgeID] * model.accum[k, nodeEdge.edgeID] -
                    model.nodes[nodeIndexMap[nodeEdge.ground  ], unitEdgeID] * model.accum[k, nodeEdge.edgeID]
                    for nodeEdge in eachrow(model.edges)
                )
        end
    end
end