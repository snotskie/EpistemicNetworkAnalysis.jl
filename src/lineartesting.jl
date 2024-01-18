function test!(
        ::Type{M},
        ::AbstractLinearENAModel, # ignores the trainmodel, since we're just reporting summary statistics
        model::AbstractLinearENAModel
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}
    
    edgeIDs = model.edges.edgeID
    unitIDs = model.accum.unitID
    total_variance = sum(var.(eachcol(model.accum[!, edgeIDs])))
    for i in 1:nrow(model.embedding)
        points = Vector(model.points[i, unitIDs])
        pointsHat = Vector(model.pointsHat[i, unitIDs])
        model.embedding[i, :variance_explained] = var(points) / total_variance
        # model.embedding[i, :pearson] = cor(points, pointsHat)
        pointsDiffs = [
            a - b
            for a in points
            for b in points
        ]

        pointsHatDiffs = [
            a - b
            for a in pointsHat
            for b in pointsHat
        ]

        model.embedding[i, :coregistration] = cor(pointsDiffs, pointsHatDiffs)
    end

    # NOTE this fails for CopyRotation
    # @assert sum(model.embedding.variance_explained) ≈ 1.0 "Var Exp does not add up to 100% as expected"
end

struct GroupwiseCoregistrationTest end # empty type for specifying when to use the groupwise method below
function test!(
        ::Type{M},
        ::AbstractLinearENAModel, # ignores the trainmodel, since we're just reporting summary statistics
        model::AbstractLinearENAModel,
        test::Type{GroupwiseCoregistrationTest};
        dim::Int=1,
        groupVar::Symbol=:Missing,
        groups::Array=[],
        kwargs...
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    if !(:coregistration_GroupBy in Symbol.(names(model.embedding)))
        model.embedding[!, :coregistration_GroupBy] = repeat(Union{Symbol,Missing}[missing], nrow(model.embedding))
        for (i, g) in enumerate(groups)
            model.embedding[!, Symbol("coregistration_Name_$(i)")] = repeat(Union{String,Missing}[missing], nrow(model.embedding))
            model.embedding[!, Symbol("coregistration_$(i)")] = repeat(Union{Float64,Missing}[missing], nrow(model.embedding))
        end
    end

    model.embedding[dim, :coregistration_GroupBy] = groupVar
    for (i, g) in enumerate(groups)
        groupIDs = model.metadata[model.metadata[!, groupVar] .== g, :unitID]
        groupPoints = Vector(model.points[dim, groupIDs])
        groupPointsHat = Vector(model.pointsHat[dim, groupIDs])
        groupPointsDiffs = [
            a - b
            for a in groupPoints
            for b in groupPoints
        ]

        groupPointsHatDiffs = [
            a - b
            for a in groupPointsHat
            for b in groupPointsHat
        ]

        model.embedding[dim, Symbol("coregistration_$(i)")] = cor(groupPointsDiffs, groupPointsHatDiffs)
        model.embedding[dim, Symbol("coregistration_Name_$(i)")] = string(g)
    end
end

function test!(
        ::Type{M},
        ::AbstractLinearENAModel, # ignores the trainmodel, since all the train info we need is given in the kwargs
        model::AbstractLinearENAModel,
        test::Type{KruskalWallisTest};
        dim::Int=1,
        groupVar::Symbol=:Missing,
        groups::Array=[],
        kwargs...
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}
    
    mdn_cols = [Symbol("KruskalWallis_Median_$(i)") for i in 1:length(groups)]
    name_cols = [Symbol("KruskalWallis_Name_$(i)") for i in 1:length(groups)]
    n_cols = [Symbol("KruskalWallis_N_$(i)") for i in 1:length(groups)]

    if !(:KruskalWallis_H in Symbol.(names(model.embedding)))
        filler = repeat(Union{Float64,Missing}[missing], nrow(model.embedding))
        filler2 = repeat(Union{Symbol,Missing}[missing], nrow(model.embedding))
        filler3 = repeat(Any[missing], nrow(model.embedding))
        model.embedding[!, :KruskalWallis_H] = copy(filler)
        model.embedding[!, :KruskalWallis_p] = copy(filler)
        model.embedding[!, :KruskalWallis_GroupBy] = copy(filler2)
        for mdn_col in mdn_cols
            model.embedding[!, mdn_col] = copy(filler)
        end

        for name_col in name_cols
            model.embedding[!, name_col] = copy(filler3)
        end

        for n_col in n_cols
            model.embedding[!, n_col] = copy(filler)
        end
    end

    kw_values = [
        Vector(model.points[dim,
            model.accum[
                model.metadata[!, groupVar] .== g,
                :unitID
            ]
        ])
        for g in groups
    ]

    kw_test = KruskalWallisTest(kw_values...)
    model.embedding[dim, :KruskalWallis_H] = kw_test.H
    model.embedding[dim, :KruskalWallis_p] = pvalue(kw_test)
    model.embedding[dim, :KruskalWallis_GroupBy] = groupVar
    model.embedding[dim, mdn_cols] .= median.(kw_values)
    model.embedding[dim, name_cols] .= groups
    model.embedding[dim, n_cols] .= length.(kw_values)
end

function test!(
        ::Type{M},
        ::AbstractLinearENAModel, # ignores the trainmodel, since all the train info we need is given in the kwargs
        model::AbstractLinearENAModel,
        test::Type{HypothesisTests.VarianceEqualityTest}; # BrownForsythe wraps this type
        dim::Int=1,
        groupVar::Symbol=:Missing,
        groups::Array=[],
        kwargs...
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}
    
    mdn_cols = [Symbol("BrownForsythe_Median_$(i)") for i in 1:length(groups)]
    name_cols = [Symbol("BrownForsythe_Name_$(i)") for i in 1:length(groups)]
    n_cols = [Symbol("BrownForsythe_N_$(i)") for i in 1:length(groups)]

    if !(:BrownForsythe_W in Symbol.(names(model.embedding)))
        filler = repeat(Union{Float64,Missing}[missing], nrow(model.embedding))
        filler2 = repeat(Union{Symbol,Missing}[missing], nrow(model.embedding))
        filler3 = repeat(Any[missing], nrow(model.embedding))
        model.embedding[!, :BrownForsythe_W] = copy(filler)
        model.embedding[!, :BrownForsythe_p] = copy(filler)
        model.embedding[!, :BrownForsythe_GroupBy] = copy(filler2)
        for mdn_col in mdn_cols
            model.embedding[!, mdn_col] = copy(filler)
        end

        for name_col in name_cols
            model.embedding[!, name_col] = copy(filler3)
        end

        for n_col in n_cols
            model.embedding[!, n_col] = copy(filler)
        end
    end

    bf_values = [
        Vector(model.points[dim,
            model.accum[
                model.metadata[!, groupVar] .== g,
                :unitID
            ]
        ])
        for g in groups
    ]

    bf_test = BrownForsytheTest(bf_values...)
    model.embedding[dim, :BrownForsythe_W] = HypothesisTests.teststatistic(bf_test)
    model.embedding[dim, :BrownForsythe_p] = pvalue(bf_test)
    model.embedding[dim, :BrownForsythe_GroupBy] = groupVar
    model.embedding[dim, mdn_cols] .= median.(bf_values)
    model.embedding[dim, name_cols] .= groups
    model.embedding[dim, n_cols] .= length.(bf_values)
end

function test!(
        ::Type{M},
        ::AbstractLinearENAModel, # ignores the trainmodel, since all the info we need is given in the kwargs
        model::AbstractLinearENAModel,
        test::Type{<:RegressionModel};
        dim::Int=1,
        formula::FormulaTerm=@formula(y ~ 1),
        contrasts::Union{Nothing,Dict},
        kwargs...
    ) where {R<:AbstractLinearENARotation, M<:AbstractLinearENAModel{R}}

    if !(:Formula_RHS in Symbol.(names(model.embedding)))
        filler = repeat(Union{Float64,Missing}[missing], nrow(model.embedding))
        filler2 = repeat(Union{String,Missing}[missing], nrow(model.embedding))
        model.embedding[!, :Formula_R2] = copy(filler)
        model.embedding[!, :Formula_AdjR2] = copy(filler)
        model.embedding[!, :Formula_RHS] = copy(filler2)
    end

    regressionData = copy(model.metadata)
    regressionData[!, :pos_x] = Vector(model.points[dim, regressionData.unitID])
    f1 = FormulaTerm(term(:pos_x), formula.rhs)
    if contrasts isa Nothing
        m1 = fit(test, f1, regressionData)
        model.embedding[dim, :Formula_R2] = r2(m1)
        model.embedding[dim, :Formula_AdjR2] = adjr2(m1)
    else
        m1 = fit(test, f1, regressionData, contrasts=contrasts)
        model.embedding[dim, :Formula_R2] = r2(m1)
        model.embedding[dim, :Formula_AdjR2] = adjr2(m1)
    end

    model.embedding[dim, :Formula_RHS] = string(formula.rhs)
end