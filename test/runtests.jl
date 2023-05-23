using EpistemicNetworkAnalysis
using Test
using DataFrames
using GLM
using Plots

@testset "EpistemicNetworkAnalysis.jl" begin

    # Test that example data loads
    data = loadExample("shakespeare")
    @test nrow(data) == 1978
    @test Set(names(data)) == Set(["X","RowID","Type","Play","Act","Scene","Utterance_Index","Line_Index","Speaker","Line","Love","Beauty","Death","Fear","Friendship","Hate","Honor","Men","Women","Pride"])

    # Test that basic model runs
    conversations = [:Play, :Act]
    units = [:Play, :Speaker]
    codes = [
        :Love,
        :Death,
        :Honor,
        :Men,
        :Women
    ]

    myENA = ENAModel(data, codes, conversations, units)
    @test typeof(myENA) == ENAModel{SVDRotation}
    p = plot(myENA)
    @test typeof(p) <: Plots.Plot
    @test length(p.subplots) in [3, 4]

    # Test that model also accepts strings
    codes = string.(codes)
    conversations = string.(conversations)
    units = string.(units)
    myENA = ENAModel(data, codes, conversations, units)
    @test typeof(myENA) == ENAModel{SVDRotation}

    # Test that each model/rotation combination runs
    models = [ENAModel, BiplotENAModel]
    rotations = Dict(
        "SVD1" => SVDRotation(),
        "Act" => FormulaRotation(
            LinearModel,
            @formula(y ~ 1 + Act),
            2,
            nothing
        ),
        "Play" => MeansRotation(:Play, "Romeo and Juliet", "Hamlet"),
        "Gender" => TopicRotation("Gender", [:Women], [:Men])
    )

    for M in models
        for (label, rotation) in rotations
            myENA = M(
                data, codes, conversations, units,
                rotateBy=rotation
            )

            @test typeof(myENA) == M{typeof(rotation)}
            @test myENA.embedding[1, :label] == label

            copiedENA = M(
                data, codes, conversations, units,
                rotateBy=ManualRotation(myENA.embedding)
            )

            @test typeof(copiedENA) == M{ManualRotation}
            @test copiedENA.embedding[1, :label] == label
        end
    end

    # TODO plotting
    # TODO to/from ODS
    # TODO means rotation
    # TODO conversions
    # TODO show
    # TODO other models and rotation types
    # TODO clean up
    # TODO auto-docs
end