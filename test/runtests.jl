using EpistemicNetworkAnalysis
using Test
using DataFrames

@testset "EpistemicNetworkAnalysis.jl" begin

    # Test that example data loads
    data = loadExample("shakespeare.data")
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

    # Test that model also accepts strings
    codes = string.(codes)
    conversations = string.(conversations)
    units = string.(units)
    myENA = ENAModel(data, codes, conversations, units)
    @test typeof(myENA) == ENAModel{SVDRotation}

    # Test that each model/rotation combination runs
    models = [ENAModel, BiplotENAModel]
    rotations = [SVDRotation(), MeansRotation(:Play, "Romeo and Juliet", "Hamlet")]
    for M in models
        for rotation in rotations
            myENA = M(
                data, codes, conversations, units,
                rotateBy=rotation
            )

            @test typeof(myENA) == M{typeof(rotation)}
        end
    end

    # TODO means rotation
    # TODO plotting
    # TODO conversions
    # TODO other models and rotation types
    # TODO clean up
    # TODO auto-docs
end