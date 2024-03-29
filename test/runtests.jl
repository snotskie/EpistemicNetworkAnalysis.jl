using EpistemicNetworkAnalysis
using Test
using DataFrames
using GLM
using Plots
using LinearAlgebra

ENV["JULIA_DEBUG"] = EpistemicNetworkAnalysis

data = loadExample("shakespeare")
@testset "Data loads" begin
    @test nrow(data) == 1978
    @test Set(names(data)) == Set(["X","RowID","Type","Play","Act","Scene","Utterance_Index","Line_Index","Speaker","Line","Love","Beauty","Death","Fear","Friendship","Hate","Honor","Men","Women","Pride"])
end

conversations = [:Play, :Act]
units = [:Play, :Speaker]
codes = [
    :Love,
    :Death,
    :Honor,
    :Men,
    :Women
]

@testset "Basic model runs and plots" begin
    myENA = ENAModel(data, codes, conversations, units)
    p = plot(myENA)
    @test typeof(myENA) == ENAModel{SVDRotation}
    @test typeof(p) <: Plots.Plot
    @test length(p.subplots) in [3, 4]
end

@testset "Model accepts strings" begin
    myENA = ENAModel(data, string.(codes), string.(conversations), string.(units))
    @test typeof(myENA) == ENAModel{SVDRotation}
end

# Test that each model/rotation combination runs
models = [ENAModel, BiplotENAModel, DigraphENAModel]
rotations = [
    # expected xaxis label => rotation
    "SVD1" => SVDRotation(),
    "Act" => FormulaRotation(
        LinearModel,
        @formula(y ~ 1 + Act),
        2,
        nothing
    ),
    "Play" => MeansRotation(:Play, "Romeo and Juliet", "Hamlet"),
    "Gender" => TopicRotation("Gender", [:Women], [:Men]),
    "Play" => MeansRotation(:Play, "Romeo and Juliet", "Hamlet", :Act, 1, 5, moderated=false),
    "Play" => MeansRotation(:Play, "Romeo and Juliet", "Hamlet", :Act, 1, 5),
    "LDA1" => LDARotation(:Play),
    "LDA1" => LDARotation(:Act),
    "MCMR1" => MulticlassRotation(:Play),
    "MCMR1" => MulticlassRotation(:Act)
]

for M in models
    for (label, rotation) in rotations
        myENA = M(
            data, codes, conversations, units,
            rotateBy=rotation
        )

        @testset "$(nameof(M)){$(nameof(typeof(rotation)))} runs" begin
            @test typeof(myENA) == M{typeof(rotation)}
        end

        @testset "$(nameof(M)){$(nameof(typeof(rotation)))} is correctly labeled" begin
            @test myENA.embedding[1, :label] == label
        end

        @testset "$(nameof(M)){$(nameof(typeof(rotation)))} has not too many dimensions" begin
            @test nrow(myENA.embedding) <= nrow(myENA.edges)
        end

        total_variance_explained = sum(myENA.embedding.variance_explained)
        @testset "$(nameof(M)){$(nameof(typeof(rotation)))} has correct variance: $(total_variance_explained) <= 1.0" begin
            @test isapprox(total_variance_explained, 1.0, atol=1e-8) || total_variance_explained <= 1.0
        end
        
        for i in 1:nrow(myENA.embedding)
            for j in (i+1):nrow(myENA.embedding)
                @testset "$(nameof(M)){$(nameof(typeof(rotation)))} axes $(i) and $(j) are orthogonal" begin
                    edgeIDs = myENA.edges.edgeID
                    @test isapprox(dot(myENA.embedding[i, edgeIDs], myENA.embedding[j, edgeIDs]), 0.0, atol=1e-8)
                end
            end
        end

        for i in 1:nrow(myENA.embedding)
            if startswith(myENA.embedding.label[i], "SVD")
                for j in (i+1):nrow(myENA.embedding)
                    if startswith(myENA.embedding.label[j], "SVD")
                        @testset "$(nameof(M)){$(nameof(typeof(rotation)))} $(myENA.embedding.label[i]) and $(myENA.embedding.label[j]) in proper order" begin
                            @test myENA.embedding.eigen_value[i] >= myENA.embedding.eigen_value[j]
                        end
                    end
                end
            end
        end

        @testset "$(nameof(M)){$(nameof(typeof(rotation)))} displays" begin
            display(myENA)
        end

        p = plot(myENA)
        @testset "$(nameof(M)){$(nameof(typeof(rotation)))} plots" begin
            @test typeof(p) <: Plots.Plot
        end

        @testset "$(nameof(M)){$(nameof(typeof(rotation)))} reconstructs" begin
            reENA = M(
                myENA;
                windowSize=4
            )

            @test typeof(reENA) == M{typeof(rotation)}
            @test reENA.embedding[1, :label] == label
            @test reENA.config.windowSize == 4
        end

        @testset "$(nameof(M)){$(nameof(typeof(rotation)))} trains" begin
            trainedENA = M(
                data, codes, conversations, units;
                rotateBy=TrainedRotation(myENA)
            )

            tp = plot(trainedENA)
            @test typeof(trainedENA) == M{TrainedRotation{typeof(myENA)}}
            @test isequal(trainedENA.nodes, myENA.nodes)
            @test trainedENA.embedding[1, :label] == label
            @test length(p.subplots) == length(tp.subplots)
        end

        @testset "$(nameof(M)){$(nameof(typeof(rotation)))} rerotates" begin
            trainedENA = M(
                myENA;
                rotateBy=TrainedRotation(myENA)
            )

            tp = plot(trainedENA)
            @test typeof(trainedENA) == M{TrainedRotation{typeof(myENA)}}
            @test trainedENA.embedding[1, :label] == label
            @test length(p.subplots) == length(tp.subplots)
        end

        @testset "$(nameof(M)){$(nameof(typeof(rotation)))} reconstructs and rerotates" begin
            trainedENA = M(
                myENA;
                windowSize=4,
                rotateBy=TrainedRotation(myENA)
            )

            tp = plot(trainedENA)
            @test typeof(trainedENA) == M{TrainedRotation{typeof(myENA)}}
            @test trainedENA.embedding[1, :label] == label
            @test trainedENA.config.windowSize == 4
            @test length(p.subplots) == length(tp.subplots)
        end
    end
end