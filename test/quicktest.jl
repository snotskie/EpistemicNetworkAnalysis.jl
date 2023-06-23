# From the REPL, include this file to run quick tests, without having to be as thorough as runtests.jl
# and without having to deal with restarting the REPL between module loads

using DataFrames
using GLM
using Serialization
include("../src/EpistemicNetworkAnalysis.jl")

# data = EpistemicNetworkAnalysis.loadExample("toy")
# conversations = [:Convo]
# units = [:Unit]
# codes = [:A, :B, :C]
# group = :Group
# rotation = EpistemicNetworkAnalysis.MeansRotation(
#     :Group, "Red", "Blue"
# )

data = EpistemicNetworkAnalysis.loadExample("shakespeare")
data[!, :rand] = rand(1:20, nrow(data))
conversations = [:Play, :Act, :Scene]
units = [:Play, :Speaker]
codes = [
   :Love,
   :Death,
   :Honor,
   :Men,
   :Women
]
group = :Play
# rotation = EpistemicNetworkAnalysis.MeansRotation(
#     :Play, "Romeo and Juliet", "Hamlet",
#     # :Act, 1, 5,
#     # moderated=true
# )

rotation = EpistemicNetworkAnalysis.MulticlassRotation(
    :Play
    # :Act
)

# rotation = EpistemicNetworkAnalysis.TopicRotation(
#     "My Topic",
#     [:Women],
#     [:Men]
# )

# rotation = EpistemicNetworkAnalysis.FormulaRotation(
#     LinearModel, @formula(y ~ 1 + rand), 2, nothing
# )

model = EpistemicNetworkAnalysis.ENAModel(
    data, codes, conversations, units,
    windowSize=4,
    rotateBy=rotation,
    # recenterEmpty=true,
    dropEmpty=true,
)

# model2 = EpistemicNetworkAnalysis.DigraphENAModel(
#     data, codes, conversations, units,
#     windowSize=10,
#     # rotateBy=EpistemicNetworkAnalysis.ManualRotation(model.embedding),
#     # recenterEmpty=true,
#     rotateBy=rotation,
#     dropEmpty=true,
# )

@show(EpistemicNetworkAnalysis.summary(model))
@show(model)
EpistemicNetworkAnalysis.to_xlsx(model, "test/temp.xlsx")
serialize("test/temp.ena", model)
modeldes = deserialize("test/temp.ena")
@show(modeldes)

# modelxlsx = EpistemicNetworkAnalysis.from_xlsx("test/temp.xlsx")
# @show(modelxlsx)

# model2 = EpistemicNetworkAnalysis.BiplotENAModel(model)
# model2 = EpistemicNetworkAnalysis.BiplotENAModel(model, rotateBy=rotation)
# model2 = EpistemicNetworkAnalysis.ENAModel(model, rotateBy=rotation)

p = EpistemicNetworkAnalysis.plot(
    model,
    showWeakEdges=false,
    zoom=.75,
    # trajectoryBy=:Act,
    # trajectoryBy=:rand,
    # groupBy=group,
    # x=3,
    # y=4,
)

display(p)

#=
TODO:

- [X] base models
- [X] base rotations
- [X] base plotting
- [X] trajectories and spectral
- [X] statistical tests
- [X] textual summaries
- [X] digraph model and plotting
- [X] xlsx export
- [~] xlsx import - for filters, do something like x -> x.id in model.whatever.ids ?
- [X] LDA and Multiclass rotations, plots, and stats tests
- [ ] ManualRotation/CopyRotation statistical tests
- [ ] module exports
- [ ] plot tweaks (eg. arrow sizes, cutoff label bugfix when many groups, etc.)
- [ ] volunteer testing
- [ ] auto-docs
- [ ] cleaner errors and warnings (see volunteer testing results)
- [ ] front-facing documentation and citations
- [ ] Pluto integration/demos
- [ ] developer demos
- [ ] test compatibility with other plot backends
- [ ] test interaction with ML ecosystem
- [ ] clean up old files
- [ ] clean up comments
- [ ] pin version
- [ ] julia package
- [ ] web presence and qehub
- [ ] code-wise
- [ ] nonlinear
=#