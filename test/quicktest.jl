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
# group = :Play
# rotation = EpistemicNetworkAnalysis.MeansRotation(
#     :Play, "Romeo and Juliet", "Hamlet",
#     # :Act, 1, 5,
#     # moderated=true
# )

# rotation = EpistemicNetworkAnalysis.LDARotation(
#     :Play
#     # :Act
# )

rotation = EpistemicNetworkAnalysis.MulticlassRotation(
    # :Play
    :Act
)

# rotation = EpistemicNetworkAnalysis.TopicRotation(
#     "Women-Death vs. Honor",
#     [:Women, :Death],
#     [:Honor]
# )

# rotation = EpistemicNetworkAnalysis.FormulaRotation(
#     LinearModel, @formula(y ~ 1 + rand), 2, nothing
# )

model = EpistemicNetworkAnalysis.ENAModel(
    data, codes, conversations, units,
    windowSize=4,
    rotateBy=rotation,
    # recenterEmpty=true,
    # unitFilter=row->row.Act in [1, 2, 5],
    dropEmpty=true,
)

model2 = EpistemicNetworkAnalysis.ENAModel(
    data, codes, conversations, units,
    windowSize=10,
    rotateBy=EpistemicNetworkAnalysis.TrainedRotation(model),
    # recenterEmpty=true,
    # rotateBy=rotation,
    # unitFilter=row->row.Act in [1, 2, 5],
    dropEmpty=true,
)

@show(EpistemicNetworkAnalysis.summary(model2))
@show(model2)
EpistemicNetworkAnalysis.to_xlsx("test/temp.xlsx", model2)
serialize("test/temp.ena", model)
modeldes = deserialize("test/temp.ena")
@show(modeldes)

# modelxlsx = EpistemicNetworkAnalysis.from_xlsx("test/temp.xlsx")
# @show(modelxlsx)

# model2 = EpistemicNetworkAnalysis.BiplotENAModel(model)
# model2 = EpistemicNetworkAnalysis.BiplotENAModel(model, rotateBy=rotation)
# model2 = EpistemicNetworkAnalysis.ENAModel(model, rotateBy=rotation)

p = EpistemicNetworkAnalysis.plot(
    model2,
    # showWeakEdges=false,
    # zoom=.6,
    # trajectoryBy=:Act,
    # trajectoryBy=:rand,
    # groupBy=group,
    # x=3,
    # y=4,
)

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
- [X] CopyRotation statistical tests
- [~] plot tweaks (eg. arrow sizes, cutoff label bugfix when many groups, etc.)
- [X] module exports
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