# From the REPL, include this file to run quick tests, without having to be as thorough as runtests.jl
# and without having to deal with restarting the REPL between module loads

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
#     :Act, 1, 5,
#     moderated=true
# )

rotation = EpistemicNetworkAnalysis.TopicRotation(
    "My Topic",
    [:Women],
    [:Men]
)

model = EpistemicNetworkAnalysis.ENAModel(
    data, codes, conversations, units,
    windowSize=4,
    # rotateBy=rotation,
    recenterEmpty=true,
    # dropEmpty=true,
)

# model2 = EpistemicNetworkAnalysis.BiplotENAModel(model)
model2 = EpistemicNetworkAnalysis.BiplotENAModel(model, rotateBy=rotation)
# model2 = EpistemicNetworkAnalysis.ENAModel(model, rotateBy=rotation)

p = EpistemicNetworkAnalysis.plot(model2, showWeakEdges=false) #, groupBy=group)
# p = EpistemicNetworkAnalysis.plot(model, groupBy=group, lims=2)
# p = EpistemicNetworkAnalysis.plot(model, groupBy=group, x=3, y=4)

#=
TODO:

- [X] base models
- [X] base rotations
- [X] base plotting
- [ ] trajectories and spectral
- [ ] statistical tests
- [ ] textual summaries
- [ ] digraph model and plotting
- [ ] LDA and Multiclass
- [ ] default exports
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