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
rotation = EpistemicNetworkAnalysis.MeansRotation(
    :Play, "Romeo and Juliet", "Hamlet",
    :Act, 1, 5,
)

myENA = EpistemicNetworkAnalysis.ENAModel(
    data, codes, conversations, units,
    windowSize=4,
    rotateBy=rotation,
    unitFilter=(row -> row.Act in [1, 5])
)

p = EpistemicNetworkAnalysis.plot(myENA, weakLinks=false, flipY=true)
# p = EpistemicNetworkAnalysis.plot(myENA, groupBy=group, lims=2)
# p = EpistemicNetworkAnalysis.plot(myENA, groupBy=group, x=3, y=4)