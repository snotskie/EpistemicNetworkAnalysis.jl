include("../src/EpistemicNetworkAnalysis.jl")
data = EpistemicNetworkAnalysis.loadExample("shakespeare.data")
conversations = [:Play, :Act]
units = [:Play, :Speaker]
codes = [
    :Love,
    :Death,
    :Honor,
    :Men,
    :Women
]

myENA = EpistemicNetworkAnalysis.ENAModel(data, codes, conversations, units)
p = EpistemicNetworkAnalysis.plot(myENA, groupBy=:Play)