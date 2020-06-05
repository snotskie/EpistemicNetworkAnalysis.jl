using EpistemicNetworkAnalysis

RSdata = ena_dataset("RS.data")
display(first(RSdata, 6))

codes = [:Data,
    :Technical_Constraints,
    :Performance_Parameters,
    :Client_and_Consultant_Requests,
    :Design_Reasoning,
    :Collaboration]

conversations = [:Condition, :GameHalf, :GroupName]
units = [:Condition, :GameHalf, :UserName]
myRotation = MeansRotation(:Condition, "FirstGame", "SecondGame")
myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=myRotation)
display(myENA)

myArtist = MeansArtist(:GameHalf, "First", "Second")
p = plot(myENA,
    showprojection=true,
    artist=myArtist
)

display(p)