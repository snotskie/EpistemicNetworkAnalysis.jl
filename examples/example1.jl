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
myRotation = MeansRotation(groupVar, controlGroup, treatmentGroup)
myENA = ENAModel(RSdata, codes, conversations, units, rotateBy=myRotation)
display(myENA)

myArtist = MeansArtist(:GameHalf, "First", "Second")
scene = plot(myENA,
    showprojection=true,
    unitscale=0,
    artist=myArtist
)
display(scene)