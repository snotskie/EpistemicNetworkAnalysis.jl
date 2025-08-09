include("../src/EpistemicNetworkAnalysis.jl")
using Plots

data = EpistemicNetworkAnalysis.loadExample("efm")
codes = Symbol.([
    "Right Image",
    "Expect Ethics",
    "Ethics is Blocked",
    "Identity Formation",
    "Belonging",
    "Sociotechnical Change",
    "Advancements",
    "Act with Integrity",
    "Reason",
    "Social and Collaborative",
    "Issues",
    "Having a Job",
    "Put into Practice",
    "Hype and Purpose",
    "System Stakeholders",
    "Hands-on and Exposure",
    "Technical Knowledge",
    "Luck",
    "Add Value",
    "Law and Policy",
    "Driven by Capitalism",
    "Virtue",
    "Clarity",
    "Reflexive",
    "Internal Motivation",
    "Scope of Work",
    "Teach and Mentor",
    "Shared Language",
])
convos = [:MemoID]
units = [:Kind, :CodewiseUnit]
rotation = EpistemicNetworkAnalysis.MulticlassRotation(:Kind)
model = EpistemicNetworkAnalysis.CodewiseENAModel(data, codes, convos, units, windowSize=4, rotateBy=rotation)

plot(model, fitNodesToCircle=true, zoom=0.7)