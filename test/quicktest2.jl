# From the REPL, include this file to run quick tests, without having to be as thorough as runtests.jl
# and without having to deal with restarting the REPL between module loads

include("../src/EpistemicNetworkAnalysis.jl")

# Load sample dataset, codes from my first year on hormone replacement therapy
data = EpistemicNetworkAnalysis.loadExample("transitions") # NOTE: To load your own data, see DataFrame(CSV.File(...))

# Derive some new codes based on old ones
EpistemicNetworkAnalysis.deriveAnyCode!(data, :BODY, :Changes, :Mood, :Oily, :Dysphoria, :Cry)
EpistemicNetworkAnalysis.deriveAnyCode!(data, :REFLECT, :Identity, :Longing, :Dream, :Childhood, :Family, :Name, :Letter, :Doubt, :Religion)
EpistemicNetworkAnalysis.deriveAnyCode!(data, :LEARN, :WWW, :Experiment, :Recipe)
EpistemicNetworkAnalysis.deriveAnyCode!(data, :PROGRESS, :Strangers, :Passed, :Out, :Affirmation)

# Add new columns for splitting the year's data in half, third, ...
# data[!, :Half] .= "First"
# data[183:end, :Half] .= "Second"
# data[!, :Third] .= "First"
# data[122:243, :Third] .= "Second"
# data[244:end, :Third] .= "Third"
# data[!, :Fourth] .= "First"
# data[92:183, :Fourth] .= "Second"
# data[184:275, :Fourth] .= "Third"
# data[276:end, :Fourth] .= "Fourth"

# List columns to use as codes, convos, and units
codes = [:DoseTracking, :SkippedDose, :Happy, :NonHappy, :Sweets, :BODY, :REFLECT, :LEARN, :PROGRESS]
conversations = []
units = [:Date]

# Run the model and plot it
model = EpistemicNetworkAnalysis.ENAModel(data, codes, conversations, units)
p = EpistemicNetworkAnalysis.plot(model)
display(p)