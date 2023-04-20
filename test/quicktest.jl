# From the REPL, include this file to run quick tests, without having to be as thorough as runtests.jl
# and without having to deal with restarting the REPL between module loads

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

#=

TODO figure this out

In the main branch:

julia> myENA.accumModel[1, :]
DataFrameRow
 Row │ ENA_UNIT         pos_x     pos_y     Love_Death  Love_Honor  Death_Men  Men_Women  Love_Women  Love_Men  Death_Women  Honor_Men  Death_Honor  Honor_Women 
     │ String           Float64   Float64   Real        Real        Real       Real       Real        Real      Real         Real       Real         Real        
─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ Hamlet.BERNARDO  -0.20005  0.739359         0.0         0.0   0.948683        0.0         0.0       0.0          0.0   0.316228          0.0          0.0

julia> myENA.accumModel[2, :]
DataFrameRow
 Row │ ENA_UNIT          pos_x     pos_y      Love_Death  Love_Honor  Death_Men  Men_Women  Love_Women  Love_Men  Death_Women  Honor_Men  Death_Honor  Honor_Women 
     │ String            Float64   Float64    Real        Real        Real       Real       Real        Real      Real         Real       Real         Real        
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   2 │ Hamlet.FRANCISCO  0.481711  0.0758556           0           0          0          0           0         0            0          0            0            0

julia> myENA.accumModel[3, :]
DataFrameRow
 Row │ ENA_UNIT        pos_x      pos_y     Love_Death  Love_Honor  Death_Men  Men_Women  Love_Women  Love_Men  Death_Women  Honor_Men  Death_Honor  Honor_Women 
     │ String          Float64    Float64   Real        Real        Real       Real       Real        Real      Real         Real       Real         Real        
─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   3 │ Hamlet.HORATIO  -0.461365  0.215635   0.0505076         0.0   0.757614   0.454569         0.0  0.454569          0.0   0.101015          0.0          0.0

julia> myENA.accumModel[4, :]
DataFrameRow
 Row │ ENA_UNIT          pos_x      pos_y     Love_Death  Love_Honor  Death_Men  Men_Women  Love_Women  Love_Men  Death_Women  Honor_Men  Death_Honor  Honor_Women 
     │ String            Float64    Float64   Real        Real        Real       Real       Real        Real      Real         Real       Real         Real        
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   4 │ Hamlet.MARCELLUS  -0.331359  0.638816         0.0         0.0   0.948683        0.0         0.0  0.316228          0.0        0.0          0.0          0.0


In the refactor branch:

julia> myENA.accum[1, :]
DataFrameRow
 Row │ unitID           Love_Death  Love_Honor  Love_Men  Love_Women  Death_Honor  Death_Men  Death_Women  Honor_Men  Honor_Women  Men_Women 
     │ Symbol           Real        Real        Real      Real        Real         Real       Real         Real       Real         Real      
─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ Hamlet.BERNARDO         0.0         0.0       0.0         0.0          0.0        1.0          0.0        0.0          0.0        0.0

julia> myENA.accum[2, :]
DataFrameRow
 Row │ unitID            Love_Death  Love_Honor  Love_Men  Love_Women  Death_Honor  Death_Men  Death_Women  Honor_Men  Honor_Women  Men_Women 
     │ Symbol            Real        Real        Real      Real        Real         Real       Real         Real       Real         Real      
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   2 │ Hamlet.FRANCISCO           0           0         0           0            0          0            0          0            0          0

julia> myENA.accum[3, :]
DataFrameRow
 Row │ unitID          Love_Death  Love_Honor  Love_Men  Love_Women  Death_Honor  Death_Men  Death_Women  Honor_Men  Honor_Women  Men_Women 
     │ Symbol          Real        Real        Real      Real        Real         Real       Real         Real       Real         Real      
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   3 │ Hamlet.HORATIO         0.0         0.0       0.0         0.0          0.0   0.948683          0.0   0.316228          0.0        0.0

julia> myENA.accum[4, :]
DataFrameRow
 Row │ unitID            Love_Death  Love_Honor  Love_Men  Love_Women  Death_Honor  Death_Men  Death_Women  Honor_Men  Honor_Women  Men_Women 
     │ Symbol            Real        Real        Real      Real        Real         Real       Real         Real       Real         Real      
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   4 │ Hamlet.MARCELLUS           0           0         0           0            0          0            0          0            0          0


The only difference I can find is that the refactor uses Inf as the default windowSize, but that shouldn't cause this difference.

Could it be accumulating to the wrong person somehow? Mixing up values?

The plots are fairly similar, which is consistent with row shuffling

=#