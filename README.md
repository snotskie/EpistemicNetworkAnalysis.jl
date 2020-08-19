# EpistemicNetworkAnalysis.jl

A port of [rENA](https://rdrr.io/cran/rENA/) version 0.2.0.1 into native Julia, with substantial API changes to fit Julia style, and with addition of new rotation methods.

Original R package by [http://www.epistemicnetwork.org/](http://www.epistemicnetwork.org/).

ALPHA

## Installation

```julia
using Pkg
Pkg.install("https://github.com/snotskie/EpistemicNetworkAnalysis.jl")
```

## Usage

```julia
using Colors
using EpistemicNetworkAnalysis

# 1. Load your data into a DataFrame
data = ena_dataset("RS.data")
display(data)

# 2. Specify your model
## 2a. What columns categorize each row? (numeric columns, usually containing only 0 or 1)
codes = [
    :Data,
    :Technical_Constraints,
    :Performance_Parameters,
    :Client_and_Consultant_Requests,
    :Design_Reasoning,
    :Collaboration
]

## 2b. How do we specify unique conversations in the df?
conversations = [:Condition, :GameHalf, :GroupName]

## 2c. How do we specify unique speakers in the df?
units = [:Condition, :GameHalf, :UserName]

# 3. Specify the rotation method (optional, defaults to SVD)
rotation = MeansRotation(:Condition, "FirstGame", "SecondGame")

# 4. Run the model
ena = ENAModel(data, codes, conversations, units, rotateBy=rotation)
display(ena)

# 5. Plot the results
p = plot(ena, title="Example", ylabel="SVD", xlabel="Condition", minColor=colorant"blue", maxColor=colorant"red")
display(p)
```