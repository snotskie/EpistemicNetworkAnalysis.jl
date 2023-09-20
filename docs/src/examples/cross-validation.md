# Cross Validation

Using a [TrainedRotation](../rotations.md#TrainedRotation), we can perform [k-fold cross-validation](https://neptune.ai/blog/cross-validation-in-machine-learning-how-to-do-it-right) on an ENA model

First, we'll load a few libraries:

```@example crossValidation
using EpistemicNetworkAnalysis
using DataFrames
using Statistics
using GLM
nothing # hide
```

Second, we'll load our data and prepare our model config. We'll be using the `FormulaRotation` example from the [ICQE23](../icqe23.md) workshop:

```@example crossValidation
data = loadExample("transitions")

deriveAnyCode!(data, :BODY, :Changes, :Mood, :Oily, :Dysphoria, :Cry)
deriveAnyCode!(data, :REFLECT, :Identity, :Longing, :Dream, :Childhood, :Family, :Name, :Letter, :Doubt, :Religion)
deriveAnyCode!(data, :LEARN, :WWW, :Experiment, :Recipe)
deriveAnyCode!(data, :PROGRESS, :Strangers, :Passed, :Out, :Affirmation)

data[!, :All] .= "All"
codes = [:DoseTracking, :SkippedDose, :Happy, :NonHappy, :Sweets, :BODY, :REFLECT, :LEARN, :PROGRESS]
conversations = [:All]
units = [:Date]
rotation = FormulaRotation(
    LinearModel, @formula(y ~ 1 + Day), 2, nothing
)
nothing # hide
```

Now we can start setting up our cross-validation. We'll give each row a random number from 1 to 5, setting us up for a 5-fold cross-validation.

```@example crossValidation
k_folds = 5
data[!, :Fold] .= rand(1:k_folds, nrow(data))
nothing # hide
```

Then, we'll iterate. We'll create a `trainmodel` with a `unitFilter`, using the logic `row.Fold != i` to select all units except our hold out set. After that, we'll create a `testmodel` with the opposite `unitFilter` and rotate it using `TrainedRotation(trainmodel)`. That will project our hold out units into our trained embedding. The last thing we'll do in this loop is grab a statistic to add to a `results` list:

```@example crossValidation
results = Real[]
for i in 1:k_folds
    trainmodel = ENAModel(
        data, codes, conversations, units,
        windowSize=4,
        recenterEmpty=true,
        rotateBy=rotation,
        unitFilter=(row)->(row.Fold != i)
    )

    testmodel = ENAModel(
        data, codes, conversations, units,
        windowSize=4,
        recenterEmpty=true,
        rotateBy=TrainedRotation(trainmodel),
        unitFilter=(row)->(row.Fold == i)
    )

    result = testmodel.embedding[1, :Formula_AdjR2]
    push!(results, result)
end
nothing # hide
```

Finally, we'll display the results and their mean:

```@example crossValidation
println(results)
println(mean(results))
```

Putting it all together, here is a helper function you should be able to drop-in and apply to your own data:

```@example crossValidation
# Helper
function kfoldcv(wholemodel, k_folds, statistic)
    results = Real[]
    wholemodel.data[!, :Fold] .= rand(1:k_folds, nrow(data))
    for i in 1:k_folds
        trainmodel = ENAModel(
            wholemodel,
            unitFilter=(row)->(row.Fold != i)
        )

        testmodel = ENAModel(
            wholemodel,
            rotateBy=TrainedRotation(trainmodel),
            unitFilter=(row)->(row.Fold == i)
        )

        result = testmodel.embedding[1, statistic]
        push!(results, result)
    end
    
    return results
end

# Example usage
wholemodel = ENAModel(
    data, codes, conversations, units,
    windowSize=4,
    recenterEmpty=true,
    rotateBy=rotation
)

results = kfoldcv(wholemodel, 5, :Formula_AdjR2)
println(results)
println(mean(results))
```