# Rotations

The following documentation describes available dimension reduction, or "rotations," available in this package.

Which rotation you choose should be informed by your sense of your research story:

- Want to get a high level view of your data? `SVDRotation` will capture the most variance and is good for getting an initial sense of your data's major features
- Want to compare groups? Use `MeanRotation` when you have exactly two groups. Otherwise, consider `LDARotation` or `MulticlassRotation`. The distinction between these is spelled out in ["Multiclass Rotations in Epistemic Network Analysis"](https://link.springer.com/chapter/10.1007/978-3-031-31726-2_5)
- Want to see how your qualitative data relates to continuous or hierarchical variables? Use `FormulaRotation`, a rotation based on a regression framework described in ["Hierarchical Epistemic Network Analysis"](https://www.qesoc.org/images/pdf/ICQE20_Proceedings_Supplement_Final_web.pdf#page=37)
- Want to model a theoretical topic directly? Use `TopicRotation` to force certain connections to the left and others to the right. This is useful when you have suspicions that certain codes account for group differences and want to test those suspicions directly instead of inferring them by other means.
- Want to test the generlizability or cross-validation of your model? Use `TrainedRotation` to create a new model that uses the same embedding and runs the same tests as an existing model.

Once you're familiar with them, continue to learn more about [plotting](plots.md).

## SVDRotation

```@docs
EpistemicNetworkAnalysis.SVDRotation
```

## MeansRotation

```@docs
EpistemicNetworkAnalysis.MeansRotation
```

## LDARotation

```@docs
EpistemicNetworkAnalysis.LDARotation
```

## MulticlassRotation

```@docs
EpistemicNetworkAnalysis.MulticlassRotation
```

## FormulaRotation

```@docs
EpistemicNetworkAnalysis.FormulaRotation
```

## TopicRotation

```@docs
EpistemicNetworkAnalysis.TopicRotation
```

## TrainedRotation

```@docs
EpistemicNetworkAnalysis.TrainedRotation
```