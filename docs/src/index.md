# EpistemicNetworkAnalysis.jl

**Author: Mariah A. Knowles (@snotskie)**

A port of [rENA](https://rdrr.io/cran/rENA/) version 0.2.0.1 into native Julia, with substantial API changes to fit Julia style, and with addition of new rotation methods.

Original R package by [http://www.epistemicnetwork.org/](http://www.epistemicnetwork.org/).

## Recommended Citation

If you use this Julia implementation of ENA in your research, please cite it, the Julia language, and the R version it is based on, as shown below; describe it as "...a port of ENA from the original R to Julia (Bezanson et al, 2017; Knowles, 2021; Marquart et al 2019)"; and justify your choice of Julia over R. For example, your research team may be more familiar with Julia, you may need interoperability with other existing Julia libraries, or you may need access to dimension reductions ("rotations") not defined in rENA.

> Bezanson, J., Edelman, A., Karpinski S., & Shah V. B. (2017). Julia: A Fresh Approach to Numerical Computing. SIAM Review, 59: 65-98.
> 
> Knowles, M. A. (2023). EpistemicNetworkAnalysis.jl (Version 0.5.0) \[Software\] Available from https://github.com/snotskie/EpistemicNetworkAnalysis.jl
> 
> Marquart, C. L., Swiecki, Z., Collier, W., Eagan, B., Woodward, R., & Shaffer, D. W. (2019). rENA: Epistemic Network Analysis (Version 0.2.0.1) \[Software\] Available from https://cran.r-project.org/web/packages/rENA/index.html

## Getting Started

I recommend installing [Julia for Visual Code](https://code.visualstudio.com/docs/languages/julia) first.

Then to install this package, run the following in Julia:

```julia
using Pkg
Pkg.add(url="https://github.com/snotskie/EpistemicNetworkAnalysis.jl")
using EpistemicNetworkAnalysis
```

Once you've done that, test that everything works by running the following:

```julia
using EpistemicNetworkAnalysis

data = loadExample("shakespeare")
codes = [:Love, :Death, :Honor, :Men, :Women]
conversations = [:Play, :Act, :Scene]
units = [:Play, :Speaker]
rotation = MeansRotation(:Play, "Romeo and Juliet", "Hamlet")

model = ENAModel(
    data, codes, conversations, units,
    windowSize=4,
    rotateBy=rotation
)

show(model)
show(plot(model))
```

If you run into any issues, don't hesitate to [file an issue or ask for help](https://github.com/snotskie/EpistemicNetworkAnalysis.jl/issues).

Then once you're ready, continue to learn more about [available models](models.md).