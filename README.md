# EpistemicNetworkAnalysis.jl

[!["Buy Me A Coffee"](https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png)](https://www.buymeacoffee.com/snotskie)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8401215.svg)](https://doi.org/10.5281/zenodo.8401215)

**Author: Mariah A. Knowles (@snotskie)**

A port of [rENA](https://rdrr.io/cran/rENA/) version 0.2.0.1 into native Julia, with substantial API changes to fit Julia style, and with addition of new rotation methods.

Original R package by [http://www.epistemicnetwork.org/](http://www.epistemicnetwork.org/).

**This package is in BETA -- It is stable enough for exploratory work and ready for user feedback**

## Recommended Citation

If you use this Julia implementation of ENA in your research, please cite it, the Julia language, and the R version it is based on, as shown below; describe it as "...a port of ENA from the original R to Julia (Bezanson et al, 2017; Knowles, 2023; Marquart et al 2019)"; and justify your choice of Julia over R. For example, your research team may be more familiar with Julia, you may need interoperability with other existing Julia libraries, or you may need access to dimension reductions ("rotations") not defined in rENA.

> Bezanson, J., Edelman, A., Karpinski S., & Shah V. B. (2017). Julia: A Fresh Approach to Numerical Computing. SIAM Review, 59: 65-98.
> 
> Knowles, M. A. (2023). EpistemicNetworkAnalysis.jl (Version 0.7.0) \[Software\] Available from https://github.com/snotskie/EpistemicNetworkAnalysis.jl
> 
> Marquart, C. L., Swiecki, Z., Collier, W., Eagan, B., Woodward, R., & Shaffer, D. W. (2019). rENA: Epistemic Network Analysis (Version 0.2.0.1) \[Software\] Available from https://cran.r-project.org/web/packages/rENA/index.html

You may copy the recommended citation in the format of your choice from <https://zenodo.org/doi/10.5281/zenodo.8401214>

## Getting Started

I recommend installing [Julia for Visual Code](https://code.visualstudio.com/docs/languages/julia) first. If you are in a country where downloading resources from GitHub is difficult, you may install Julia using the [Python jill package](https://pypi.org/project/jill/) and following the advice linked in the [Julia forums for users in China](https://discourse.julialang.org/t/why-is-it-so-hard-to-install-packages-in-julia-from-china/87566)

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
display(plot(model))
```

If you run into any issues, don't hesitate to [file an issue or ask for help](https://github.com/snotskie/EpistemicNetworkAnalysis.jl/issues).

Then once you're ready, continue to learn more about [available models](https://snotskie.github.io/EpistemicNetworkAnalysis.jl/latest/models/).