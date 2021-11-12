# EpistemicNetworkAnalysis.jl

A port of [rENA](https://rdrr.io/cran/rENA/) version 0.2.0.1 into native Julia, with substantial API changes to fit Julia style, and with addition of new rotation methods.

Original R package by [http://www.epistemicnetwork.org/](http://www.epistemicnetwork.org/).

ALPHA

## Recommended Citation

If you use this Julia implementation of ENA in your research, please cite it, the Julia language, and the R version it is based on, as shown below; describe it as "...an experimental port of ENA from the original R to Julia (Bezanson et al, 2017; Knowles, 2021; Marquart et al 2019)"; and justify your choice of Julia over R. For example, your research team may be more familiar with Julia or you may need interoperability with other existing Julia libraries.

> Bezanson, J., Edelman, A., Karpinski S., & Shah V. B. (2017). Julia: A Fresh Approach to Numerical Computing. SIAM Review, 59: 65-98.
> 
> Knowles, M. A. (2021). EpistemicNetworkAnalysis.jl (Version 0.1.0) \[Software\] Available from https://github.com/snotskie/EpistemicNetworkAnalysis.jl
> 
> Marquart, C. L., Swiecki, Z., Collier, W., Eagan, B., Woodward, R., & Shaffer, D. W. (2019). rENA: Epistemic Network Analysis (Version 0.2.0.1) \[Software\] Available from https://cran.r-project.org/web/packages/rENA/index.html

## Getting Started

Dive in with examples by launching this project on binder! Click here: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/snotskie/EpistemicNetworkAnalysis.jl/HEAD)

## Documentation

See latest docs here: [https://snotskie.github.io/EpistemicNetworkAnalysis.jl/latest/](https://snotskie.github.io/EpistemicNetworkAnalysis.jl/latest/)

## Installation

To install this package locally, run the following in Julia:

```julia
using Pkg
Pkg.add(url="https://github.com/snotskie/EpistemicNetworkAnalysis.jl")
using EpistemicNetworkAnalysis
```
