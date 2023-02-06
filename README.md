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

## Installation

### Minimal Installation

To install this package locally, run the following in Julia:

```julia
using Pkg
Pkg.add(url="https://github.com/snotskie/EpistemicNetworkAnalysis.jl")
using EpistemicNetworkAnalysis
```

You may also want to install the following packages too:

```julia
Pkg.add([
    "CSV",
    "DataFrames",
    "Plots",
    "Statistics"
])
```

### Maximal Installation

Assuming you have Python installed, but nothing else, to install this package and configure everything to work in Jupyter, from the command line run:

```sh
JILL_INSTALL_DIR=$HOME/.opt/julias # edit if you want your julia binaries installed elsewhere
export JILL_INSTALL_DIR
python3 -m pip install jill jupyter
python3 -m jill install --confirm 1.8
PATH=$HOME/.local/bin:$PATH
export PATH
julia-1.8 -e 'using Pkg; Pkg.add(url="https://github.com/snotskie/EpistemicNetworkAnalysis.jl")'
julia-1.8 -e 'using Pkg; Pkg.add(["DataFrames", "Plots", "CSV", "Statistics", "IJulia"])'
```

Then launch Jupyter with:

```sh
python3 -m jupyter notebook
```

## Getting Started

To test your installation of this package, in Julia run:

```julia
using EpistemicNetworkAnalysis

data = ena_dataset("shakespeare.data")

conversations = [:Play, :Act]
units = [:Play, :Speaker]
codes = [
    :Love,
    :Death,
    :Honor,
    :Men,
    :Women
]

rotation = MeansRotation(:Play, "Romeo and Juliet", "Hamlet")

myENA = ENAModel(
    data, codes, conversations, units,
    rotateBy=rotation
)

p = plot(myENA)
display(p)
```

## Usage

### Loading Data

There are two built-in datasets for running examples:

```julia
data = ena_dataset("shakespeare.data")
data = ena_dataset("transitions.data")
```

To load your own data, do:

```julia
using CSV
using DataFrames

filename = "example.csv" # the name of the file to load
missingtext = "#N/A" # the text used to represent missing values in your CSV file
data = DataFrame(CSV.File(filename, missingstring=missingtext))
```

### ENA Model

TODO

### Rotations

TODO

### Digraph ENA Model

TODO

### Nonlinear ENA Model

TODO

### Plotting

TODO

