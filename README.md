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

```julia
function mySubsetFilter(unit)
    if 1 <= unit[:Act] <= 3 # select only acts between 1 and 3
        return true
    else
        return false
    end
end

function myRelationshipFilter(i, j, ci, cj)
    if ci == :Love && cj == :Death
        return false # don't connect love and death
    elseif cj == :Honor && ci == :Women
        return false # don't connect honor and women
    else
        return i < j # connect everything else, just one direction for each pair of codes
    end
end

myENA = ENAModel(
    data, codes, conversations, units,
    windowSize=4,
    sphereNormalize=true,
    dropEmpty=false,
    meanCenter=true,
    deflateEmpty=false,
    subspaces=0,
    fitNodesToCircle=false,
    subsetFilter=mySubsetFilter,
    relationshipFilter=myRelationshipFilter,
    rotateBy=rotation
)
```

#### Required Parameters

The `ENAModel` constructor has four required parameters:

1. A DataFrame containing the data to be analyzed
2. An array of symbols specifying the column names of the qualitative codes to include in the model. The order of the array makes no difference.
3. An array of symbols used to divide the rows into units of analysis. The order makes no difference.
4. Similarly, an array of symbols used to divide the rows into conversations. The order makes no difference.

#### Optional Parameters

The `ENAModel` constructor has multiple optional named parameters that are primitive types:

- `windowSize`, an integer used to specify the size of the sliding stanza window. By default, this is 4. To specify an infinite stanza window, use a very large integer value.
- `sphereNormalize`, a boolean that tells the model whether to normalize units to the sphere. By default, this is true.
- `dropEmpty`, a boolean that tells the model whether to drop units with empty networks. By default, this is false.
- `meanCenter`, a boolean that tells the model whether to center the units such that the mean lies at the origin. By default, this is true. When `dropEmpty=true` and `meanCenter=false`, the zero-network will lie at the origin instead.
- `recenterEmpty`, a boolean that tells the model whether to recenter units with empty networks at the mean. By default, this is false.
- `deflateEmpty`, a boolean that tells the model whether to deflate the "umbrella handle" from the high dimensional space used for rotation. By default, this is `false`. When it is true, the variance that runs from the zero-network to the mean-network is removed. When `deflateEmpty=true`, changing the value of `meanCenter` will have no effect.
- `subspaces`, an integer `N >= 2` that tells the model whether to project the high dimensional space used for rotation to the subspace formed by the first `N` SVD dimensions. By default, this is `0`, meaning to use the entire high dimensional space for rotation.
- `fitNodesToCircle`, a boolean that tells the model whether to degrade the optimized positions of the code positions by forcing them to the unit circle instead, which sacrifices accuracy of the model in favor of readability

#### Subset Filtering

The `ENAModel` constructor takes another optional named parameter, `subsetFilter`.

This parameter is a lambda function that is used by the model to remove units from the model after accumulation.

By default, a lambda is used that keeps all units.

#### Relationship Filtering

The `ENAModel` constructor takes another optional named parameter, `relationshipFilter`.

This parameter is a lambda funciton that is used by the model to decide which pairs of codes to include as relationships in the model.

By default, a lambda is used that keeps all pairs of codes where the first code occurs before the second code in the `codes` parameter.

Because `ENAModel` is undirected, one should avoid including both orderings of any pair of codes in the list of relationships, for example `X_Y` and `Y_X`.

#### Rotation Parameter

Finally, the `ENAModel` constructor has one more optional named parameter, `rotateBy`.

See the Rotations section for more detail.

### Rotations

TODO

### Plotting

TODO

### Exporting Data

To export results, you can save the `accumModel`, `codeModel`, and `networkModel` fields of your ENA model object to a CSV.

```julia
using CSV
CSV.write(myENA.accumModel)
CSV.write(myENA.codeModel)
CSV.write(myENA.networkModel)
```

You can also save your plot results using `savefig`, which supports multiple file types.

```julia
using Plots
savefig(p, "my-output.png")
savefig(p, "my-output.jpeg")
savefig(p, "my-output.svg")
```

### Digraph ENA Model

By default, `ENAModel` is configured for undirected data: the order of which code response to which other code does not matter.

For seeing the directional effects, consider using a `DigraphENAModel` instead. This alternative has the same configuration options as `ENAModel`. The difference is how it embeds and visualizes the network so as to get a quick overview of how directed connections between codes "drive" the space. 

```julia
myENA = DigraphENAModel(data, codes, conversations, units)
```

### Nonlinear ENA Model

TODO

