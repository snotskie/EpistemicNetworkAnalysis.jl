# EpistemicNetworkAnalysis.jl

A port of [rENA](https://rdrr.io/cran/rENA/) version 0.2.0.1 into native Julia.

Original R package by [http://www.epistemicnetwork.org/](http://www.epistemicnetwork.org/).

IN DEVELOPMENT

## Math

### Measure

TODO

### Rotation

![x_i \propto \sum_j w_j r_{ij}](https://render.githubusercontent.com/render/math?math=x_i%20%5Cpropto%20%5Csum_j%20w_j%20r_%7Bij%7D)

#### Regression Rotation

TODO

#### Means Rotation

TODO

#### Moderated Means Rotation

TODO

#### SVD Rotation

TODO

### Model

TODO

## Tests

(as of [969eb82](https://github.com/snotskie/EpistemicNetworkAnalysis.jl/commit/969eb822c7c8e420f0459c154a30ad2043062a42))

Means rotation:

![](examples/images/mr1.png)

Means rotation, group alignment plot:

![](examples/images/mr1-sub.png)

Means rotation, group alignment plot, detail showing target projections:

![](examples/images/mr1-sub-detail.png)

SVD rotation:

![](examples/images/svd.png)

SVD rotation, group alignment plot:

![](examples/images/svd-sub.png)

SVD rotation, group alignment plot, detail showing target projections:

![](examples/images/svd-sub-detail.png)

Regression rotation:

![](examples/images/rr1.png)

Regression rotation, group alignment plot:

![](examples/images/rr1-sub.png)

Regression rotation, group alignment plot, detail showing target projections:

![](examples/images/rr1-sub-detail.png)

Moderated means rotation:

![](examples/images/mmr1.png)

Moderated means rotation, group alignment plot:

![](examples/images/mmr1-sub.png)

Moderated means rotation, group alignment plot, detail showing target projections:

![](examples/images/mmr1-sub-detail.png)

AC-SVD rotation:

![](examples/images/acsvd.png)

AC-SVD rotation, group alignment plot:

![](examples/images/acsvd-sub.png)

AC-SVD rotation, group alignment plot, detail showing target projections:

![](examples/images/acsvd-sub-detail.png)
