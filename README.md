# EpistemicNetworkAnalysis.jl

A port of [rENA](https://rdrr.io/cran/rENA/) version 0.2.0.1 into native Julia.

Original R package by [http://www.epistemicnetwork.org/](http://www.epistemicnetwork.org/).

IN DEVELOPMENT

## Math

### Measure

TODO

### Rotation

The x-axis values for each unit (x) are proportional to the sum of that unit's connections (c), weighted (w) to give a "rotation" of the high dimensional connection space that highlights some feature of interest.

![x_i \propto \sum_j w_j c_{ji}](https://render.githubusercontent.com/render/math?math=x_i%20%5Cpropto%20%5Csum_j%20w_j%20c_%7Bji%7D)

The y-axis values are then generally chosen to explain the most variance of the high dimensional space while being orthogonal to the x-axis (ie, orthogonal svd).

#### Regression Rotation

Given an explanatory variable of interest (u) and a set of other explanatory variables to control for (v), regress each of the connections (c) on those explanatory variables:

![c_j \sim 1 + u + v_1 + v_2 + ...](https://render.githubusercontent.com/render/math?math=c_j%20%5Csim%201%20%2B%20u%20%2B%20v_1%20%2B%20v_2%20%2B%20...)

The weight for a connection is then chosen to be proportional to the effect size of the variable of interest in that regression:

![w_j \propto \beta_u](https://render.githubusercontent.com/render/math?math=w_j%20%5Cpropto%20%5Cbeta_u)

This gives a rotation such that, holding v constant, there is a stronger positive effect between the connections on the right side of the x-axis, and and a stronger negative effect for those on the left.

#### Means Rotation

Use a regression rotation, where the explantory variable of interest is a dummy group variable (g), where 0 represents the control group and 1 the treatment group.

![c_j \sim 1 + g](https://render.githubusercontent.com/render/math?math=c_j%20%5Csim%201%20%2B%20g)

![w_j \propto \beta_g](https://render.githubusercontent.com/render/math?math=w_j%20%5Cpropto%20%5Cbeta_g)

This gives a rotation such that there is a stronger correlation between the connections on the right and being in the treatment group, and a stronger correlation between the connections on the left and being in the control group.

#### Moderated Means Rotation

Given a set of rival explanatory variables ("confounds", v), use a regression rotation, where the explantory variable of interest is a dummy group variable (g), where 0 represents the control group and 1 the treatment group; and also include in the regression the confounds and the interaction between the grouping variable and the confounds.

![c_j \sim 1 + g + v_1 + v_2 + ... + gv_1 + gv_2 + ...](https://render.githubusercontent.com/render/math?math=c_j%20%5Csim%201%20%2B%20g%20%2B%20v_1%20%2B%20v_2%20%2B%20...%20%2B%20gv_1%20%2B%20gv_2%20%2B%20...)

![w_j \propto \beta_g](https://render.githubusercontent.com/render/math?math=w_j%20%5Cpropto%20%5Cbeta_g)

This gives a rotation such that, holding the confounds and interaction with the confounds constant, there is a stronger correlation between the connections on the right and being in the treatment group, and a stronger correlation between the connections on the left and being in the control group.

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

https://alexanderrodin.com/github-latex-markdown/
