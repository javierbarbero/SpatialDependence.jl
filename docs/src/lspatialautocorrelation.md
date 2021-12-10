```@meta
CurrentModule = SpatialDependence
```

# Local Spatial Autocorrelation

Local spatial autocorrelation statistics are used to discover clusters (hot spots or cold spots) and spatial outliers. Local Indicators of Spatial Association (LISA), introduced by Anselin (1995), are local spatial statistics that provide a statistic for each location and assess its contribution to the corresponding global spatial autocorrelation statistic.

Inference is performed through a computational conditional permutation approach by randomly reshuffling the values of the variable to different locations while the value of the observation under consideration is held fixed at its location. The observed statistic for each location is compared to a reference distribution under the null hypothesis of spatial randomness.

All the examples on this page assume that the Guerry dataset is loaded and a polygon contiguity spatial weights object has been built and row-standardized.
```@example lscor
using SpatialDatasets
using SpatialDependence
using StableRNGs
using Plots

guerry = sdataset("Guerry")
W = polyneigh(guerry) 
nothing # hide
```

## Local Moran

Local Moran (Anselin, 1995) is the most used local spatial autocorrelation statistic. It is computed as:
```math
I = \frac{z_i}{m_2} \sum_{j}w_{ij}z_{j}
```
where $z$ is the variable of interest in deviations from the mean, and $m_2 = \sum_{i}z_i / (n - 1)$ or  $m_2 = \sum_{i}z_i / n$ is the scaling factor.

Local Moran I can be computed with the `localmoran` function. By default, $9,999$ permutations are calculated for the inference. It is possible to specify a different number of permutations with the `permutations` optional parameter. For reproduciibility, it is possible to specify a custom random number generator with the `rng` optional parameter. If `corrected` is set to `false` the scaling factor is divided by $n$ instead of $n - 1$.
```@example lscor
lmguerry = localmoran(guerry.Litercy, W, permutations = 9999, rng = StableRNG(1234567))
```

In the Local Moran, observations are classified in four categories, depending on the value of the attribute and the value of their neighbors with respect to the mean:

| Code   | Category    | Cluster or Outliers | Interpretation                         |
|:-------|:------------|:--------------------|:---------------------------------------|
| `:HH`  | High-High   | Cluster: hot spot   | High values surrounded by high values  |
| `:LL`  | Low-Low     | Cluster: cold spot  | Low values surrounded by low values    |
| `:LH`  | Low-High    | Outlier: doughnut   | High values surrounded by low values   |
| `:HL`  | High-Low    | Outlier: diamond    | Low values surrounded by high values   |

The category for whicheach observation is assigned can be retrieved with the `assignments` function:
```@example lscor
assignments(lmguerry)
```

## Significance

The conditional randomization procedure returns a pseudo p-value that can be used to assed the significance of the identified clusters and spatial outliers. The function `issignificant` returns a vector of boolean values indicating if the local statistics are significant at the desired threshold level:
```@example lscor
issignificant(lmguerry, 0.05)
```

As multiple tests are performed on the same dataset, the p-values suffer from the problem of multiple comparisons, leading to many false positives (exceeding the nominal Type I error rate). Multiple solutions have been suggested in the literature to address this issue. Some of them are implemented through the `adjust` parameter. Set the parameter to `:bonferroni` for the conservative Bonferroni approach. For the *False Discovery Rate*, set the parameter to `:fdr`. 

```@example lscor
issignificant(lmguerry, 0.05, adjust = :fdr)
```

## LISA Cluster Map

A cluster map with the significant locations can be plotted with the `plot` function if the [Plots.jl](http://docs.juliaplots.org) package is loaded:
```@example lscor
plot(guerry, lmguerry)
```

The threshold significance value and the adjustment can also be set when plotting a cluster map with the `sig` and `adjust` parameters:
```@example lscor
plot(guerry, lmguerry, sig = 0.05, adjust = :fdr)
```
