```@meta
CurrentModule = SpatialDependence
```

# Global Spatial Autocorrelation

Global spatial autocorrelation statistics are used to test for clustering of the data and identify global positive or negative spatial autocorrelation patterns.

Inference is performed through a computational permutation approach by randomly reshuffling the values of the variable to different locations. The observed statistic is compared to a reference distribution under the null hypothesis of spatial randomness.

All the examples on this page assume that the Guerry dataset is loaded and a polygon contiguity spatial weights object has been built and row-standardized.
```@example gscor
using SpatialDatasets
using SpatialDependence
using StableRNGs
using Plots

guerry = sdataset("Guerry")
W = polyneigh(guerry) 
nothing # hide
```

## Moran's I

Moran's I (Moran, 1948) is the most used global spatial autocorrelation statistic. It is computed as:
```math
I = \frac{N}{S_0}\frac{\sum_{i}\sum_{j}w_{ij}(x_{i} - \bar{x})(x_{j} - \bar{x})}{\sum_{i}(x_{i} - \bar{x})^{2}}
```
with $S_0$ being the sum of the spatial weights, $S_0 = \sum_{i}\sum_{j} w_{ij}$.

Moran's I can be computed with the `moran` function. By default, $9,999$ permutations are calculated for the inference. It is possible to specify a different number of permutations with the `permutations` optional parameter. For reproduciibility, it is possible to specify a custom random number generator with the `rng` optional parameter.
```@example gscor
moran(guerry.Litercy, W, permutations = 9999, rng = StableRNG(1234567))
```

The interpretation of Moran's I depends on its value and significance:

| Moran's I | z-value                 | Interpretation                   |
|:----------|:------------------------|:---------------------------------|
| $> 0$     | $> 0$ and significant   | Positive spatial autocorrelation |
| $< 0$     | $< 0$ and significant   | Negative spatial autocorrelation |
| Any       | Any and non-significant | Spatially random                 |

## Geary's c

Geary's c (Geary, 1954) is a global spatial autocorrelation statistic that focuses on dissimilarity. It is computed as:
```math
c = \frac{N - 1}{2S_0}\frac{\sum_{i}\sum_{j}w_{ij}(x_{i} -x_{j})^2}{\sum_{i}(x_{i} - \bar{x})^{2}}
```

Geary's I can be computed with the `geary` function. By default, $9,999$ permutations are calculated for the inference. It is possible to specify a different number of permutations with the `permutations` optional parameter. For reproduciibility, it is possible to specify a custom random number generator with the `rng` optional parameter.
```@example gscor
geary(guerry.Litercy, W, permutations = 9999, rng = StableRNG(1234567))
```

The interpretation of Gery's c depends on its value and significance:

| Geary's c | z-value                 | Interpretation                   |
|:----------|:------------------------|:---------------------------------|
| $< 1$     | $< 0$ and significant   | Positive spatial autocorrelation |
| $> 1$     | $> 0$ and significant   | Negative spatial autocorrelation |
| Any       | Any and non-significant | Spatially random                 |

## Reference distribution

The random permutation operation results in a reference distribution for the statistic under the null hypothesis of spatial randomness. If the [Plots.jl](http://docs.juliaplots.org) package is loaded, it is possible to plot the reference distribution together with the observed statistic (vertical red line) using the `plot` function:

```@example gscor
Ilitercy = moran(guerry.Litercy, W, permutations = 9999, rng = StableRNG(1234567))
plot(Ilitercy)
```

## Moran Scatter Plot

The Moran Scatter Plot (Anselin, 1996) is a scatterplot with the variable in the horizontal axis and its spatial lag on the vertical axis. 

If the [Plots.jl](http://docs.juliaplots.org) package is loaded, the `plot` function can be used to plot a Moran scatter plot. By default, the values of the variable are z-standardized, but it is possible to build the plot without standardizing by setting the optional parameter `standardize` to `false`.

```@example gscor
plot(guerry.Litercy, W)
```

In the Moran scatter plot, observations are located in four quadrants, depending on the value of the attribute and the value of their neighbors with respect to the mean:

| Quadrants    | Spatial Autocorrelation | Interpretation                         |
|:-------------|:------------------------|:---------------------------------------|
| Upper right  | Positive: High-high     | High values surrounded by high values  |
| Lower left   | Positive: Low-low       | Low values surrounded by low values    |
| Lower right  | Negative: High-Low      | High values surrounded by low values   |
| Upper left   | Negative: Low-high      | Low values surrounded by high values   |
