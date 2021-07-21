```@meta
CurrentModule = SpatialDependence
```

# Choropleth Maps

The package contains functions to plot Choropleth Maps using different classification methods.

All the examples on this page asummes that the Guerry dataset and the [Plots.jl](http://docs.juliaplots.org) package is loaded.
```@example choropleth
using SpatialDependence
using SpatialDatasets
using Plots

guerry = sdataset("Guerry"); 
nothing # hide
```

Choropleth maps can be made with `plot` function, with the first parameter being a vector of geometries, the second parameter the variable to classify, and the third parameter the classification method to use:
```@example choropleth
plot(guerry.geometry, guerry.Litercy, EqualIntervals())
nothing # hide
```

It is also possible to pass a Table or a DataFrame as the first parameter and a symbol as the second parameter:
```@example choropleth
plot(guerry, :Litercy, EqualIntervals())
nothing # hide
```

## Graduated

### Equal Intervals

Spatial observations are grouped into categories of the same length:
```@example choropleth
plot(guerry, :Litercy, EqualIntervals())
```

It is possible to specify a different number of intervals:
```@example choropleth
plot(guerry, :Litercy, EqualIntervals(3))
nothing # hide
```

### Quantiles

Spatial observations are grouped according to quantiles:
```@example choropleth
plot(guerry, :Litercy, Quantiles())
```

By default, spatial observations are categorized in quintiles (5 quantiles). It is possible to specify a different number of quantiles:
```@example choropleth
plot(guerry, :Litercy, Quantiles(4))
nothing # hide
```

### Natural Breaks (Jenks)

The Natural Breaks (Jenks) algorithm classify observations into groups with low heterogeneity within groups and high heterogeneity between groups:

```@example choropleth
plot(guerry, :Litercy, NaturalBreaks())
```

Is it possible to specify a different number of groups:
```@example choropleth
plot(guerry, :Litercy, NaturalBreaks(7))
nothing # hide
```

### Custom Breaks

Users can specify the bins to generate custom breaks. The lower bound of the first class is set at the minimum value and the upper bound of the last class at the maximum value.
```@example choropleth
plot(guerry, :Litercy, CustomBreaks([20, 30, 60]))
```

## Statistical

### Box Plot (Box map)

This method classify observations as in a box plot:

|Category | Lower | Upper |
|---------|-------|-------|
| 1       | -$\infty$ | $Q_{25} - h * IQR$ |
| 2       | $Q_{25} - h * IQR$ | $Q_{25}$ |
| 3       | $Q_{25}$ | $Q_{50}$ |
| 4       | $Q_{50}$ | $Q_{75}$ |
| 5       | $Q_{75}$ | $Q_{75} + h * IQR$ |
| 6       | $Q_{75} + h * IQR$ | $\infty$ |

where $Q_{25}$, $Q_{50}$ and $Q_{75}$, are the first, second and third quartile respectively. $IQR$ is the interquartile range, and $h$ is the hinge.

```@example choropleth
plot(guerry, :Litercy, BoxPlot())
```

By default, the hinge is set to $1.5$, but the user can choose a different value:

```@example choropleth
plot(guerry, :Litercy, BoxPlot(3.0))
nothing # hide
```

### Standard Deviation

In this method, observations are z-standardization and classified as standard deviations from the mean:

|Category | Lower     | Upper    |
|---------|-----------|----------|
| 1       | -$\infty$ | $-2$     |
| 2       | $-2$      | $-1$     |
| 3       | $-1$      | $0$      |
| 4       | $0$       | $1$      |
| 5       | $1$       | $2$      |
| 6       | $2$       | $\infty$ |

```@example choropleth
plot(guerry, :Litercy, StdMean())
```


### Percentiles

Observations are classified according to the following percentiles:

|Category | Lower     | Upper    |
|---------|-----------|----------|
| 1       | minimum   | $1\%$    |
| 2       | $1\%$     | $10\%$   |
| 3       | $10\%$    | $50\%$   |
| 4       | $50\%$    | $90\%$   |
| 5       | $90\%$    | $99\%$   |
| 6       | $99\%$    | maximum  |

```@example choropleth
plot(guerry, :Litercy, Percentiles())
```

It is possible to specify different percentiles:
```@example choropleth
plot(guerry.geometry, guerry.Litercy, Percentiles([10, 30, 50, 70, 90]))
nothing # hide
```

## Unique Values

In a Unique Values map, each unique value is assigned to a different category.

```@example choropleth
plot(guerry, :Region, Unique())
```

## Co-location

Co-location maps are used to identify spatial observations where the categories of two or more variables match.

```@example choropleth
plot(guerry, [:Litercy, :Donatns], CoLocation(BoxPlot()))
```

## Customizing Maps

By default, the lower interval is **open** (value not included), and the upper interval is **closed** (value is included). It is possible to change this behavior by setting the `lower` and `upper` parameters to `:open` or `:closed`.

```@example choropleth
plot(guerry, :Litercy, Quantiles(), lower = :closed, upper = :open)
```

**Total counts** are automatically added to the labels of the categories in the map. It is possible to suppress counts by setting the `counts` parameter to `false`.

```@example choropleth
plot(guerry, :Litercy, Quantiles(), counts = false)
```

In graduated maps, the lower and upper bound are separated with a comma. It is possible to specify a different **separator** with the `sep` parameter:

```@example choropleth
plot(guerry, :Litercy, Quantiles(), sep = " -- ")
```

It is possible to change the **legend** position with the `legend` parameter. Possible values are: :best, `:topleft`, `:top`, `:topright`, `:right`, `:bottomright`, `:bottom`,  `:bottomleft`, `:left`. It also possible to place the legend outside the plot area with: `:outertopleft`, `:outertop`, `:outertopright`, `:outerright`, `:outerbottomright`, `:outerbottom`, `:outerbottomleft`. The legend can be removed by setting the parameter to `false`.

```@example choropleth
plot(guerry, :Litercy, Quantiles(), legend = :topleft)
```

**Color palette** is set by default to `:YlOrBr` for graduated maps, reverse `:RdBu` for statistical maps, and  `:Paired` for unique values maps. Users can change the palette with the `palette` parameter and reverse the colors by setting `rev` to `true`. A list of the color shcemes available in the Plots.jl package can be found [here](http://docs.juliaplots.org/latest/generated/colorschemes/#Pre-defined-ColorSchemes).

```@example choropleth
plot(guerry, :Litercy, Quantiles(), palette = :greens, rev = true)
```

## Classify without mapping

Users interested only in the classification and not in the choropleth map can use the functions `mapclassify` to obtain the classification of observations according to the criteria listed before.

```@example choropleth
mc = mapclassify(Quantiles(), guerry.Litercy);
nothing # hide
```

It is possible to obtain the number of observations on each class with the `counts` function:
```@example choropleth
counts(mc)
```

The lower and upper bounds can be retrieved with the `bounds` function:
```@example choropleth
lbound, ubound = bounds(mc)
```

With the `maplabels` function, it is possible to generate the labels of the classes:
```@example choropleth
labs = maplabels(mc)
```

We can obtain the class to which each observation is assigned with the `assignments` function:
```@example choropleth
assignments(mc)
```

For unique values and Co-location classifications, the levels of the classes can be obtained with the `levels` class:
```@example choropleth
mc = mapclassify(Unique(), guerry.Region)
levels(mc)
```
