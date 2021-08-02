```@meta
CurrentModule = SpatialDependence
```

# Spatial lags

The spatial lag of a variable is calculated with the values of the neighboring observations. If the spatial weights matrix is row-standardized, the spatial lag is the average value of the neighbors.

Spatial lags are calculated with the `slag` function or using the `*` operator with a spatial weights objects in the first position.
```@example slag
using SpatialDependence # hide
using SpatialDatasets # hide
guerry = sdataset("Guerry") # hide
W = polyneigh(guerry)
slag(W, guerry.Litercy)
```

```@example slag
W * guerry.Litercy
```
