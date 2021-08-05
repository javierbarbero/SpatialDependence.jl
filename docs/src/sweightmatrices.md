```@meta
CurrentModule = SpatialDependence
```

# Spatial Weight Matrices

Spatial Weight Matrices can be created from different geometries (polygons, multi polygons, or points) or raw data in Matrix format. By default, spatial weights will be row-standardized so that all rows sum to 1.

## Polygon contiguity

The `polyneigh` function creates a spatial weights object that identifies as neighbors those polygons that are physically contiguous.
```@example polyW
using SpatialDependence
using SpatialDatasets

guerry = sdataset("Guerry")
W = polyneigh(guerry.geometry)
```

By default, the queen contiguity criterion will be used for contiguity. It is possible to specify a Rook criterion by setting the parameter `criterion` to `:Rook`.

```@example polyW
using Plots # hide
pQueen = plot(reggeomlattice(5, 5), [1, 1, 1, 1, 1,  1,  2, 2, 2, 1,  1,  2, 3, 2, 1,  1, 2, 2, 2, 1,  1, 1, 1, 1, 1], Unique(), title = "Queen", legend = false) # hide
pRook  = plot(reggeomlattice(5, 5), [1, 1, 1, 1, 1,  1,  1, 2, 1, 1,  1,  2, 3, 2, 1,  1, 1, 2, 1, 1,  1, 1, 1, 1, 1], Unique(), title = "Rook", legend = false) # hide
plot(pQueen, pRook, layout = @layout([a b])) # hide
```

```@example polyW
W = polyneigh(guerry.geometry, criterion = :Rook)
nothing # hide
```

When the geometry has imperfections and polygons do not touch completely, a tolerance parameter for contiguity can be set with `tol`.

## Points distance threshold

For points geometry, a spatial weights object that identifies as neighbors those points that are below or equal to a specific threshold is computed with the `dnearneigh` function:

```@example distW
using SpatialDependence
using SpatialDatasets

boston = sdataset("Bostonhsg")
W = dnearneigh(boston.geometry, threshold = 4.0)
```

The `dnearneigh` function can also be used with two vectors of coordinates:
```@example distW
W = dnearneigh(boston.x, boston.y, threshold = 4.0)
```

## Points K nearest neighbors

For points geometry, a spatial weights object that identifies the k-nearest neighbors is computed with the `knearneigh` function:

```@example knnW
using SpatialDependence
using SpatialDatasets

boston = sdataset("Bostonhsg")
W = knearneigh(boston.geometry, k = 5)
```

The `knearneigh` function can also be used with two vectors of coordinates:
```@example knnW
W = knearneigh(boston.x, boston.y, k = 5)
```

## Polygon centroids

Polygon centroids and mean centers are computed with the `centroid` and `meancenter` functions.
```@example polyW
cx, cy = centroid(guerry.geometry);
nothing # hide
```

```@example polyW
cmx, cmy = meancenter(guerry.geometry);
nothing # hide
```

Polygon centroids can be used to calculate neighbors based on distance thresholds or k nearest neighbors.
```@example polyW
knearneigh(cx, cy, k = 5)
```

## Spatial Weights from a Matrix

It is possible to create a spatial weights object by supplying a square Matrix. All values that are different from 0 will be neighbors. 

```@example matW
using SpatialDependence # hide
W = SpatialWeights([0 1 0; 1 0 1; 0 1 0])
```

By default, the spatial weights will be row-standardized. However, it is possible to keep the original values in the Matrix as weighs by setting the `standardize` optional parameter to `false`.
```@example matW
W = SpatialWeights([0 1 0; 1 0 1; 0 1 0], standardize = false)
nothing # hide
```

## Spatial Weights transformation

By default, Spatial weight objects are row-standardized, with rows summing 1. It is possible to specify a different transformation with the `wtransform!` and `wtransform` functions. `wtransform!` applies an in-place transformation to the spatial weights object. In contrast,`wtransform` returns a transformed copy of the spatial weights object.

Available transoformations are:

| Code      | Transformation      |
|:----------|:--------------------|
| `:binary` | Binary coding       |
| `:row`    | Row standardization |

```@example polyW
wtransform!(W, :binary)
```

```@example polyW
Wbin = wtransform(W, :binary)
```

```@example polyW
wtransform!(W, :row) # hide
```

## Spatial Weights information

There is a set of functions to obtain information about spatial weights objects.

The function `nobs` returns the number of observations in the spatial weights object. 
```@example polyW
nobs(W)
```

The number of neighbors of each observation is returned with the `cardinalities` function:
```@example polyW
cardinalities(W)
```

A vector of neighbors for particular observations is obtained with the `neighbors` function:
```@example polyW
neighbors(W, 1)
```
And a vector of weights with the `weights` function:
```@example polyW
weights(W, 1)
```

The number of islands - observatiosn with no neighbors - is obtained with the `nislands` function, and the vector of islands with the `islands` function.
```@example polyW
nislands(W)
```

```@example polyW
islands(W)
```

Spatial weights objects implement the following methods to obtain descriptive statistics of the number of neighbors: `minimum`, `maximum`, `mean`, `median`.

## Connectivity Histogram

The `plot` function with a spatial weights object plots a connectivity histogram for the number of neighbors.

```@example polyW
plot(W)
```


## Convert to other objects

Spatial weights objects can be converted to a Matrix using the `Matrix` constructor:
```@example polyW
Matrix(W)
```

and to a sparse matrix with the `sparse` function:
```@example polyW
sparse(W)
```
