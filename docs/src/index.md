```@meta
CurrentModule = SpatialDependence
```

# SpatialDependence

The package [SpatialDependence.jl](https://github.com/javierbarbero/SpatialDependence.jl) is a Julia package for spatial dependence (spatial autocorrelation), spatial weights matrices, and exploratory spatial data analysis (ESDA).

## Example

The following example reads Guerry's Moral statistics of France data and builds a spatial contiguity matrix from the polygons. The spatial weights matrix is row standardized, and the Morans' I index is calculated:

```julia 1
julia> using SpatialDependence
julia> using SpatialDatasets
julia> using StableRNGs

julia> guerry = sdataset("Guerry");

julia> W = polyneigh(guerry.geometry);

julia> wtransform!(W, :row);

julia> moran(guerry.Litercy, W, permutations = 9999, rng = StableRNG(1234567))
Global Moran test of Spatial Autocorrelation
--------------------------------------------

Moran's I: 0.7176053
Expectation:-0.0119048

Randomization test with 9999 permutations.
Standard Error: 0.0707896
zscore: 10.3150637
p-value: 0.0001
```

## Documentation index

```@contents
Pages = ["sweightmatrices.md", "spatiallag.md", "gspatialautocorrelation.md"]
Depth = 3
```
