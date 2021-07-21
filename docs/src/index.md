```@meta
CurrentModule = SpatialDependence
```

# SpatialDependence

The package [SpatialDependence.jl](https://github.com/javierbarbero/SpatialDependence.jl) is a Julia package for spatial dependence (spatial autocorrelation), spatial weights matrices, and exploratory spatial data analysis (ESDA).

## Installation

The package can be installed with the Julia package manager:
```julia
julia> using Pkg; Pkg.add("SpatialDependence")
```

## Example

The following example reads Guerry's Moral statistics of France data and builds a spatial contiguity matrix from the polygons. The spatial weights matrix is row standardized, and the Morans' I index is calculated:

```@example intro
# Load packages
using Plots
using SpatialDependence
using SpatialDatasets
using StableRNGs

# Guerry's Moral statistics of France data from the SpatialDatasets.jl package
guerry = sdataset("Guerry");

# Plot Litercy variable
plot(guerry, :Litercy, NaturalBreaks(), legend = :topleft, title = "Litercy")
```

```@example intro
# Build polygon contiguity matrix
W = polyneigh(guerry.geometry);
```

```@example intro
# Row-standardize the spatial weights matrix
wtransform!(W, :row);

# Global Moran test of Spatial Autocorrelation of the Litercy variable
moran(guerry.Litercy, W, permutations = 9999, rng = StableRNG(1234567))
```

```@example intro
# Moran Scatterplot of the Litercy variable
plot(guerry.Litercy, W, true, xlabel = "Litercy")
```

## Documentation index

```@contents
Pages = ["sweightmatrices.md", "spatiallag.md", "choropleth.md", "gspatialautocorrelation.md"]
Depth = 2
```

## Authors

SpatialDependence.jl is being developed by [Javier Barbero](http://www.javierbarbero.net).
