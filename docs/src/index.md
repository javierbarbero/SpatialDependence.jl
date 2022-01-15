```@meta
CurrentModule = SpatialDependence
```

# SpatialDependence

The package [SpatialDependence.jl](https://github.com/javierbarbero/SpatialDependence.jl) is a Julia package for exploratory spatial data analysis (ESDA), including functions for spatial weights matrices creation, testing for spatial dependence (spatial autocorrelation), and choropleth mapping.

## Installation

The package can be installed with the Julia package manager:
```julia
julia> using Pkg; Pkg.add("SpatialDependence")
```

## Example

The following example reads Guerry's Moral statistics of France data, builds a spatial contiguity matrix from the polygons, and calculates the Morans' I global spatial autocorrelation statistic:

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
# Build polygon contiguity matrix
W = polyneigh(guerry.geometry);
```

```@example intro
# Global Moran test of Spatial Autocorrelation of the Litercy variable
moran(guerry.Litercy, W, permutations = 9999, rng = StableRNG(1234567))
```

```@example intro
# Moran Scatterplot of the Litercy variable
plot(guerry.Litercy, W, xlabel = "Litercy")
```

```@example intro
# Local Indicators of Spatial Association (LISA) - Local Moran
lmguerry = localmoran(guerry.Litercy, W, permutations = 9999, rng = StableRNG(1234567))
```

```@example intro
# LISA Cluster Map
plot(guerry, lmguerry, sig = 0.05, adjust = :fdr)
```

## Documentation index

```@contents
Pages = ["sweightmatrices.md", "spatiallag.md", "gspatialautocorrelation.md", "lspatialautocorrelation.md", "choropleth.md", "parallelcomputing.md", "bibliography.md"]
Depth = 2
```

## Authors

SpatialDependence.jl is being developed by [Javier Barbero](http://www.javierbarbero.net).
