# SpatialDependence.jl

A Julia package for exploratory spatial data analysis (ESDA), including functions for spatial weights matrices creation, testing for spatial dependence (spatial autocorrelation), and choropleth maps.

| Documentation | Build Status      | Coverage    |
|:-------------:|:-----------------:|:-----------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] |  [![][githubci-img]][githubci-url] | [![][codecov-img]][codecov-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://javierbarbero.github.io/SpatialDependence.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://javierbarbero.github.io/SpatialDependence.jl/dev

[githubci-img]: https://github.com/javierbarbero/SpatialDependence.jl/workflows/CI/badge.svg
[githubci-url]: https://github.com/javierbarbero/SpatialDependence.jl/actions

[codecov-img]: https://codecov.io/gh/javierbarbero/SpatialDependence.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/javierbarbero/SpatialDependence.jl

The package **SpatialDependence.jl** contains functions to create and handle spatial weights matrices from polygon and point geometries. It also has functions for calculating spatial lags, test for spatial autocorrelation, and plotting choropleth maps.

The package takes advantage of Julia multi-threading features to increase performance if Julia is started with multiple threads.

See the last stable version [documentation][docs-stable-url] to learn how to use it.

The package is currently under heavy development, with more functionality coming soon.

## Installation

The package can be installed with the Julia package manager:
```julia
julia> using Pkg; Pkg.add("SpatialDependence")
```

## Main functionality

**Spatial Weights Matrices**

- Polygon contiguity: Queen and Rook.
- Points distance threshold.
- Points K nearest neighbors.

**Global Spatial Autocorrelation**

- Moran's I.
- Geary's c.

**Choropleth Maps**

- Graduated: Equal Intervals, Quantiles, Natural Breaks, Custom Breaks.
- Statistical: Box Plot, Standard Deviation, Percentiles.
- Unique Values.
- Co-location.

## Author

SpatialDependence.jl is being developed by [Javier Barbero](http://www.javierbarbero.net).
