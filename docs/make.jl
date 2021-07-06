using Documenter
using Plots
using SpatialDatasets
using SpatialDependence
using StableRNGs

DocMeta.setdocmeta!(SpatialDependence, :DocTestSetup, :(using SpatialDependence); recursive=true)

makedocs(;
    modules=[SpatialDependence],
    authors="Javier Barbero and contributors",
    repo="https://github.com/javierbarbero/SpatialDependence.jl/blob/{commit}{path}#{line}",
    sitename="SpatialDependence.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://javierbarbero.github.io/SpatialDependence.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Spatial Weight Matrices" => "sweightmatrices.md",
        "Spatial Lag" => "spatiallag.md",
        "Global Spatial Autocorrelation" => "gspatialautocorrelation.md",
    ],
)

deploydocs(;
    repo="github.com/javierbarbero/SpatialDependence.jl",
    devbranch = "main"
)
