using SpatialDependence
using GeoInterface
using GeoInterface
using RecipesBase
using SpatialDatasets
using StableRNGs
using Test

@testset "SpatialDependence.jl" begin
    
    include("sweights.jl")
    include("guerry.jl")
    include("bostonhsg.jl")

end