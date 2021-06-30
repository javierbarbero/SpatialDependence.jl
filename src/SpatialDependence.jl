module SpatialDependence

    """
    SpatialDependence
    A Julia package for spatial dependence (spatial autocorrelation), spatial weights matrices, and exploratory spatial data analysis (ESDA).
    [SpatialDependence repository](https://github.com/javierbarbero/SpatialDependence.jl).
    """    

    using GeoInterface: Point, coordinates
    using NearestNeighbors: KDTree, knn, inrange
    using Random: shuffle, AbstractRNG, default_rng
    using RecipesBase

    import SparseArrays: SparseMatrixCSC, sparse
    import StatsBase: weights, nobs, minimum, maximum, median, score, zscore, standardize, ZScoreTransform
    import Statistics: mean, median, std
    
    export 
        # Types
        SpatialWeights,

        # Spatial Weights Creation functions
        dnearneigh,
        knearneigh,
        polyneigh,

        # Spatial Weights functions
        neighbors, 
        weights, 
        nislands, 
        islands,
        wtransform, 
        wtransform!, 
        wtransformation,
        slag,
        sparse,

        # Spatial autocorrelation
        GlobalMoran,
        moran,
        
        # Statistics
        nobs,
        mean,
        median,
        std,
        score,
        zscore,
        pvalue

    include("sweights/sweights.jl")
    include("sweights/dnearneigh.jl")
    include("sweights/knearneigh.jl")
    include("sweights/polyneigh.jl")

    include("slag.jl")

    include("scor/moran.jl")
    include("scor/moranplot.jl")

end
