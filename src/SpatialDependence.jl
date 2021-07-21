module SpatialDependence

    """
    SpatialDependence
    A Julia package for spatial dependence (spatial autocorrelation), spatial weights matrices, and exploratory spatial data analysis (ESDA).
    [SpatialDependence repository](https://github.com/javierbarbero/SpatialDependence.jl).
    """    

    using GeoInterface: AbstractGeometry, AbstractPoint, AbstractPolygon, AbstractMultiPolygon, coordinates, geotype, Point, shapecoords
    using NearestNeighbors: KDTree, knn, inrange
    using PlotUtils: palette
    using Random: shuffle, AbstractRNG, default_rng
    using RecipesBase
    using Tables: istable

    import Base: length
    import SparseArrays: SparseMatrixCSC, sparse
    import StatsBase: counts, levels, maximum, median, minimum, nobs, percentile, score, standardize, weights, zscore, ZScoreTransform
    import Statistics: mean, median, std, quantile
    
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
        pvalue,

        # Choropleth Maps and Classification<
        AbstractMapClassificator,
        AbstractGraduatedMapClassificator,
        EqualIntervals,
        Quantiles,
        CustomBreaks,
        GraduatedMapClassification,

        AbstractStatisticalMapClassificator,
        NaturalBreaks,
        BoxPlot,
        StdMean,
        Percentiles,

        AbstractUniqueMapClassificator,
        Unique,
        UniqueMapClassification,

        AbstractCoLocationMapClassificator,
        CoLocation,
        CoLocationMapClassification,

        mapclassify,
        maplabels,

        assignments,
        counts,
        bounds,
        levels

    include("sweights/sweights.jl")
    include("sweights/dnearneigh.jl")
    include("sweights/knearneigh.jl")
    include("sweights/polyneigh.jl")

    include("slag.jl")

    include("scor/moran.jl")

    include("maps/classification.jl")
    include("maps/graduated.jl")
    include("maps/statistical.jl")
    include("maps/unique.jl")
    include("maps/colocation.jl")

    include("plotrecipes.jl")
    include("plotchoropleth.jl")

end
