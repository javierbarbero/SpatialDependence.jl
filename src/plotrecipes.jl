# Plot recipes

# Recipe for Moran plot
@recipe function f(x::Vector{Float64}, W::SpatialWeights)

    Wx = slag(W, x)

    zstandardize  = get(plotattributes, :standardize, true)
    if zstandardize
        z = standardize(ZScoreTransform, x)
        Wz = standardize(ZScoreTransform, Wx)
    else
        z = x
        Wz = Wx
    end

    title --> "Moran Scatterplot"
    xguide --> "Attribute"
    yguide --> "Spatial Lag of " * plotattributes[:xguide]
    legend --> false
    grid --> false

    # Variable and spatial lag of variable
    @series begin
        seriestype := :scatter
        smooth := true
        z, Wz
    end

    # Vertical line at variable mean
    zmean = mean(z)
    @series begin
        seriestype := :vline
        linestyle := :dash
        seriescolor --> :red
        [zmean]
    end

    # Horizontal line at W * variable mean
    Wzmean = mean(Wz)
    @series begin
        seriestype := :hline
        linestyle := :dash
        seriescolor --> :red
        [Wzmean]
    end

end

# Recipe for Global Spatial Autocorrelation reference distribution
@recipe function plot(x::AbstractGlobalSpatialAutocorrelation)

    S = score(x)
    Sperms = scoreperms(x)

    legend --> false

    # Histogram of permutated values
    @series begin
        seriestype := :histogram
        Sperms
    end

    # Vertical line at Moran's i
    @series begin
        seriestype := :vline
        seriescolor --> :red
        [S]
    end

end
