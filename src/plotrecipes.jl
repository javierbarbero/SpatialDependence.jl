# Plot recipes

# Recipe for Moran plot
@recipe function f(x::Vector{Float64}, W::SpatialWeights, zstandardize::Bool = false)

    Wx = slag(W, x)

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

# Recipe for Moran's I distribution
@recipe function plot(x::GlobalMoran)

    I = x.I
    Iperms = x.Iperms

    legend --> false

    # Histogram of permutated values
    @series begin
        seriestype := :histogram
        Iperms
    end

    # Vertical line at Moran's i
    @series begin
        seriestype := :vline
        seriescolor --> :red
        [I]
    end

end
