# Test using the Guerry's Moral statistics of France dataset
@testset "Guerry" begin
    
    guerry = sdataset("Guerry")

    # In this dataset there is no difference between Queen or Rook contiguity
    W = polyneigh(guerry)

    @test nobs(W) == 85
    @test mean(W) ≈ 4.941176 atol = 1e-5
    @test median(W) ≈ 5.0
    @test minimum(W) == 2
    @test maximum(W) == 8
    @test nislands(W) == 0
    @test islands(W) == []

    # Spatial lag
    @test W * guerry.Litercy == slag(W, guerry.Litercy)

    # Moran's I
    mguerry = moran(guerry.Litercy, W, rng = StableRNG(1234567), permutations = 9999)

    @test score(mguerry) ≈ 0.717605 atol = 1e-5    
    @test length(scoreperms(mguerry)) == 9999
    @test mean(mguerry) ≈ -0.012594 atol = 1e-5
    @test std(mguerry) ≈ 0.070790 atol = 1e-5
    @test zscore(mguerry) ≈ 10.315064 atol = 1e-5
    @test pvalue(mguerry) ≈ 0.0001 atol = 1e-5   
    @test SpatialDependence.expected(mguerry) == - 1 / (length(guerry.Litercy) - 1)
    @test SpatialDependence.testname(mguerry) == "Moran's I"

    @test_nowarn show(IOBuffer(), mguerry)

    # Geary's c
    cguerry = geary(guerry.Litercy, W, rng = StableRNG(1234567), permutations = 9999)

    @test score(cguerry) ≈ 0.250202 atol = 1e-5
    @test length(scoreperms(cguerry)) == 9999
    @test mean(cguerry) ≈ 1.000411 atol = 1e-5
    @test std(cguerry) ≈ 0.072280 atol = 1e-5
    @test zscore(cguerry) ≈ -10.379150 atol = 1e-5
    @test pvalue(cguerry) ≈ 0.0001 atol = 1e-5
    @test SpatialDependence.expected(cguerry) == 1
    @test SpatialDependence.testname(cguerry) == "Geary's c"

    @test_nowarn show(IOBuffer(), cguerry)

    # Polygon mean center
    cx, cy = meancenter(guerry.geometry)
    @test cx[1] ≈ 836732.867469879565760 atol = 1e-5
    @test cy[1] ≈ 2125844.293172690551728 atol = 1e-5
    @test cx[50] ≈ 814647.356215213309042 atol = 1e-5
    @test cy[50] ≈ 2350156.217068645637482 atol = 1e-5

    # Polygon centroid
    cx, cy = centroid(guerry.geometry)
    @test cx[1] ≈ 832852.278780025895685 atol = 1e-5
    @test cy[1] ≈ 2126600.575969800818712 atol = 1e-5
    @test cx[50] ≈ 815193.599669872317463 atol = 1e-5
    @test cy[50] ≈ 2349544.092820507008582 atol = 1e-5

    # Moran with 5 k nearest neighbors
    Wk5 = knearneigh(cx, cy, k = 5)
    morank5 = moran(guerry.Litercy, Wk5)
    @test score(morank5) ≈ 0.678857 atol = 1e-5

    # ----------------
    #  Plot Recipes
    # ----------------

    # Test Moran's I distribution
    @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), mguerry)

    # Test Geary's c distribution
    @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), cguerry)
end
