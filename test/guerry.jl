# Test using the Guerry's Moral statistics of France dataset
@testset "Guerry" begin
    
    guerry = sdataset("Guerry")

    # In this dataset there is no difference between Queen or Rook contiguity
    W = polyneigh(guerry)
    wtransform!(W, :row)

    @test nobs(W) == 85
    @test mean(W) ≈ 4.941176 atol = 1e-5
    @test median(W) ≈ 5.0
    @test minimum(W) == 2
    @test maximum(W) == 8
    @test nislands(W) == 0
    @test islands(W) == []

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

    # ----------------
    #  Plot Recipes
    # ----------------

    # Test Moran's I distribution
    @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), mguerry)

    # Test Geary's c distribution
    @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), cguerry)
end
