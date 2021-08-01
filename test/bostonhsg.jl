# Test using the Boston housing and neighborhood dataset
@testset "Bostonhsg" begin
    
    boston = sdataset("Bostonhsg")

    # ----------------
    #  Distance  with 3
    # ----------------
    Wdist3 = dnearneigh(boston, threshold = 3)
    wtransform!(Wdist3, :row)   

    @test nobs(Wdist3) == 506
    @test mean(Wdist3) ≈ 44.711462 atol = 1e-5
    @test median(Wdist3) ≈ 29.0
    @test minimum(Wdist3) == 0
    @test maximum(Wdist3) == 122
    @test nislands(Wdist3) == 4
    @test islands(Wdist3) == [55, 65, 286, 287]

    # Moran's I
    mboston = moran(boston.MEDV, Wdist3, permutations = 0)

    # ----------------
    #  Distance  with 4
    # ----------------
    Wdist4 = dnearneigh(boston, threshold = 4)
    wtransform!(Wdist4, :row)

    @test nobs(Wdist4) == 506
    @test mean(Wdist4) ≈ 72.225296 atol = 1e-5
    @test median(Wdist4) ≈ 56.0
    @test minimum(Wdist4) == 1
    @test maximum(Wdist4) == 178
    @test nislands(Wdist4) == 0
    @test islands(Wdist4) == []

    # Moran's I with 999 permutations
    mboston = moran(boston.MEDV, Wdist4, rng = StableRNG(1234567), permutations = 999)

    @test score(mboston) ≈ 0.267418 atol = 1e-5
    @test pvalue(mboston) ≈ 0.001 atol = 1e-5
    @test std(mboston) ≈ 0.013811 atol = 1e-5
    @test zscore(mboston) ≈ 19.493602 atol = 1e-5

    # ----------------
    #  KNN with k = 10
    # ----------------
    Wknn10 = knearneigh(boston, k = 10)
    wtransform!(Wknn10, :row)

    @test nobs(Wknn10) == 506
    @test mean(Wknn10) == 10
    @test median(Wknn10) ≈ 10.0
    @test minimum(Wknn10) == 10
    @test maximum(Wknn10) == 10
    @test nislands(Wknn10) == 0

    # Moran I with 999 permutations
    mboston = moran(boston.MEDV, Wknn10, rng = StableRNG(1234567), permutations = 999)

    @test score(mboston) ≈ 0.530743 atol = 1e-5
    @test std(mboston) ≈ 0.018956 atol = 1e-5
    @test zscore(mboston) ≈ 28.080397 atol = 1e-5
    @test pvalue(mboston) ≈ 0.001 atol = 1e-5

    # Moran's I with 9999 permutations
    mboston = moran(boston.MEDV, Wknn10, rng = StableRNG(1234567), permutations = 9999)

    @test score(mboston) ≈ 0.530743 atol = 1e-5
    @test std(mboston) ≈ 0.018541 atol = 1e-5
    @test zscore(mboston) ≈ 28.727399 atol = 1e-5
    @test pvalue(mboston) ≈ 0.0001 atol = 1e-5

    # Geary's c I with 9999 permutations
    cboston = geary(boston.MEDV, Wknn10, rng = StableRNG(1234567), permutations = 9999)

    @test score(cboston) ≈ 0.432582 atol = 1e-5
    @test std(cboston) ≈ 0.022467 atol = 1e-5
    @test zscore(cboston) ≈ -25.249835 atol = 1e-5
    @test pvalue(cboston) ≈ 0.0001 atol = 1e-5

    # ----------------
    #  Plot Recipes
    # ----------------

    # Test Moran's I distribution
    @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), mboston)

    # Test Moran's Plot
    @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), boston.MEDV, Wknn10)
    @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}([(:standardize, false)]), boston.MEDV, Wknn10)

end
