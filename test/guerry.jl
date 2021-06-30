# Test using the Guerry's Moral statistics of France dataset
@testset "Guerry" begin
    
    guerry = sdataset("Guerry")

    # In this dataset there is no difference between Queen or Rook contiguity
    W = polyneigh(guerry.geometry)
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
    @test pvalue(mguerry) ≈ 0.0001 atol = 1e-5
    @test std(mguerry) ≈ 0.070790 atol = 1e-5
    @test zscore(mguerry) ≈ 10.315064 atol = 1e-5

end
