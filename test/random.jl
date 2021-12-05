# Test using the Random data
@testset "Random" begin

    @testset "Negative Moran" begin
        geom = reggeomlattice(10, 10)
        W = polyneigh(geom)
        X = rand(StableRNG(1234567), 100)
        Wmat = Matrix(W)
        rho = -0.9

        Y = X + rho * Wmat * X
        mnegative = moran(Y, W, rng = StableRNG(1234567))

        @test score(mnegative) ≈ -0.197453 atol = 1e-5    
        @test mean(mnegative) ≈ -0.010786 atol = 1e-5
        @test std(mnegative) ≈ 0.053380 atol = 1e-5
        @test zscore(mnegative) ≈ -3.496951 atol = 1e-5
        @test pvalue(mnegative) ≈ 0.0001 atol = 1e-5   
    end

end
