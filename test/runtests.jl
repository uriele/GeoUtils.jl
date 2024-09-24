using GeoUtils
using Test
using Aqua

@testset "GeoUtils.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(GeoUtils)
    end
    # Write your tests here.
end
