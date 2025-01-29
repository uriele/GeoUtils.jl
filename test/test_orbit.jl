#using LinearAlgebra: Diagonal
using StaticArrays
using CoordinateTransformations

@testset "Creation Rays and Ellipsoid" begin
  @testset "Creation of rays" begin
    r=Ray2D(Vec2(1.0,1.0),Vec2(0.0,2.0))
    @test r.origin==Vec2(1.0,1.0)
    @test r.direction==Vec2(0.0,1.0)

    @test_throws ArgumentError r_broken=Ray2D(Vec2(1.0,1.0),Vec2(0.0,0.0))
  end
  @testset "Creation of ellipsoids" begin
    e=Ellipsoid(2.0,3.0)
    @test e.radius==Vec2(2.0,3.0)
    @test e.scale([1. 0;0. 1.])==[1/2 0.0;0.0 1/3]
    @test e.unscale([1. 0;0. 1.])==[2.0 0.0;0.0 3.0]
  end
end
#@testset "Distance" begin
  @testset "Distance from unit circle" begin
    r=Ray2D(Vec2(1.0,1.0),Vec2(-1.0,-1.0))
    @test distance_from_unit_circle(r)≈(sqrt(2)-1)
    r=Ray2D(Vec2(1.0,1.0),Vec2(1.0,-1.0))
    @test distance_from_unit_circle(r)==Inf
  end
  @testset "Distance from segment" begin
    r=Ray2D(Vec2(1.0,1.0),Vec2(1.0,0.0))
    @test distance_from_segment(r,Vec2(2.0),Vec2(2.0,2.0))≈1.0
    @test distance_from_segment(r,Vec2(1.01,0.0),Vec2(1.01,1.0))≈0.01
    @test distance_from_segment(r,Vec2(0.0,0.0),Vec2(0.0,2.0))==Inf
  end
  @testset "Distante from radius" begin
    for θ in 1:60
      match_at_1=convert(ECEF2D{NormalizedEarth},LLA2D{NormalizedEarth}(1.,θ)) |> x-> Vec2([x.w,x.z])
      direction_to_radii=Vec2(-1.0,1.0) |> x->x/hypot(x...)

      r=Ray2D(match_at_1.-direction_to_radii,direction_to_radii)
      my_radius_θ= create_radii_from_θ(θ)
      @test distance_from_radii(r,my_radius_θ)≈1.0
    end
  end
