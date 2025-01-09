using LinearAlgebra
using StaticArrays

@testset "Creation Rays and Ellipsoid" begin
  @testset "Creation of rays" begin
    r=Ray2D(SVec2(1.0,1.0),SVec2(0.0,2.0))
    @test r.origin==SVec2(1.0,1.0)
    @test r.direction==SVec2(0.0,1.0)

    @test_throws ArgumentError r_broken=Ray2D(SVec2(1.0,1.0),SVec2(0.0,0.0))
  end
  @testset "Creation of ellipsoids" begin
    radii=[3.0,2.0]
    e=Ellipsoid(radii...)

    @test e.scale==1.0./radii
    @test e.unscale==radii


    e=Ellipsoid(radii[2],radii[1])

    @test e.scale==1.0./radii
    @test e.unscale==radii
  end
end
#@testset "Distance" begin
  @testset "Distance from unit circle" begin
    r=Ray2D(SVec2(1.0,1.0),SVec2(-1.0,-1.0))
    @test distance_from_unit_circle(r)≈(sqrt(2)-1)
    r=Ray2D(SVec2(1.0,1.0),SVec2(1.0,-1.0))
    @test distance_from_unit_circle(r)==Inf
  end
  @testset "Distance from segment" begin
    r=Ray2D(SVec2(1.0,1.0),SVec2(1.0,0.0))
    @test distance_from_segment(r,SVec2(2.0),SVec2(2.0,2.0))≈1.0
    @test distance_from_segment(r,SVec2(1.01,0.0),SVec2(1.01,1.0))≈0.01
    @test distance_from_segment(r,SVec2(0.0,0.0),SVec2(0.0,2.0))==Inf
  end
  @testset "Distante from radius" begin
    for θ in 1:60
      match_at_1=convert(ECEF2D{NormalizedEarth},LLA2D{NormalizedEarth}(1.,θ)) |> x-> SVec2([x.w,x.z])
      direction_to_radii=SVec2(-1.0,1.0) |> x->x/hypot(x...)

      r=Ray2D(match_at_1.-direction_to_radii,direction_to_radii)
      my_radius_θ= create_radii_from_θ(θ)
      @test distance_from_radii(r,my_radius_θ)≈1.0
    end
  end
