@testset "Test Angle Conversion at the Surface" begin
  vθc_deg=Float32.(LinRange(0,360,361))
  vθc_rad=Float32.(LinRange(0,2*pi,361))
  testθ_deg=copy(vθc_deg)
  testθ_rad=copy(vθc_rad)
  θc45_deg=45
  θc135_deg=135
  θc225_deg=225
  θc315_deg=315

  θc45_rad=pi/4
  θc135_rad=3pi/4
  θc225_rad=5pi/4
  θc315_rad=7pi/4

  # check that the angles return the same value after conversion
  θ45_deg=convert_surface_angle_geocentric_to_geodesic_deg(θc45_deg)
  θ135_deg=convert_surface_angle_geocentric_to_geodesic_deg(θc135_deg)
  θ225_deg=convert_surface_angle_geocentric_to_geodesic_deg(θc225_deg)
  θ315_deg=convert_surface_angle_geocentric_to_geodesic_deg(θc315_deg)

  @test θ45_deg≠θc45_deg
  @test θ135_deg≠θc135_deg
  @test θ225_deg≠θc225_deg
  @test θ315_deg≠θc315_deg

  @test convert_surface_angle_geodesic_to_geocentric_deg(θ45_deg)≈θc45_deg
  @test convert_surface_angle_geodesic_to_geocentric_deg(θ135_deg)≈θc135_deg
  @test convert_surface_angle_geodesic_to_geocentric_deg(θ225_deg)≈θc225_deg
  @test convert_surface_angle_geodesic_to_geocentric_deg(θ315_deg)≈θc315_deg

  θ45_rad=convert_surface_angle_geocentric_to_geodesic_rad(θc45_rad)
  θ135_rad=convert_surface_angle_geocentric_to_geodesic_rad(θc135_rad)
  θ225_rad=convert_surface_angle_geocentric_to_geodesic_rad(θc225_rad)
  θ315_rad=convert_surface_angle_geocentric_to_geodesic_rad(θc315_rad)

  @test θ45_rad≠θc45_rad
  @test θ135_rad≠θc135_rad
  @test θ225_rad≠θc225_rad
  @test θ315_rad≠θc315_rad

  @test convert_surface_angle_geodesic_to_geocentric_rad(θ45_rad)≈θc45_rad
  @test convert_surface_angle_geodesic_to_geocentric_rad(θ135_rad)≈θc135_rad
  @test convert_surface_angle_geodesic_to_geocentric_rad(θ225_rad)≈θc225_rad
  @test convert_surface_angle_geodesic_to_geocentric_rad(θ315_rad)≈θc315_rad

  convert_surface_angle_geocentric_to_geodesic_deg!(testθ_deg)
  convert_surface_angle_geocentric_to_geodesic_rad!(testθ_rad)


  @test all(vθc_deg.≠testθ_deg)==false
  @test all(vθc_rad.≠testθ_rad)==false

  convert_surface_angle_geodesic_to_geocentric_deg!(testθ_deg)
  convert_surface_angle_geodesic_to_geocentric_rad!(testθ_rad)

  @test vθc_deg≈ testθ_deg
  @test vθc_rad≈testθ_rad

end
