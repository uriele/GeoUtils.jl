using GeoUtils
using Test
using Aqua
const alt_99= [120.0, 119.0, 118.0, 117.0, 116.0, 115.0, 114.0, 113.0, 112.0, 111.0,
  110.0, 109.0, 108.0, 107.0, 106.0, 105.0, 104.0, 103.0, 102.0, 101.0, 100.0, 99.0, 98.0,
  97.0, 96.0, 95.0, 94.0, 93.0, 92.0, 91.0, 90.0, 89.0, 88.0, 87.0, 86.0, 85.0, 84.0,
  83.0, 82.0, 81.0, 80.0, 79.0, 78.0, 77.0, 76.0, 75.0, 74.0, 73.0, 72.0, 71.0, 70.0,
  69.0, 68.0, 67.0, 66.0, 65.0, 64.0, 63.0, 62.0, 61.0, 60.0, 59.0, 58.0, 57.0, 56.0,
  55.0, 54.0, 53.0, 52.0, 51.0, 50.0, 49.0, 48.0, 47.0, 46.0, 45.0, 44.0, 43.0, 42.0,
  41.0, 40.0, 39.0, 38.0, 37.0, 36.0, 35.0, 34.0, 33.0, 32.0, 31.0, 30.0, 29.0, 28.0,
  27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0, 20.0, 19.02887, 16.73447, 15.08326, 13.78967,
  12.73618, 11.86036, 11.12134, 10.46806, 9.81183, 9.07048, 8.43950, 7.94138, 7.40065,
  6.84547, 6.29096, 5.74693, 5.21813, 4.71070, 4.23089, 3.77917, 3.35240, 2.95164, 2.58209,
  2.24559, 1.93902, 1.66020, 1.41141,1.19548, 1.01039, 0.85022, 0.70959, 0.58788, 0.48020,
  0.38196, 0.29295, 0.21291, 0.14163, 0.07936, 0.00000];
# Do this test only if the data folder is present in the root directory
if isdir("data")
  @testset "Test utils on data" begin
    @testset "Test all_files_are_same" begin
      @test GeoUtils.all_files_are_same("data","in_lat.dat")==true
      @test GeoUtils.all_files_are_same("data","in_alt.dat")==true
      @test GeoUtils.all_files_are_same("data","in_lon.dat")==false
    end

    @testset "Test convert_to_array" begin
      @test GeoUtils.convert_to_array("data/99/in_alt.dat") == alt_99
    end
  end
end



@testset "GeoUtils.jl" begin
  @testset "Code quality (Aqua.jl)" begin
    Aqua.test_ambiguities(GeoUtils)
    Aqua.test_unbound_args(GeoUtils)
    Aqua.test_undefined_exports(GeoUtils)
  end
  # Write your tests here.
end
