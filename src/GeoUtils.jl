module GeoUtils
  using Reexport
  using Unitful
  using CoordRefSystems
  using DataFrames
  using UnitfulData
  using Dates
  using IsacBinaryReader
  using IsacFileReader
  using StaticArrays
  using Interpolations
  using LinearAlgebra: dot,qr,Diagonal
  using CoordinateTransformations: LinearMap
  using Core.Intrinsics: sqrt_llvm
  using ScopedValues
  using SatelliteToolboxTransformations
  using StructArrays
  using Accessors # @reset
  # Used for define and convert from LLA to ECEF and ECI
  using CoordRefSystems: Deg,Rad
  using CoordRefSystems: Geographic
  using CoordRefSystems: Datum
  using CoordRefSystems: raw,constructor,reconstruct,units
  using CoordRefSystems: fixlon
  using Unitful:Quantity,°
  using Unitful:𝐋
  using CoordRefSystems:ellipfromab
  using LinearAlgebra: ⋅
  using Polyester: @batch  # for batch processing
  #using Mmap: mmap
  ####
  #using GeoUtils
  #using CoordRefSystems
  #using StructArrays

  #using WGLMakie,Makie
  ####
  import Base.==
  import Base.convert
  import CoordRefSystems:ellipsoidparams,ellipsoid

  using CoordRefSystems
  using Random

  @reexport using Unitful: s,g,kg
  @reexport using Unitful:°C,K,°F # temperature units
  @reexport using Unitful:Pa,atm,bar # pressure units
  @reexport using Unitful:μm,nm,cm,m,km # length units
  @reexport using Unitful: uconvert
  import UnitfulData.Byte as byte
  import Unitful.° as deg
  import Unitful.μs as us
  import Base:IEEEFloat
  import Base.==
  import Unitful.Length as ULength
  @reexport using CoordRefSystems: majoraxis,minoraxis,ellipsoid,eccentricity²,flattening,eccentricity
  include("Utils.jl")
  const _NormalizedEarth🌎= Ref(ellipsfrome²(eccentricity²(CoordRefSystems.ellipsoid(WGS84Latest))))


  include("SemiCircularMatrix.jl")
  include("EarthCenteredInertial.jl")
  include("ReadData.jl")
  include("RefractionIndex.jl")
  include("EarthStratification.jl")
  include("SatelliteOrbit.jl")
  include("IOUtils.jl")
  export get_data,convert_to_array,fix_latitudes
  export RealNumber,isRealNumber,isNotRealNumber
  export latitude,longitude,altitude
  export unique_altitudes,unique_latitudes,unique_longitudes
  export AirModel
  export Mathar,Mathar1,Mathar2,Mathar3,Mathar4
  export Ciddor
  export refractive_index
  export Vec2,Vec3
  export Ray2D,NormalizedRay2D
  export Ellipsoid
  export distance_from_unit_circle,distance_from_segment
  export distance_from_radii
  export h20_ppmv_to_rh
  export SemiCircularMatrix,SemiCircularArray,SemiCircularVector
  export mjd2000_to_jd
  #export EarthCenteredInertial,ECI
  export EarthCenteredEarthFixed,ECEF
  export ECEF2D,LLA2D
  export  read_local_atmosphere, read_orbit,discretize_atmosphere
  export NormalizedEarth,ellipsfrome²
  export setNormalizedEarth,getNormalizedEarth
  export LocalAtmosphere2D,LocalAtmosphereECEF2D,LocalAtmosphereLLA2D
  export SatOrbit,normalize_datum,normalize_datum!
  export IntersectionStyle,NoIntersection,IsIntersection
  export LevelIntersection,RadiusIntersection,RadiusLevelIntersection,LevelRadiusIntersection
  export advance,bend,Interface
  export create_rays
  export Radius,get_direction,get_origin
  export origin,direction,isnormalized,islimb,isnadir
  export getIntersectionObjects

  export LeftIntersection,RightIntersection
  export TopIntersection,BottomIntersection

  export LeftTopIntersection,LeftBottomIntersection
  export RightTopIntersection,RightBottomIntersection

  export TopLeftIntersection,TopRightIntersection
  export BottomLeftIntersection,BottomRightIntersection

  export new_intersection

  export LogarithmicPressure,LinearPressure,AbstractPressureInterpolation
  export Carlotti,NoAtmosphere
  export geocentric_xy_to_geodesic_θ

  export convert_surface_angle_geocentric_to_geodesic_deg,convert_surface_angle_geocentric_to_geodesic_rad
  export convert_surface_angle_geodesic_to_geocentric_deg,convert_surface_angle_geodesic_to_geocentric_rad
  export convert_surface_angle_geocentric_to_geodesic_deg!,convert_surface_angle_geocentric_to_geodesic_rad!
  export convert_surface_angle_geodesic_to_geocentric_deg!,convert_surface_angle_geodesic_to_geocentric_rad!
  export create_radii_from_θ,scale_earth_by_h
  export setDebugIntersection,getDebugIntersection

  export create_bundle_rays,rays_from_orbit!
  #export AzimuthElevationRange,AER
  #export EastNorthUp,ENU
end
