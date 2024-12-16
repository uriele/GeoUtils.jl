module GeoUtils
  using Reexport
  using Unitful
  using Unitful: ustrip,unit,°,m,km
  using CoordRefSystems
  using DataFrames
  using UnitfulData
  using Dates
  using IsacBinaryReader
  using IsacFileReader
  using StaticArrays
  using LinearAlgebra: dot,qr,Diagonal
  using CoordinateTransformations: LinearMap
  using Core.Intrinsics: sqrt_llvm
  using ScopedValues
  @reexport using Unitful: s,km,m,g,kg
  @reexport using Unitful:°C,K,°F # temperature units
  @reexport using Unitful:Pa,atm,bar # pressure units
  @reexport using Unitful:μm,nm,cm,m,km # length units
  import UnitfulData.Byte as byte
  import Unitful.° as deg
  import Unitful.μs as us
  using Base:IEEEFloat
  include("Utils.jl")
  include("ReadData.jl")
  include("RefractionIndex.jl")
  include("Orbit.jl")
  export get_data,convert_to_array,fix_latitudes
  export RealNumber,isRealNumber,isNotRealNumber
  export latitude,longitude,altitude
  export unique_altitudes,unique_latitudes,unique_longitudes
  export AirModel
  export Mathar,Mathar1,Mathar2,Mathar3,Mathar4
  export Ciddor
  export refractive_index
  export Vec2,Vec3
  export Ray2D,Ellipsoid
  export distance_from_unit_circle,distance_from_segment,distance_from_radius,distance_from_radius_new
  export h20_ppmv_to_rh
end
