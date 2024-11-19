module GeoUtils
  using Reexport
  using Unitful
  using Unitful: ustrip,unit,°,m,km
  using CoordRefSystems
  using DataFrames
  @reexport using Unitful:°C,K,°F # temperature units
  @reexport using Unitful:Pa,atm,bar # pressure units
  @reexport using Unitful:μm,nm,cm,m,km # length units
  using Base:IEEEFloat
  include("Utils.jl")
  include("ReadData.jl")
  include("RefractionIndex.jl")
  export get_data,convert_to_array,fix_latitudes
  export RealNumber,isRealNumber,isNotRealNumber
  export latitude,longitude,altitude
  export unique_altitudes,unique_latitudes,unique_longitudes
  export AirModel
  export Mathar,Mathar1,Mathar2,Mathar3,Mathar4
  export Ciddor
  export refractive_index

end
