
"""
  files_are_same(file1::String, file2::String)::Bool

Check if two files are the same by comparing their contents and sizes.
"""
function files_are_same(file1::String, file2::String)::Bool
  # Check if both files exist
  if !isfile(file1) || !isfile(file2)
      return false
  end

  # Compare file sizes first for a quick check
  if stat(file1).size != stat(file2).size
      return false
  end

  # Read and compare file contents
  open(file1, "r") do f1
      open(file2, "r") do f2
          return read(f1) == read(f2)
      end
  end
end

"""
  all_files_are_same(directory::String, file::String="in_lat.dat")::Bool

Check if all files with the same name in subdirectories are the same by comparing their contents and sizes.
"""
function all_files_are_same(directory::String,file::String="in_lat.dat")::Bool
  cd(directory) do
    flag = true
    x=readdir()

    for i in 2:length(x)
      file1=x[i]*"/"*file
      file2=x[i-1]*"/"*file
      files_are_same(file1,file2) || return false
    end
    return true
  end
end


"""
  convert_to_array([::Type{T}=Float64],file::String;skip=15,sink::A=Array{T})::A where {A}

Convert a file to an array of Float64 values, skipping the first `skip` lines.

"""
function convert_to_array(Float64::Type{T},file::String;skip=15,sink::A=Array) where {T, A}
  open(file, "r") do f
    data=readlines(f)[skip+1:end] |>
    x-> join(x," ") |> x-> String.(split(x)) |> x-> parse.(T,x) |> sink
  end
end

convert_to_array(file::String;kwargs...)=convert_to_array(Float64,file;kwargs...)


abstract type RealNumber end
struct isRealNumber end
struct isNotRealNumber end


RealNumber(x::T) where T = isNotRealNumber()
RealNumber(::Type{<:Real}) = isRealNumber()
RealNumber(x::T) where T<:Real = RealNumber(typeof(x))
RealNumber(x::U) where U<:Unitful.Quantity = RealNumber(Unitful.ustrip(x))

isRealNumber(x) = isa(x,Number) && !isa(x,Complex)

"""
  fix_latitudes(lat::T,orbital_coordinates=true) where T<:Number
  fix_latitudes(lat::AbstractArray{T},orbital_coordinates=true) where T<:Number

Fix the latitudes so that the north pole is at the top of the raster.
"""
function fix_latitudes(lat::AbstractArray{T},orbital_coordinates=true) where T<:Number
  orbital_coordinates == false && return lat;
  # if the north pole is zero, we need to fix the latitudes
  # so that the north pole is at the top of the raster
  # also if the latitude is smaller than 180 => 90-lat
  # if the latitude is larger than 180 => lat-270
  GeoUtils._fix_latitudes(RealNumber(first(lat)),lat)
end

function fix_latitudes(lat::T,orbital_coordinates=true) where T<:Number
  orbital_coordinates == false && return lat;
  # if the north pole is zero, we need to fix the latitudes
  # so that the north pole is at the top of the raster
  # also if the latitude is smaller than 180 => 90-lat
  # if the latitude is larger than 180 => lat-270
  GeoUtils._fix_latitudes(RealNumber(lat),lat)
end

@inline _fix_latitudes(::isNotRealNumber,lat)= throw(ArgumentError("The latitude should be a real number"))

@inline function _fix_latitudes(::isRealNumber,lat::AbstractArray{T}) where T<:Real
  return [GeoUtils._conversion_latitude(_lat) for _lat in @. mod(lat,360)]
end

@inline function _fix_latitudes(::isRealNumber,lat::T) where T<:Real
  return GeoUtils._conversion_latitude(lat )
end

@inline function _fix_latitudes(::isRealNumber,lat::AbstractArray{U}) where U<:Unitful.Quantity
  _unit=Unitful.unit(first(lat));
  return GeoUtils._fix_latitudes(isRealNumber(),ustrip(lat)).*_unit
end

@inline function _fix_latitudes(::isRealNumber,lat::T) where T<:Unitful.Quantity
  _unit=Unitful.unit(lat);
  return GeoUtils._fix_latitudes(isRealNumber(),ustrip(lat)).*_unit
end

@inline _conversion_latitude(x::T) where T<:Real= x <= 180 ? (90 - x) : (x -270)

"""
  latitude(x::CoordRefSystems.Geographic)::Float64

Get the latitude from a `CoordRefSystems.Geographic` object.

See also [`longitude`](@ref), [`altitude`](@ref).
"""
latitude(x::CoordRefSystems.Geographic)= x.lat;
"""
  longitude(x::CoordRefSystems.Geographic)::Float64

Get the longitude from a `CoordRefSystems.Geographic` object.

See also [`latitude`](@ref), [`altitude`](@ref).
"""
longitude(x::CoordRefSystems.Geographic)=x.lon;

"""
  altitude(x::CoordRefSystems.Geographic)::Float64

Get the altitude from a `CoordRefSystems.Geographic` object that has an altitude.
By default only GeocentricLatLonAlt and LatLonAlt have an altitude.

See also [`latitude`](@ref), [`longitude`](@ref).
"""
altitude(x::GeocentricLatLonAlt)=x.alt;
altitude(x::LatLonAlt)=x.alt;


function parse_vrm_profile(directory::String,file::String)
  _files=[];
  cd(directory) do
    open(file, "r") do f
      text=readlines(f)
      regex_dividers=regex=r"^\s{2,4}[0-9]{1,3}";
      regex_title=r"[0-9]{1,3}\s([a-z,A-Z,0-9\+]*)\s?";
      dividers=findall([!isnothing(match.(regex,line)) for line in text])
      @debug dividers

      titles=[string(title.captures[1])  for title in match.(regex_title,text[dividers])]

      push!(dividers,length(text)+1)

      for (i,t) in enumerate(titles)
        push!(_files,"vrm_"*string(t)*".dat");
        if !isfile("vrm_"*string(t)*".dat")
          open("vrm_"*string(t)*".dat","w") do f
            for str in text[dividers[i]+1:dividers[i+1]-1]
              write(f,str)
            end
          end
        end
      end

    end
  end
  return _files

end


"""
  h20_ppmv_to_rh(ppmv::T; pressure::T=101325.0,temperature::T=15.0) where T<:IEEEFloat

Convert the water vapor concentration in parts per million by volume to relative humidity[^1].
Default values for pressure and temperature are 101325 Pa and 15°C respectively.

  RH=1000× Pd(H₂0)/Ps(H₂0)

where Pd(H₂0) is the partial pressure of water vapor and Ps(H₂0) is the saturation pressure
of water vapor in the mixture at the given temperature.

# inputs
- `ppmv::T`: water vapor concentration in parts per million by volume
- `pressure::T=101325.0`: pressure in Pascal
- `temperature::T=15.0`: temperature in °C
# outputs
- `RH::T`: relative humidity in %

[^1]: Huang, J. (2018). A simple accurate formula for calculating saturation vapor pressure of water and ice. Journal of Applied Meteorology and Climatology, 57(6), 1265-1272.
"""
function h20_ppmv_to_rh(ppmv::T; pressure::T=101325.0,temperature::T=15.0)::T where T<:IEEEFloat
pressure_h20=ppmv*pressure*PPMV_TO_PARTS

saturated_pressure_h20= (temperature>0)*_saturated_pressure_water(pressure,temperature)+
                        (temperature<=0)*_saturated_pressure_ice(pressure,temperature)
return min(100*pressure_h20/saturated_pressure_h20,100)
end
const PPMV_TO_PARTS=1e-6
const A_HUNG_WATER= 34.494
const B_HUNG_WATER= 4924.99
const C_HUNG_WATER= 237.1
const D_HUNG_WATER= 105.0
const E_HUNG_WATER= 1.5
const F_HUNG_WATER= 1.00071
const P_HUNG_WATER= 0.000000045

const A_HUNG_ICE= 43.497
const B_HUNG_ICE= 6545.8
const C_HUNG_ICE= 278.0
const D_HUNG_ICE= 868.0
const E_HUNG_ICE= 2.0
const F_HUNG_ICE= 0.99882
const P_HUNG_ICE= 0.00000008

for species in (:water,:ice)
  lower_species_str=String(species)
  upper_species_str=uppercase(lower_species_str)
  upper_species_symbol=Symbol(upper_species_str)
  fn = Symbol("_saturated_pressure_"*lower_species_str)
  hung_const_sym=Array{Symbol,1}(undef,7)
  for (i,hung_const) in enumerate(("A","B","C","D","E","F","P"))
    hung_const_sym[i]=Symbol(hung_const*"_HUNG_"*upper_species_str)
  end


  @eval @inline function $fn(pressure::T,temperature::T)::T where T<:IEEEFloat
    $species=exp($(hung_const_sym[1])-($(hung_const_sym[2])/(temperature+$(hung_const_sym[3]))))/(temperature+$(hung_const_sym[4]))^$(hung_const_sym[5])
    return $species*$(hung_const_sym[6])*exp($(hung_const_sym[7])*pressure)
  end
end
