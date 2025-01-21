
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


"""
    MJD2000::Int

The Modified Julian Date for the year 2000 at 00:00:00 UTC, used by the ENVISAT mission.
"""
const MJD2000=51544 +2400000.5
const HOURS2DAY=1/24
const MINUTES2HOURS=1/60
const SECONDS2MINUTES=1/60
const MICROSECONDS2SECOND=1e-6

"""
  mjd2000_to_jd(days::Int=0,hours::Int=0,minutes::Int=0,seconds::Int=0,microseconds::Number=0.0)

Transform the input time given in days, hours, minutes, seconds and microseconds from the Modified Julian Date 2000 to the Julian Date.
"""
function mjd2000_to_jd(days::Int=0,hours::Int=0,minutes::Int=0,seconds::Int=0,microseconds::Int=0)::Float64
      float(MJD2000+days+(HOURS2DAY*(hours+MINUTES2HOURS*(minutes+SECONDS2MINUTES*(seconds+MICROSECONDS2SECOND*microseconds)))))
end


function ellipsfrome²(squared_eccentricity_earth)
  minoraxis_earth=(1-squared_eccentricity_earth)^(1/2)
  CoordRefSystems.ellipfromab(1,minoraxis_earth)
end

@inline function _NN(θrad::T,squared_eccentricity_earth::T)::T where T
    sinθ=sin(θrad)
    return 1/sqrt(1-squared_eccentricity_earth*sinθ*sinθ)
end

function geocentric_xy_to_geodesic_θ(datum::Datum,x::T,y::T)::Tuple{T,T} where {Datum,T}
  squared_eccentricity_earth=eccentricity²(ellipsoid(datum))
  majoraxis_earth=majoraxis(ellipsoid(datum))
  minoraxis_earth=minoraxis(ellipsoid(datum))
  minoraxis_earth/=majoraxis_earth
  x/=majoraxis_earth
  y/=majoraxis_earth
  θ=atan(y,x)



  @inline function _hh(x,y,θ)
    cosθ=cos(θ)
    if cosθ==0
      return abs(y)-minoraxis_earth
    end


    return x/cosθ-_NN(θ,squared_eccentricity_earth)
  end

  @inline function _θ(θ,x,z)
    R=_NN(θ)

    h=_hh(x,z,θ)


    return atan(z,x*(1-squared_eccentricity_earth*(R/(R+h))))
  end

  maxerror=1e-15

  θnew=_θ(θ,x,y)


  while abs(θnew-θ)>maxerror
    θ=θnew
    θnew=_θ(θ,x,y)
  end
  return (mod(rad2deg(θnew),360),_hh(x,y,θnew)*majoraxis_earth)
end

geocentric_xy_to_geodesic_θ(x,y)=geocentric_xy_to_geodesic_θ(NormalizedEarth,x,y)

#ref: https://www.oc.nps.edu/oc2902w/coord/coordcvt.pdf



# this finds the equivalent latitude on the surface of the earth h=0
# given the geocentric latitude θ
# unlike the angle at an altitude h>0, the result is exact and It
# does not require an iterative method
# tanθ_c= (1-e²)tanθ
@inline function _surface_from_geocentric_θdeg_to_geodesic_θdeg(datum::Datum, θc::T) where {Datum, T}
  θc==360 && return θc   # angle is the same at 0,90,180,270,360, but the code would return 0 for 360
  squared_eccentricity_earth = eccentricity²(ellipsoid(datum))
  y,x= sind(θc),cosd(θc)
  return atand(y, x*(1 - squared_eccentricity_earth)) |> x-> ifelse(x==360,x,mod(x,360))
end

# Function to convert geodesic angle to geocentric angle with quadrant information
@inline function _surface_from_geodesic_θdeg_to_geocentric_θdeg(datum::Datum, θ::T) where {Datum, T}
  θ==360 && return θ
  squared_eccentricity_earth = eccentricity²(ellipsoid(datum))
  y,x= sind(θ),cosd(θ)
  mod(atand(y*(1 - squared_eccentricity_earth), x),360)

end

@inline function _surface_from_geocentric_θrad_to_geodesic_θrad(datum::Datum, θc::T) where {Datum, T}
  geodesic_angle = _surface_from_geocentric_θdeg_to_geodesic_θdeg(datum, rad2deg(θc))
  return deg2rad(geodesic_angle)
end

@inline function _surface_from_geodesic_θrad_to_geocentric_θrad(datum::Datum, θ::T) where {Datum, T}
  geocentric_angle = _surface_from_geodesic_θdeg_to_geocentric_θdeg(datum, rad2deg(θ))
  return deg2rad(geocentric_angle)
end

_surface_from_geocentric_θdeg_to_geodesic_θdeg(θc::T) where T = _surface_from_geocentric_θdeg_to_geodesic_θdeg(NormalizedEarth, θc)
_surface_from_geocentric_θrad_to_geodesic_θrad(θc::T) where T = _surface_from_geocentric_θrad_to_geodesic_θrad(NormalizedEarth, θc)
_surface_from_geodesic_θdeg_to_geocentric_θdeg(θ::T) where T = _surface_from_geodesic_θdeg_to_geocentric_θdeg(NormalizedEarth, θ)
_surface_from_geodesic_θrad_to_geocentric_θrad(θ::T) where T = _surface_from_geodesic_θrad_to_geocentric_θrad(NormalizedEarth, θ)





"""
  convert_surface_angle_geocentric_to_geodesic_deg([::Datum=NormalizedEarth],θ::T)::T where where T<:IEEEFloat

Convert a geocentric angle in degrees to the geodesic angle in degrees.
"""
convert_surface_angle_geocentric_to_geodesic_deg(datum::Datum,θ::T) where {Datum,T} = _surface_from_geocentric_θdeg_to_geodesic_θdeg(datum, θ)
convert_surface_angle_geocentric_to_geodesic_deg(θ::T) where T = _surface_from_geocentric_θdeg_to_geodesic_θdeg(θ)
"""
  convert_surface_angle_geocentric_to_geodesic_rad([::Datum=NormalizedEarth],θ::T)::T where where T<:IEEEFloat

Convert a geocentric angle in radians to the geodesic angle in radians.
"""
convert_surface_angle_geocentric_to_geodesic_rad(datum::Datum,θ::T) where {Datum,T} = _surface_from_geocentric_θrad_to_geodesic_θrad(datum, θ)
convert_surface_angle_geocentric_to_geodesic_rad(θ::T) where T = _surface_from_geocentric_θrad_to_geodesic_θrad(θ)
"""
  convert_surface_angle_geodesic_to_geocentric_deg([::Datum=NormalizedEarth],θ::T)::T where where T<:IEEEFloat

Convert a geodesic angle in degrees to the geocentric angle in degrees.
"""
convert_surface_angle_geodesic_to_geocentric_deg(datum::Datum,θ::T) where {Datum,T} = _surface_from_geodesic_θdeg_to_geocentric_θdeg(datum, θ)
convert_surface_angle_geodesic_to_geocentric_deg(θ::T) where T = _surface_from_geodesic_θdeg_to_geocentric_θdeg(θ)
"""
  convert_surface_angle_geodesic_to_geocentric_rad([::Datum=NormalizedEarth],θ::T)::T where where T<:IEEEFloat

Convert a geodesic angle in radians to the geocentric angle in radians.
"""
convert_surface_angle_geodesic_to_geocentric_rad(datum::Datum,θ::T) where {Datum,T} = _surface_from_geodesic_θrad_to_geocentric_θrad(datum, θ)
convert_surface_angle_geodesic_to_geocentric_rad(θ::T) where T = _surface_from_geodesic_θrad_to_geocentric_θrad(θ)



"""
  convert_surface_angle_geocentric_to_geodesic_deg!([::Datum=NormalizedEarth],θ::A) where A<:AbstractArray{T} where T

Convert an array of geocentric angle in degree to the geodesic angle in degrees in-place.
"""
convert_surface_angle_geocentric_to_geodesic_deg!(θ::A) where A<:AbstractArray{T} where T = map!(_surface_from_geocentric_θdeg_to_geodesic_θdeg, θ, θ)
convert_surface_angle_geocentric_to_geodesic_deg!(datum::Datum,θ::A) where {Datum,A<:AbstractArray{T}} where T =
  map!(x-> _surface_from_geocentric_θdeg_to_geodesic_θdeg(datum,x), θ, θ)

"""
  convert_surface_angle_geocentric_to_geodesic_rad!([::Datum=NormalizedEarth],θ::A) where A<:AbstractArray{T} where T

Convert an array of geocentric angle in radians to the geodesic angle in radians in-place.
"""
convert_surface_angle_geocentric_to_geodesic_rad!(θ::A) where A<:AbstractArray{T} where T = map!(_surface_from_geocentric_θrad_to_geodesic_θrad, θ, θ)
convert_surface_angle_geocentric_to_geodesic_rad!(datum::Datum,θ::A) where {Datum,A<:AbstractArray{T}} where T =
  map!(x-> _surface_from_geocentric_θrad_to_geodesic_θrad(datum,x), θ, θ)
"""
  convert_surface_angle_geodesic_to_geocentric_deg!([::Datum=NormalizedEarth],θ::A) where A<:AbstractArray{T} where T

Convert an array of geodesic angle in degrees to the geocentric angle in degrees in-place.
"""
convert_surface_angle_geodesic_to_geocentric_deg!(θ::A) where A<:AbstractArray{T} where T = map!(_surface_from_geodesic_θdeg_to_geocentric_θdeg, θ, θ)
convert_surface_angle_geodesic_to_geocentric_deg!(datum::Datum,θ::A) where {Datum,A<:AbstractArray{T}} where T =
  map!(x-> _surface_from_geodesic_θdeg_to_geocentric_θdeg(datum,x), θ, θ)

"""
  convert_surface_angle_geodesic_to_geocentric_rad!([::Datum=NormalizedEarth],θ::A) where A<:AbstractArray{T} where T
Convert an array of geodesic angle in radians to the geocentric angle in radians in-place.
"""
convert_surface_angle_geodesic_to_geocentric_rad!(θ::A) where A<:AbstractArray{T} where T = map!(_surface_from_geodesic_θrad_to_geocentric_θrad, θ, θ)
convert_surface_angle_geodesic_to_geocentric_rad!(datum::Datum,θ::A) where {Datum,A<:AbstractArray{T}} where T =
  map!(x-> _surface_from_geodesic_θrad_to_geocentric_θrad(datum,x), θ, θ)
