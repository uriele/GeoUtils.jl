using StructArrays
using Interpolations
@kwdef struct LocalAtm{T}
  pressure::T=101325.0
  temperature::T=15.0
  height::T=0.0
  latitude::T=0.0
  x::T=0.0
  y::T=0.0
end


abstract type LatitudeType end

"""
  GeocentricLatitude <: LatitudeType

Define the latitude to be expressed in geocentric coordinates.
"""
struct GeocentricLatitude<:LatitudeType end

"""
  GeodesicLatitude <: LatitudeType

Define the latitude to be expressed in geodesic coordinates.
"""
struct GeodesicLatitude<:LatitudeType end

function read_local_atmosphere(::Type{T},folder::String;
  skip_altitude=16,skip_latitude=15,skip_temperature=13,skip_pressure=13,
  units_pressure::Unitful.PressureUnits=u"mbar",
  units_temperature::Unitful.TemperatureUnits=u"K",
  units_altitude::Unitful.LengthUnits=u"km",
  units_latitude::Unitful.FreeUnits=u"°",
  datum::Datum=WGS84Latest,
  normalize::Bool=true,
  latitude_type::LatitudeType=GeocentricLatitude()
  )::StructArray where {T,Datum}
  @assert(units_latitude==u"°" || units_latitude==u"rad", "units_latitude must be either u\"°\" or u\"rad\"")

  # Read the data and convert to degree
  θ=convert_to_array(T,"$(folder)/in_lat.dat";skip=skip_latitude) |>
  x-> uconvert.(u"°",x.*units_latitude) |> x-> ustrip.(x) |>
  x-> x.+90

  # Read the data and convert to km
  h=convert_to_array(T,"$(folder)/in_alt.dat";skip=skip_altitude) |>
  x-> uconvert.(u"km",x.*units_altitude) |> x-> ustrip.(x)

  # Read the data and convert to °C
  temperature=convert_to_array(T,"$(folder)/in_temp.dat";skip=skip_pressure) |>
  x-> uconvert.(u"°C",x.*units_temperature) |> x-> ustrip.(x) |>
  x-> reshape(x,(length(h),length(θ)))

  # Read the data and convert to Pa
  pressure=convert_to_array(T,"$(folder)/in_pres.dat";skip=skip_pressure) |>
  x-> uconvert.(u"Pa",x.*units_pressure) |> x-> ustrip.(x) |>
  x-> reshape(x,(length(h),length(θ)))

  a=majoraxis(ellipsoid(datum)) |> x-> uconvert(km,x) |> ustrip
  b=minoraxis(ellipsoid(datum)) |> x-> uconvert(km,x) |> ustrip

  if normalize
    (b,a)=_normalize_ellipse!(h,b,a)
  end

  StructArray(
  [
      begin val = convertellipse2cartesian(latitude_type,θ[j],h[i],a,b)
        LocalAtm(;
        pressure=pressure[i,j],
        temperature=temperature[i,j],
        height=h[i],
        latitude=θ[j],
        x=val[1],
        y=val[2]
        )
      end
    for i in eachindex(h), j in eachindex(θ)
  ]
  )
end

"""
  discretize_atmosphere(atmosphere::A,levels::Int,radii::Int)::B where {A<:AbstractArray{T},B<:AbstractArray{T}} where T<:LocalAtm
  discretize_atmosphere(atmosphere::A,levels::AbstractVector{V},radii::Int)::B where {A<:AbstractArray{T},B<:AbstractArray{T}} where T<:LocalAtm
discretize_atmosphere(atmosphere::A,levels::Int,radii::AbstractVector{V})::B where {A<:AbstractArray{T},B<:AbstractArray{T}} where T<:LocalAtm
discretize_atmosphere(atmosphere::A,levels::AbstractVector{V},radii::AbstractVector{V})::B where {Int<V<:Real A<:AbstractArray{T},B<:AbstractArray{T}} where T<:LocalAtm

Discretize the atmosphere into levels and radii.
"""
discretize_atmosphere(atmosphere)
"""
  convertellipse2cartesian([GeocentricLatitude::LatitudeType],θ::T,h::T,a::T,b::T)::NTuple{2,T} where T

Convert the coordinates from the geographic coordinates to the Cartesian coordinates. By default, the latitude
is assumed to be geocentric. The latitude can be geodesic by specifying the type of latitude GeodesicLatitude.

See: [`GeodesicLatitude`](@ref), [`GeocentricLatitude`](@ref).
"""
convertellipse2cartesian(::GeocentricLatitude,θ::T,h::T,a::T,b::T) where T = ((a+h)*cosd(θ),(b+h)*sind(θ))
function convertellipse2cartesian(::GeodesicLatitude,θ::T,h::T,a::T,b::T)::NTuple{2,T} where T
  a²=a*a
  b²=b*b
  cosx=cosd(θ)
  sinx=sind(θ)
  N=a²/sqrt(a²*cosx*cosx+b²*sinx*sinx)
  f=b²/a²
  return convertellipse2cartesian(GeocentricLatitude(),θ,h,N,f*N)
end
convertellipse2cartesian(θ::T,h::T,a::T,b::T) where {T} =convertellipse2cartesian(GeocentricLatitude(),θ,h,a,b)


read_local_atmosphere(folder::String;kwargs...)=read_local_atmosphere(Float64,folder;kwargs...)

atm=read_local_atmosphere("data_atm/INP_FILES";normalize=false)


scatter(atm.x,atm.y,markersize=1,legend=false)


@inline function _normalize_ellipse!(h::A,b::T,a::T)::NTuple{2,T} where A<:AbstractArray{T} where T
  h./=a
  return b./a,T(1)
end
