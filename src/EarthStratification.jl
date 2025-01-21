abstract type AbstractLocalAtmosphere{Datum,T<:IEEEFloat} end
const KM{T}=Quantity{T,𝐋,typeof(km)}

"""

   LocalAtmosphere2D(;pressure::T=101325.0,temperature::T=15.0,height::T=0.0,latitude::T=0.0,x::T=0.0,y::T=0.0) where T<: IEEEFloat
Define the local atmosphere at a point in the Earth. The default values are the standard values at sea level.
"""
struct LocalAtmosphere2D{Datum,T}<:AbstractLocalAtmosphere{Datum, T}
  pressure::T
  temperature::T
  h::T
  θ::T
  w::T
  z::T
end

LocalAtmosphere2D{Datum}(pressure::T,temperature::T,h::T,θ::T,w::T,z::T) where {Datum,T} =
  LocalAtmosphere2D{Datum,T}(pressure,temperature,h,θ,w,z)
LocalAtmosphere2D(pressure,temperature,h,θ,w,z)=
  LocalAtmosphere2D{WGS84Latest}(pressure,temperature,h,θ,w,z)

function LocalAtmosphereLLA2D(datum::Type{Datum},pressure::T,temperature::T,h::T,θ::T) where {T,Datum}
  convert(ECEF2D{datum},LLA2D{datum}(h,θ)) |> x-> (x.w,x.z) |>
  x-> [ustrip(x) for x in x] |> x-> LocalAtmosphere2D{datum}(T(pressure),T(temperature),T(h),T(θ),T(x[1]),T(x[2]))
end
LocalAtmosphereLLA2D(pressure::T,temperature::T,h::T,θ::T) where T = LocalAtmosphereLLA2D(WGS84Latest,pressure,temperature,h,θ)

function LocalAtmosphereECEF2D(datum::Type{Datum},pressure::T,temperature::T,w::T,z::T) where {T,Datum}
  convert(LLA2D{datum},ECEF2D{datum}(w,z)) |> x-> (x.h,x.θ) |>
  x-> [ustrip(x) for x in x] |> x-> LocalAtmosphere2D{datum}(T(pressure),T(temperature),T(x[1]),T(x[2]),T(w),T(z))
end
LocalAtmosphereECEF2D(pressure::T,temperature::T,w::T,z::T) where T = LocalAtmosphereECEF2D(WGS84Latest,pressure,temperature,w,z)

LocalAtmosphere2D(pressure::T,temperature::T,h::T,θ::T) where T = LocalAtmosphereLLA2D(pressure,temperature,h,θ)


@inline function _normalize_ellipse!(h::A,majoraxis_earth::T) where A<:AbstractArray{T} where T
  h./=majoraxis_earth
  return nothing
end


function read_local_atmosphere(::Type{T},folder::String;
  skip_altitude=16,skip_latitude=15,skip_temperature=13,skip_pressure=13,
  units_pressure::Unitful.PressureUnits=u"mbar",
  units_temperature::Unitful.TemperatureUnits=u"K",
  units_altitude::Unitful.LengthUnits=u"km",
  units_latitude::Unitful.FreeUnits=u"°",
  datum::Datum=WGS84Latest,
  normalize::Bool=true
  )::StructArray where {T,Datum}
  @assert(units_latitude==u"°" || units_latitude==u"rad", "units_latitude must be either u\"°\" or u\"rad\"")

  # Read the data and convert to degree
  θ=convert_to_array(T,"$(folder)/in_lat.dat";skip=skip_latitude) |>
  x-> uconvert.(u"°",x.*units_latitude) |> x-> ustrip.(x)

  # Read the data and convert to km
  h=convert_to_array(T,"$(folder)/in_alt.dat";skip=skip_altitude) |>
  x-> uconvert.(u"km",x.*units_altitude) |> x-> ustrip.(x)

  # Read the data and convert to °C
  temperature=convert_to_array(T,"$(folder)/in_temp.dat";skip=skip_temperature) |>
  x-> uconvert.(u"°C",x.*units_temperature) |> x-> ustrip.(x) |>
  x-> reshape(x,(length(h),length(θ)))

  # Read the data and convert to Pa
  pressure=convert_to_array(T,"$(folder)/in_pres.dat";skip=skip_pressure) |>
  x-> uconvert.(u"Pa",x.*units_pressure) |> x-> ustrip.(x) |>
  x-> reshape(x,(length(h),length(θ)))

  majoraxis_earth=majoraxis(ellipsoid(datum)) |> x-> uconvert(km,x) |> ustrip
  squared_eccentricity_earth=eccentricity²(ellipsoid(datum))

  if normalize

    _normalize_ellipse!(h,majoraxis_earth)
    setNormalizedEarth(squared_eccentricity_earth)
    datum=NormalizedEarth
end

  StructArray(
  [
      LocalAtmosphereLLA2D(datum,
        pressure[i,j],
        temperature[i,j],
        h[i],
        θ[j]
        )
    for i in eachindex(h), j in eachindex(θ)
  ]
  )
end

read_local_atmosphere(folder::String;kwargs...)=read_local_atmosphere(Float64,folder;kwargs...)

@inline function _convertellipse2cartesian(θ::T,h::T,majoraxis_earth::T,squared_eccentricity_earth::T)::NTuple{2,T} where T
  N=majoraxis_earth/sqrt(1+squared_eccentricity_earth*sind(θ)*sind(θ))
  return (N+h)*cosd(θ),(N+h)*sind(θ)
end


using Interpolations
"""
  `discretize_atmosphere(atmosphere::A,levels,radii)`

Discretize the atmosphere into levels and radii.

# Input
- `atmosphere::A`: The local atmosphere
- `levels`: The number of levels (can be an integer for linear spacing or a vector for custom spacing)
- `radii`: The number of radii (can be an integer for linear spacing or a vector for custom spacing)

# Optional arguments
- `wavelength=10.0`: The wavelength of the light in μm
- `model=Ciddor()`: The model for the refractive index
- `interpolation_pressure=LinearPressure()`: The interpolation for the pressure

# Output
- `pressure`: The pressure in Pa
- `temperature`: The temperature in °C
- `refractive`: The refractive index


"""
function discretize_atmosphere(atmosphere::A,levels::Int,radii::Int;kwargs...)   where {A<:AbstractArray{Atm}} where Atm<:AbstractLocalAtmosphere{Datum,T} where {Datum,T}
  h_max=maximum(atmosphere.h)
  θ_max=max(maximum(atmosphere.θ),360)
  hᵢ=LinRange(h_max,0,levels)
  θᵢ=LinRange(0,θ_max,radii)
  discretize_atmosphere(atmosphere,hᵢ,θᵢ;kwargs...)
end


function discretize_atmosphere(atmosphere::A,hᵢ::AbstractVector{V},radii::Int;kwargs...) where {A<:AbstractArray{Atm},V} where Atm<:AbstractLocalAtmosphere{Datum,T} where {Datum,T}
  θ_max=max(maximum(atmosphere.θ),360)
  θᵢ=LinRange(0,θ_max,radii)
  discretize_atmosphere(atmosphere,hᵢ,θᵢ;kwargs...)
end

function discretize_atmosphere(atmosphere::A,levels::Int,θᵢ::AbstractVector{V};kwargs...) where {A<:AbstractArray{Atm},V} where Atm<:AbstractLocalAtmosphere{Datum,T} where {Datum,T}
  h_max=maximum(atmosphere.h)
  hᵢ=LinRange(h_max,0,levels)
  discretize_atmosphere(atmosphere,hᵢ,θᵢ;kwargs...)
end

discretize_atmosphere(atmosphere::A,hᵢ::AbstractVector{V},θᵢ::AbstractVector{V};kwargs...) where {A<:AbstractArray{Atm},V} where Atm<:AbstractLocalAtmosphere{Datum,T} where {Datum,T} =
_discretize_atmosphere_lowallocation(
  view(atmosphere.pressure,:,:),
  view(atmosphere.temperature,:,:),
  view(atmosphere.h,:,1),
  view(atmosphere.θ,1,:),
  view(hᵢ,:),view(θᵢ,:);kwargs...)


abstract type AbstractPressureInterpolation end
struct LogarithmicPressure<:AbstractPressureInterpolation end
struct LinearPressure<:AbstractPressureInterpolation end

@inline _interpolation_pt_h(::AbstractPressureInterpolation,knots,pressure,temperature,values,extrapolation_bc)=  throw(ArgumentError("Interpolation not implemented for $(typeof(interpolation))"))


@inline function _interpolate_pt_theta(knots,pressure,temperature,values,extrapolation_bc)
  p1=SemiCircularArray(linear_interpolation(knots,pressure,extrapolation_bc=extrapolation_bc)(values));
  t1=SemiCircularArray(linear_interpolation(knots,temperature,extrapolation_bc=extrapolation_bc)(values));
  return ((p1[1:end]+p1[2:end+1])./2,(t1[1:end]+t1[2:end+1])./2)
end


@inline function _interpolate_pt_theta!(
  pressure_interpolated::V,temperature_interpolated::V,
  knots::K,
  pressure_discrete::V1,temperature_discrete::V1,
  values::P,extrapolation_bc) where {V<:AbstractVector{T},V1<:AbstractVector{T},K,P} where T


  itp_p=linear_interpolation(knots,pressure_discrete,extrapolation_bc=extrapolation_bc);
  itp_t=linear_interpolation(knots,temperature_discrete,extrapolation_bc=extrapolation_bc);

  for i in eachindex(values[1:end-1])
    pressure_interpolated[i]=(itp_p(values[i])+itp_p(values[i+1]))*T(0.5)
    temperature_interpolated[i]=(itp_t(values[i])+itp_t(values[i+1]))*T(0.5)
  end
  pressure_interpolated[end]=(itp_p(values[end])+itp_p(values[1]))*T(0.5)
  temperature_interpolated[end]=(itp_t(values[end])+itp_t(values[1]))*T(0.5)
  return nothing
end





@inline function _interpolate_pt_h(::LinearPressure,knots,pressure,temperature,values,extrapolation_bc)
    pressure[pressure.<=0].=eps()  # Avoid log(0)
    p1=linear_interpolation(knots,pressure,extrapolation_bc=extrapolation_bc)(values);

    t1=linear_interpolation(knots,temperature,extrapolation_bc=extrapolation_bc)(values);
    return (
       p1[1:end-1].*exp.((log.(p1[2:end]).-log.(p1[1:end-1]))./2.0),
      (t1[1:end-1].+t1[2:end])./2
      )
end

@inline function _interpolate_pt_h(::LogarithmicPressure,knots,pressure,temperature,values,extrapolation_bc)
  pressure[pressure.<=0].=eps()  # Avoid log(0)
  logp1=linear_interpolation(knots,log.(pressure),extrapolation_bc=extrapolation_bc)(values);
  p1=exp.(logp1)
  t1=linear_interpolation(knots,temperature,extrapolation_bc=extrapolation_bc)(values);
  return (
    (p1[1:end-1]-p1[2:end])./(logp1[1:end-1]-logp1[2:end]),
    (t1[1:end-1]+t1[2:end])./2
    )
end




@inline function _interpolate_pt_h!(::LinearPressure,
  pressure_interpolated::V,temperature_interpolated::V,
  knots::K,
  pressure_discrete::V1,temperature_discrete::V1,
  values::P,extrapolation_bc) where {V<:AbstractVector{T},V1<:AbstractVector{T},K,P} where T

  itp_p=linear_interpolation(knots,pressure_discrete,extrapolation_bc=extrapolation_bc);
  itp_t=linear_interpolation(knots,temperature_discrete,extrapolation_bc=extrapolation_bc);

  for i in eachindex(values[1:end-1])
    p0=itp_p(values[i])
    p1=itp_p(values[i+1])
    pressure_interpolated[i]=p0*exp((log(p1)-log(p0))*T(0.5))
    temperature_interpolated[i]=(itp_t(values[i])+itp_t(values[i+1]))*T(0.5)
  end

  return nothing
end

@inline function _interpolate_pt_h!(::LogarithmicPressure,
  pressure_interpolated::V,temperature_interpolated::V,
  knots::K,
  pressure_discrete::V1,temperature_discrete::V1,
  values::P,extrapolation_bc) where {V<:AbstractVector{T},V1<:AbstractVector{T},K,P} where T

  itp_lp=linear_interpolation(knots,log.(pressure_discrete),extrapolation_bc=extrapolation_bc);
  itp_t=linear_interpolation(knots,temperature_discrete,extrapolation_bc=extrapolation_bc);

  for i in eachindex(values[1:end-1])
    logp1=itp_lp(values[i])
    logp2=itp_lp(values[i+1])
    pressure_interpolated[i]=(exp.(logp1)-exp.(logp2))/(logp1-logp2)
    temperature_interpolated[i]=0.5*(itp_t(values[i])+itp_t(values[i+1]))
  end

  return nothing
end



@inline function _discretize_atmosphere(atmosphere::A,hᵢ::AbstractVector{V},θᵢ::AbstractVector{V};wavelength=10.0,
  model=Ciddor(),
  interpolation_pressure=LinearPressure(),
  kwargs...) where {A<:AbstractArray{Atm},V} where Atm<:AbstractLocalAtmosphere{Datum,T} where {Datum,T}
  h_knots= atmosphere.h[:,1]
  θ_knots= atmosphere.θ[1,:]
  hᵢ=sort(hᵢ;rev=true)
  θᵢ=sort(θᵢ;rev=false)
  radii=length(θᵢ)
  levels=length(hᵢ)

  Press1= Array{T,2}(undef,size(atmosphere,1),radii)
  Temp1 = Array{T,2}(undef,size(atmosphere,1),radii)
  # Average Across Height (linear for temperature, exponential for pressure)
  Press_out= Matrix{T}(undef,levels-1,radii)
  Temp_out = similar(Press_out)
  Refr_out = similar(Press_out)

  @debug "Using $(typeof(interpolation_pressure)) interpolation for pressure along height"

  # Average Along Phi (circular interpolation)
  [
   (Press1[i,:],Temp1[i,:]) = _interpolate_pt_theta(θ_knots,atmosphere.pressure[i,:],atmosphere.temperature[i,:],θᵢ,Periodic())
   for i in axes(atmosphere,1)
  ]

    @debug "Using $(typeof(model)) model for refractive index"
    [
    begin
      (parent(Press_out)[:,j],parent(Temp_out)[:,j]) = _interpolate_pt_h(interpolation_pressure,reverse(h_knots),reverse(Press1[:,j]),reverse(Temp1[:,j]),hᵢ,Flat())
      [Refr_out[i,j]= refractive_index(model;wavelength=wavelength,temperature=Temp_out[i,j],pressure=Press_out[i,j],CO2ppm=0.0) for i in axes(Press_out,1)]
    end
      for j in axes(Press1,2)
    ]

  return SemiCircularMatrix(Press_out),
         SemiCircularMatrix(Temp_out),
         SemiCircularMatrix(Refr_out),
         SemiCircularMatrix([(h,θ) for h in hᵢ, θ in θᵢ]),
         SemiCircularMatrix([
    let
    convert(ECEF2D{Datum},LLA2D{Datum}(h,θ)) |> x-> (x.w,x.z) |>
    x-> (ustrip(x[1]),ustrip(x[2]))
    end
    for h in hᵢ, θ in θᵢ])
end


@inline function _discretize_atmosphere_lowallocation(
  pressure_knots::M,
  temperature_knots::M,
  h_knots::V1,θ_knots::V2,hᵢ::G,θᵢ::G;wavelength=10.0,
  model=Ciddor(),
  interpolation_pressure=LinearPressure(),
  kwargs...) where {M<:AbstractMatrix{T},V1<:AbstractVector{T},V2<:AbstractVector{T},G} where T

  idx_h=zeros(Int,length(h_knots))
  idx_θ=zeros(Int,length(θ_knots))
  sortperm!(idx_h,h_knots)
  sortperm!(idx_θ,θ_knots)
  h_knots.= view(h_knots,idx_h)
  θ_knots.= view(θ_knots,idx_θ)


  pressure_knots.=view(pressure_knots,idx_h,idx_θ)
  temperature_knots.=view(temperature_knots,idx_h,idx_θ)


  sort!(hᵢ;rev=true)
  sort!(θᵢ;rev=false)
  radii=length(θᵢ)
  levels=length(hᵢ)

  pressure_interpolated= SemiCircularMatrix(Matrix{T}(undef,levels-1,radii))
  temperature_interpolated =SemiCircularMatrix(similar(parent(pressure_interpolated)))
  refraction_interpolated = SemiCircularMatrix(similar(parent(pressure_interpolated)))

  pressure_interpolated_θ=Matrix{T}(undef,length(h_knots),radii)
  temperature_interpolated_θ=Matrix{T}(undef,length(h_knots),radii)
  #write over the radii to update the value
  for i in axes(pressure_interpolated_θ,1)
    @inbounds _interpolate_pt_theta!(
      view(pressure_interpolated_θ,i,:),
      view(temperature_interpolated_θ,i,:),
      view(θ_knots,:),
      view(pressure_knots,i,:),
      view(temperature_knots,i,:),
      view(θᵢ,:),
      Periodic())
  end

  # Average Across Height (linear for temperature, exponential for pressure)

    @debug "Using $(typeof(model)) model for refractive index"

  for j in axes(refraction_interpolated,2)
      @inbounds _interpolate_pt_h!(
        interpolation_pressure,
        view(pressure_interpolated,:,j),
        view(temperature_interpolated,:,j),
        view(h_knots,:),
        view(pressure_interpolated_θ,:,j),
        view(temperature_interpolated_θ,:,j),
        view(hᵢ,:),
        Flat())

      for i in axes(refraction_interpolated,1)
          @inbounds refraction_interpolated[i,j]= refractive_index(model;
          wavelength=wavelength,
          temperature=temperature_interpolated[i,j],
          pressure=pressure_interpolated[i,j],CO2ppm=0.0)
        end

    end

  return pressure_interpolated,temperature_interpolated,
        refraction_interpolated
end

struct SatOrbit{T<:IEEEFloat}
    z::T
    w::T
    ang::T
    h::T
    orbital_coordinate::T
    _normalized::Bool
    _islimb::Bool

    @inline SatOrbit{T}(z::T,w::T,ang::T,h::T,other::T,normalized::Bool=false,islimb::Bool=true) where T = new(z,w,ang,h,other,normalized,islimb)
end

"""
  SatOrbit(z,w,ang,h,orbital_coordinate)

Define an orbit for the satellite
"""
SatOrbit(z::T,w::T,ang::T,h::T,orbital_coordinate::T) where T<:IEEEFloat = SatOrbit{T}(z,w,ang,h,orbital_coordinate)
SatOrbit(z::I,w::I,ang::I,h::I,orbital_coordinate::I) where I<:Integer = SatOrbit(float.([z,w,ang,h,orbital_coordinate])...)
SatOrbit(z::Real,w::Real,ang::Real,h::Real,orbital_coordinate::Real) = SatOrbit(promote(z,w,ang,h,orbital_coordinate)...)
SatOrbit(z::KM,w::KM,ang,h::KM,orbital_coordinate) = SatOrbit(ustrip(z),ustrip(w),ang,ustrip(h),orbital_coordinate)
SatOrbit(z::Unitful.Length,w::Unitful.Length,ang,h::Unitful.Length,orbital_coordinate) =
SatOrbit(uconvert(km,z),uconvert(km,w),ang,uconvert(km,h),orbital_coordinate)

altitude(o::SatOrbit) = o.h
satellite_angle(o::SatOrbit) = o.ang
orbital_coordinate(o::SatOrbit) = o.orbital_coordinate
islimb(o::SatOrbit) = o._islimb
isnadir(o::SatOrbit) = !o._islimb



function __read_orbit(::Type{T},file::String) where T
  open(file) do f
    lines=readlines(f)

    pos=findall(x->!isnothing(match(r"^\sSEQUENCE",x)),lines);

    nseq=length(pos)
    mgeom=parse(Int,lines[pos[1]+6])
    orb=Matrix{SatOrbit{T}}(undef,nseq,mgeom)

    for j in 1:mgeom
      jj=7+j
      for i in eachindex(pos)
        sj=pos[i]+jj
        orb[i,j]=SatOrbit(parse.(T,split(lines[sj]))...)
      end
    end

    return orb
  end
end

__read_orbit(file::String)=__read_orbit(Float64,file)

function read_orbit(::Type{T},file::String) where T
  open(file) do f
    lines=readlines(f)
    (nseq,mgeom,seq)=findall(x->!isnothing(match(r"^\sSEQUENCE",x)),lines) |>
    x-> (x,length(x),parse.(Int,lines[x.+6]))|>
    x-> (x[2],x[3],
    let
      [(init_seq+2+6:init_seq+1+6+len_seq) for (init_seq,len_seq) in zip(x[1],x[3])]
    end)

    orb=Matrix{SatOrbit{T}}(undef,nseq,maximum(mgeom))


   [let
      [orb[i,j]=SatOrbit(parse.(T,split(lines[sj]))...) for (j,sj) in enumerate(seq[i])]
   end
   for i in eachindex(seq)
   ]
   return StructArray(orb)
  end
end

read_orbit(file::String)=read_orbit(Float64,file)


"""
  normalize_orbit!([datum=WGS84Latest],sorbit::StructArray{SatOrbit})

Normalize in-place the orbit with respect to the Earth's major axis.

Note: sorbit is a StructArray of SatOrbit (SatOrbit is an immutable struct).
"""
function normalize_orbit!(datum::Datum,orbit::Ao) where {Datum,Ao<:StructArray{O}} where O<:SatOrbit{T} where T
  majoraxis_earth= majoraxis(ellipsoid(datum)) |> ma-> uconvert(km,ma) |> ustrip
  for i in eachindex(orbit)
    orbit._normalized[i]==true && continue
    orbit.w[i]/=majoraxis_earth
    orbit.z[i]/=majoraxis_earth
    orbit.h[i]/=majoraxis_earth
    orbit._normalized[i]=true
  end
  return nothing
end
normalize_orbit!(orbit::Ao)  where {Ao<:StructArray{O}} where O<:SatOrbit{T} where T = normalize_orbit!(WGS84Latest,orbit)


"""
  normalize_orbit([datum=WGS84Latest],orbit::SatOrbit)

Normalize the orbit with respect to the Earth's major axis.
"""

function normalize_orbit(datum::Datum,orbit::O) where {Datum,O<:SatOrbit{T}} where T

  orbit._normalized==true && return orbit
  majoraxis_earth= majoraxis(ellipsoid(datum)) |> ma-> uconvert(km,ma) |> ustrip
  w_orbit= orbit.w/majoraxis_earth
  z_orbit= orbit.z/majoraxis_earth
  h_orbit= orbit.h/majoraxis_earth
  return SatOrbit{T}(z_orbit,w_orbit,orbit.ang,h_orbit,orbit.orbital_coordinate,true)
end
normalize_orbit(orbit::O) where O<:SatOrbit{T} where T = normalize_orbit(WGS84Latest,orbit)
