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


@inline function _normalize_ellipse!(h::A,a::T) where A<:AbstractArray{T} where T
  h./=a
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

  a=majoraxis(ellipsoid(datum)) |> x-> uconvert(km,x) |> ustrip
  e²=eccentricity²(ellipsoid(datum))

  if normalize
    @debug "a: $a, h: $h"
    _normalize_ellipse!(h,a)
    setNormalizedEarth(e²)
    datum=NormalizeEarth
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

@inline function _convertellipse2cartesian(θ::T,h::T,a::T,e²::T)::NTuple{2,T} where T
  N=a/sqrt(1+e²*sind(θ)*sind(θ))
  return (N+h)*cosd(θ),(N+h)*sind(θ)
end


"""
  discretize_atmosphere(atmosphere::A,levels::Int,radii::Int)::B where {A<:AbstractArray{T},B<:AbstractArray{T}} where T<:AbstractLocalAtmosphere
  discretize_atmosphere(atmosphere::A,levels::AbstractVector{V},radii::Int)::B where {A<:AbstractArray{T},B<:AbstractArray{T}} where T<:AbstractLocalAtmosphere
discretize_atmosphere(atmosphere::A,levels::Int,radii::AbstractVector{V})::B where {A<:AbstractArray{T},B<:AbstractArray{T}} where T<:AbstractLocalAtmosphere
discretize_atmosphere(atmosphere::A,levels::AbstractVector{V},radii::AbstractVector{V})::B where {Int<V<:Real A<:AbstractArray{T},B<:AbstractArray{T}} where T<:AbstractLocalAtmosphere

Discretize the atmosphere into levels and radii.
"""

using Interpolations
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
  _discretize_atmosphere(atmosphere,hᵢ,θᵢ;kwargs...)


abstract type AbstractPressureInterpolation end
struct LogarithmicPressure<:AbstractPressureInterpolation end
struct LinearPressure<:AbstractPressureInterpolation end

@inline _interpolation_pt_h(::AbstractPressureInterpolation,knots,pressure,temperature,values,extrapolation_bc)=  throw(ArgumentError("Interpolation not implemented for $(typeof(interpolation))"))


@inline function _interpolate_pt_theta(knots,pressure,temperature,values,extrapolation_bc)
  p1=SemiCircularArray(linear_interpolation(knots,pressure,extrapolation_bc=extrapolation_bc)(values));
  t1=SemiCircularArray(linear_interpolation(knots,temperature,extrapolation_bc=extrapolation_bc)(values));
  return ((p1[1:end]+p1[2:end+1])./2,(t1[1:end]+t1[2:end+1])./2)
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

  @info "Using $(typeof(interpolation_pressure)) interpolation for pressure along height"

  # Average Along Phi (circular interpolation)
  [
   (Press1[i,:],Temp1[i,:]) = _interpolate_pt_theta(θ_knots,atmosphere.pressure[i,:],atmosphere.temperature[i,:],θᵢ,Periodic())
   for i in axes(atmosphere,1)
  ]
  # Average Across Height (linear for temperature, exponential for pressure)
    Press_out= Matrix{T}(undef,levels-1,radii)
    Temp_out = similar(Press_out)
    Refr_out = similar(Press_out)

    @info "Using $(typeof(model)) model for refractive index"
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

struct Orbit{T<:IEEEFloat}
    z::T
    w::T
    ang::T
    h::T
    other::T
end


function read_orbit(::Type{T},file::String) where T
  open(file) do f
    lines=readlines(f)
    (nseq,mgeom,seq)=findall(x->!isnothing(match(r"^\sSEQUENCE",x)),lines) |>
    x-> (x,length(x),parse.(Int,lines[x.+6]))|>
    x-> (x[2],x[3],
    let
      [(a+2+6:a+1+6+b) for (a,b) in zip(x[1],x[3])]
    end)

    orb=Matrix{Orbit{T}}(undef,nseq,maximum(mgeom))
    @debug typeof(orb)

   [let
      [orb[i,j]=Orbit(parse.(T,split(lines[sj]))...) for (j,sj) in enumerate(seq[i])]
   end
   for i in eachindex(seq)
   ]
   return StructArray(orb)
  end
end

read_orbit(file::String)=read_orbit(Float64,file)

normalize_orbit(datum::Datum,orbit::O) where {Datum,O<:Orbit{T}} where T = majoraxis(ellipsoid(datum)) |> a-> uconvert(km,a) |> ustrip|>
a-> Orbit(orbit.z/a,orbit.w/a,orbit.ang,orbit.h/a,orbit.other)
normalize_orbit(orbit::O) where O<:Orbit{T} where T = normalize_orbit(WGS84Latest,orbit)

struct ECEF2D{Datum,T<:IEEEFloat}
  w::T
  z::T
end

ECEF2D{Datum}(x::T,y::T) where {Datum,T<:IEEEFloat} =ECEF2D{Datum,T}(x,y)
ECEF2D{Datum}(x::Number,y::Number) where Datum = ECEF2D{Datum,float{Number}}(promote(x,y)...)
ECEF2D{Datum}(x::L1,y::L2) where {L1<:ULength,L2<:ULength,Datum} = ECEF2D{Datum}(ustrip(uconvert(km,x)),ustrip(uconvert(km,y)))
ECEF2D{Datum}(x::L,y::Number) where {L<:ULength,Datum} = ECEF2D{Datum}(ustrip(uconvert(km,x)),y)
ECEF2D{Datum}(x::Number,y::L) where {L<:ULength,Datum} = ECEF2D{Datum}(x,ustrip(uconvert(km,y)))

ECEF2D(x,y) = ECEF2D{WGS84Latest}(x,y)


Base.convert(::Type{ECEF2D{T,Datum}},coords::ECEF2D{Datum}) where {T,Datum} = ECEF2D{T,Datum}(coords.w,coords.z)
CoordRefSystems.constructor(::Type{<:ECEF2D{Datum}}) where {Datum} = ECEF2D{Datum}
==(coords₁::ECEF2D{Datum},coord₂::ECEF2D{Datum}) where {Datum}=
coords₁.z==coord₂.z && coords₁.w==coord₂.w

struct LLA2D{Datum,T<:IEEEFloat}
  h::T
  θ::T
end


LLA2D{Datum}(h::T,θ::T) where {Datum,T<:IEEEFloat} =LLA2D{Datum,T}(h,mod(θ,360))
LLA2D{Datum}(h::Number,θ::Number) where Datum = LLA2D{Datum,float(Number)}(promote(h,mod(θ,360))...)
LLA2D{Datum}(h::L,θ::Number) where {L<:ULength,Datum} = LLA2D{Datum}(ustrip(uconvert(km,h)),θ)
LLA2D{Datum}(h::Number,θ::D) where {D<:Deg,Datum} = LLA2D{Datum}(h,ustrip(θ))
LLA2D{Datum}(h::L,θ::D) where {L<:ULength,D<:Deg,Datum} = LLA2D{Datum}(ustrip(uconvert(km,h)),ustrip(θ))
LLA2D{Datum}(h::Number,θ::R) where {R<:Rad,Datum} = LLA2D{Datum}(h,ustrip(uconvert(°,θ)))
LLA2D{Datum}(h::L,θ::R) where {L<:ULength,R<:Rad,Datum} = LLA2D{Datum}(ustrip(uconvert(km,h)),ustrip(uconvert(°,θ)))

LLA2D(h::T,θ::T) where T<:IEEEFloat = LLA2D{WGS84Latest,T}(h,θ)


Base.convert(::Type{LLA2D{T,Datum}},coords::LLA2D{Datum}) where {T,Datum} = LLA2D{T,Datum}(coords.h,coords.θ)
CoordRefSystems.constructor(::Type{<:LLA2D{Datum}}) where {Datum} = LLA2D{Datum}
==(coords₁::LLA2D{Datum},coord₂::LLA2D{Datum}) where {Datum}=
coords₁.h==coord₂.h && coords₁.θ==coord₂.θ

CoordRefSystems.ellipsoid(datum::Datum) = ellipsoid(datum)


function getNormalizedEarth()
  NE=ellipsoid(NormalizeEarth)
  @debug "a: $(majoraxis(NE)), b: $(minoraxis(NE))"
  @debug "flattening: $(flattening(NE)), eccentricity: $(eccentricity(NE))"
  @debug "eccentricity²: $(eccentricity²(NE))"
  return _NormalizedEarth🌎[]
end


function setNormalizedEarth(e²::T) where T<:IEEEFloat
  @assert(0<=(e²)<=1, "The eccentricity squared must be between 0 and 1")
  _NormalizedEarth🌎[]=ellipsfrome²(e²)
  @debug "Normalized Earth set to e²=$e²"
  getNormalizedEarth()
  return nothing
end
setNormalizedEarth() = setNormalizedEarth(eccentricity²(CoordRefSystems.ellipsoid(WGS84Latest)))


abstract type NormalizeEarth🌎 <: RevolutionEllipsoid end
struct NormalizeEarth<:Datum end
ellipsoidparams(::Type{NormalizeEarth🌎}) = _NormalizedEarth🌎[]
ellipsoid(::Type{NormalizeEarth}) = NormalizeEarth🌎
eccentricity(ellipsoid(NormalizeEarth))
##############
# CONVERSIONS
##############

function Base.convert(::Type{ECEF2D{Datum}},coords::LLA2D{Datum,T}) where {T,Datum}
  🌎 = ellipsoid(Datum)
  a = majoraxis(🌎) |> ustrip
  e²= eccentricity²(🌎)
  h=ustrip(coords.h)
  θ=ustrip(coords.θ)
  sinθ=sind(θ)
  cosθ=cosd(θ)
  N=a/sqrt(1-e²*sinθ*sinθ)
  w=(N+h)*cosθ
  z=(N*(1-e²)+h)*sinθ
  return ECEF2D{Datum}(w,z)
end

function Base.convert(::Type{LLA2D{Datum}},coords::ECEF2D{Datum,T}) where {T,Datum}
  🌎 = ellipsoid(Datum)
  a = T(ustrip(majoraxis(🌎)))
  b = T(ustrip(minoraxis(🌎)))
  e² = T(eccentricity²(🌎))
  z = ustrip(coords.z)
  p = ustrip(coords.w)
  e′² = e² / (1 - e²)
  ψ = atand(a * z, b * p)
  ϕ = mod1(atand(z + b * e′² * sind(ψ)^3, p - a * e² * cosd(ψ)^3),360)
  #ϕ = mod1(atand(z + b * e′² * sind(ϕ)^3, p - a * e² * cosd(ϕ)^3),360)
  #ϕ = mod1(atand(z + b * e′² * sind(ϕ)^3, p - a * e² * cosd(ϕ)^3),360)

  N = a / sqrt(1 - e² * sind(ϕ)^2)
  h = p / cosd(ϕ) - N
  return LLA2D{Datum}(h,ϕ)
end
