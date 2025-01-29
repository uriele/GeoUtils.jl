#module CoordRefSystemsExt

  fixlat(lat)=ifelse(abs(lat)==90°,lat,rem(lat,180°,RoundNearest))
  fixaltitude2D(alt)=ifelse(abs(alt)==0deg,alt,mod(alt,360°))

  # Define the KM unit for the radius that is more common than the m used in the original package for satellite orbits
  const KM{T}=Quantity{T,u"𝐋",typeof(km)}

"""
  EarthCenteredEarthFixed(altitude,azimuth,radius)
  EarthCenteredEarthFixed{Datum}(altitude,azimuth,radius)

Earth Centered Earth Fixed given as E,F,G in length units (default to km) for a given `Datum` (default to WGS84Latest).

## Examples

```julia
EarthCenteredEarthFixed(3727.195946164105,2476.2243282336854,-4561.757835246039) # add default units and Datum
EarthCenteredEarthFixed(3727km,2476km,-4561km) # integer are converted to float
EarthCenteredEarthFixed{WGS84Latest}(3727.195946164105km,2476.2243282336854km,-4561.757835246039km) # specifies the Datum
````
"""
struct EarthCenteredEarthFixed{Datum,L<:Length} <: Geographic{Datum}
  E::L
  F::L
  G::L
end

EarthCenteredEarthFixed{Datum}(E::L,F::L,G::L) where {Datum,L<:Length}=EarthCenteredEarthFixed{Datum,float(L)}(E,F,G)

EarthCenteredEarthFixed{Datum}(E::Number,F::Number,G::L) where {Datum,L<:Length}=EarthCenteredEarthFixed{Datum}(L(E),L(F),G)
EarthCenteredEarthFixed{Datum}(E::Number,F::L,G::Number) where {Datum,L<:Length}=EarthCenteredEarthFixed{Datum}(L(E),F,L(G))
EarthCenteredEarthFixed{Datum}(E::L,F::Number,G::Number) where {Datum,L<:Length}=EarthCenteredEarthFixed{Datum}(E,L(F),L(G))
EarthCenteredEarthFixed{Datum}(E::Number,F::L,G::L) where {Datum,L<:Length}=EarthCenteredEarthFixed{Datum}(L(E),F,G)
EarthCenteredEarthFixed{Datum}(E::L,F::Number,G::L) where {Datum,L<:Length}=EarthCenteredEarthFixed{Datum}(E,L(F),G)
EarthCenteredEarthFixed{Datum}(E::L,F::L,G::Number) where {Datum,L<:Length}=EarthCenteredEarthFixed{Datum}(E,F,L(G))
EarthCenteredEarthFixed{Datum}(E::Number,F::Number,G::Number) where {Datum}=EarthCenteredEarthFixed{Datum}(KM(E),KM(F),KM(G))
EarthCenteredEarthFixed(args...)=EarthCenteredEarthFixed{WGS84Latest}(args...)


"""
    ECEF(altitude,azimuth,radius)
    ECEF{Datum}(altitude,azimuth,radius)

Alias to [`EarthCenteredEarthFixed`](@ref).

## Examples

```julia
ECEF(3727.195946164105,2476.2243282336854,-4561.757835246039) # add default units and Datum
ECEF(3727km,2476km,-4561km) # integer are converted to float
ECEF{WGS84Latest}(3727.195946164105km,2476.2243282336854km,-4561.757835246039km) # specifies the Datum
````
"""
const ECEF{Datum}=EarthCenteredEarthFixed{Datum}

lentype(::Type{<:EarthCenteredEarthFixed{Datum,L}}) where {Datum,L}=L

Base.convert(::Type{EarthCenteredEarthFixed{Datum,L}},coord::EarthCenteredEarthFixed{Datum}) where {Datum,L}=EarthCenteredEarthFixed{Datum,L}(coord.E,coord.F,coord.G)

CoordRefSystems.constructor(::Type{<:EarthCenteredEarthFixed{Datum}}) where {Datum} = EarthCenteredEarthFixed{Datum}

CoordRefSystems.reconstruct(C::Type{EarthCenteredEarthFixed{Datum}}, raw) where {Datum} = begin
  (E,F,G) = raw.*units(C)
  constructor(C)(E,F,G)
end

==(coords₁::ECEF{Datum},coord₂::ECEF{Datum}) where {Datum}=coords₁.E==coord₂.E && coords₁.F==coord₂.F && coords₁.G==coord₂.G

Random.rand(rng::Random.AbstractRNG,::Type{EarthCenteredEarthFixed{Datum}}) where {Datum}=let
  majoraxis_earth=majoraxis(ellipsoid(Datum)) |> x-> uconvert(km,x) |> ustrip
  squared_eccentricity_earth=eccentricity²(ellipsoid(Datum))
  h=rand(rng)*100;
  ϕ=-90+180*rand(rng)
  γ=-180+360*rand(rng)
  N=majoraxis_earth/sqrt(1-squared_eccentricity_earth*sin(ϕ)^2)
  EarthCenteredEarthFixed{Datum}(
    (N+h)*cosd(ϕ)*cosd(γ),
    (N+h)*cosd(ϕ)*sind(γ),
    (N*(1-squared_eccentricity_earth)+h)*sind(ϕ)
  )
end

Random.rand(rng::Random.AbstractRNG,::Type{EarthCenteredEarthFixed})= rand(rng,ECEF{WGS84Latest})

"""
    EarthCenteredInertial([::Type{ReferenceFrame}=MJD2000],X,Y,Z,time)
    EarthCenteredInertial{Datum}([::Type{ReferenceFrame}=MJD2000],X,Y,Z,time)

Earth Centered Inertial given as X,Y,Z in length units (default to km) for a given `Datum` (default to WGS84Latest),
and time in days from the epoch (default to J2000).


  ## Examples
  ```julia
  EarthCenteredInertial(45,0,6371) # add default units and Datum
  EarthCenteredInertial(45°,0°,6371km) # integer are converted to float
  EarthCenteredInertial((π/4)rad,0,6371.0km) # converts altitude to degrees, azimuth assumed to be in degrees
  EarthCenteredInertial{WGS84Latest}(45,0,6371) # specifies the Datum
  ```
  """
  struct EarthCenteredInertial{Datum,L<:KM,T} <: Geographic{Datum}
    X ::L
    Y ::L
    Z ::L
    time::T
  end
  EarthCenteredInertial{Datum}(altitude::D,azimuth::D,radius::L) where {Datum,D<:Deg,L<:Length}=
    EarthCenteredInertial{Datum,float(D),float(L)}(fixlat(altitude),fixlon(azimuth),radius)
  EarthCenteredInertial{Datum}(altitude::D,azimuth::D,radius::Number) where {Datum,D<:Deg}=
    EarthCenteredInertial{Datum}(altitude,azimuth,KM(radius))
  EarthCenteredInertial{Datum}(altitude::Deg,azimuth::Deg,radius) where {Datum}=
    EarthCenteredInertial{Datum}(promote(altitude,azimuth)...,radius)
  EarthCenteredInertial{Datum}(altitude::Rad,azimuth::Rad,radius) where {Datum}=
    EarthCenteredInertial{Datum}(rad2deg(altitude),rad2deg(azimuth),radius)
  EarthCenteredInertial{Datum}(altitude::Deg,azimuth::Rad,radius) where {Datum}=
    EarthCenteredInertial{Datum}(altitude,rad2deg(azimuth),radius)
  EarthCenteredInertial{Datum}(altitude::Rad,azimuth::Deg,radius) where {Datum}=
    EarthCenteredInertial{Datum}(rad2deg(altitude),azimuth,radius)

  EarthCenteredInertial{Datum}(altitude::Number,azimuth::Number,radius) where {Datum}=
    EarthCenteredInertial{Datum}(Deg(altitude),Deg(azimuth),radius)

  EarthCenteredInertial{Datum}(altitude::Number,azimuth::Deg,radius) where {Datum}=
    EarthCenteredInertial{Datum}(Deg(altitude),azimuth,radius)

  EarthCenteredInertial{Datum}(altitude::Deg,azimuth::Number,radius) where {Datum}=
    EarthCenteredInertial{Datum}(altitude,Deg(azimuth),radius)

    EarthCenteredInertial{Datum}(altitude::Number,azimuth::Rad,radius) where {Datum}=
    EarthCenteredInertial{Datum}(Deg(altitude),azimuth,radius)

  EarthCenteredInertial{Datum}(altitude::Rad,azimuth::Number,radius) where {Datum}=
    EarthCenteredInertial{Datum}(altitude,Deg(azimuth),radius)



  EarthCenteredInertial(altitude,azimuth,radius)=
    EarthCenteredInertial{WGS84Latest}(altitude,azimuth,radius)

"""
    ECI(altitude,azimuth,radius)
    ECI{Datum}(altitude,azimuth,radius)

Alias to [`EarthCenteredInertial`](@ref).


## Examples
```julia
ECI(45,0,6371) # add default units and Datum
ECI(45°,0°,6371km) # integer are converted to float
ECI((π/4)rad,0,6371.0km) # converts altitude to degrees, azimuth assumed to be in degrees
ECI{WGS84Latest}(45,0,6371) # specifies the Datum
```
"""
const ECI{Datum}=EarthCenteredInertial{Datum}


Base.convert(::Type{EarthCenteredInertial{Datum,D,L}},coord::EarthCenteredInertial{Datum}) where {Datum,D,L}=
EarthCenteredInertial{Datum,D,L}(coord.altitude,coord.azimuth,coord.radius)

CoordRefSystems.constructor(::Type{<:EarthCenteredInertial{Datum}}) where {Datum} = EarthCenteredInertial{Datum}
CoordRefSystems.reconstruct(C::Type{EarthCenteredInertial{Datum}}, raw) where {Datum} = begin
  (altitude,azimuth,radius) = raw.*units(C)
  constructor(C)(altitude,azimuth,radius)
end

lentype(::Type{<:EarthCenteredInertial{Datum,D,L}}) where {Datum,D,L}=L

==(coords₁::ECI{Datum},coord₂::ECI{Datum}) where {Datum}=
coords₁.altitude==coord₂.altitude && coords₁.azimuth==coord₂.azimuth && coords₁.radius==coord₂.radius

Random.rand(rng::Random.AbstractRNG,::Type{EarthCenteredInertial{Datum}}) where {Datum}=let
  majoraxis_earth=majoraxis(ellipsoid(Datum)) |> x-> uconvert(km,x) |> ustrip
  EarthCenteredInertial{Datum}(
    -90+180*rand(rng),
    -180+360*rand(rng),
    majoraxis_earth+800*abs(rand(rng))
  )
end
Random.rand(rng::Random.AbstractRNG,::Type{EarthCenteredInertial})= rand(rng,ECI{WGS84Latest})


"""
    EarthCenteredInertial2D(altitude,radius)
    EarthCenteredInertial2D{Datum}(altitude,radius)

  Earth Centered Inertial latitude `lat ∈ [0,360]` in angular units (default to degrees),
  and radius `r`  from the center of the Earth in length units (default to km) for a given `Datum` (default to WGS84Latest).

  ## Examples

  ```julia
  EarthCenteredInertial2D(45,6371) # add default units and Datum
  EarthCenteredInertial2D(45°,6371km) # integer are converted to float
  EarthCenteredInertial2D((π/4)rad,6371.0km) # converts altitude to degrees
  EarthCenteredInertial2D{WGS84Latest}(45,6371) # specifies the Datum
  ```
  """
  struct EarthCenteredInertial2D{Datum,D<:Deg,L<:KM} <: Geographic{Datum}
    altitude::D
    radius ::L
  end

  EarthCenteredInertial2D{Datum}(altitude::D,radius::L) where {Datum,D<:Deg,L<:Length}=
    EarthCenteredInertial2D{Datum,float(D),float(L)}(fixaltitude2D(altitude),radius)
  EarthCenteredInertial2D{Datum}(altitude::D,radius::Number) where {Datum,D<:Deg}=
    EarthCenteredInertial2D{Datum}(altitude,KM(radius))
  EarthCenteredInertial2D{Datum}(altitude::Rad,radius) where {Datum}=
    EarthCenteredInertial2D{Datum}(deg2rad(altitude),radius)
  EarthCenteredInertial2D{Datum}(altitude::Number,radius) where {Datum}=
    EarthCenteredInertial2D{Datum}(Deg(altitude),radius)

  EarthCenteredInertial2D(altitude,radius)=
    EarthCenteredInertial2D{WGS84Latest}(altitude,radius)



"""
    ECI2D(altitude,radius)
    ECI2D{Datum}(altitude,radius)

Alias to [`EarthCenteredInertial2D`](@ref).

## Examples

```julia
ECI2D(45,6371) # add default units and Datum
ECI2D(45°,6371km) # integer are converted to float
ECI2D((π/4)rad,6371.0km) # converts altitude to degrees
ECI2D{WGS84Latest}(45,6371) # specifies the Datum
```
"""
const ECI2D{Datum}=EarthCenteredInertial2D{Datum}

#=
Base.convert(::Type{EarthCenteredInertial2D{Datum,D,L}},coord::EarthCenteredInertial2D{Datum}) where {Datum,D,L}=
EarthCenteredInertial2D{Datum,D,L}(coord.altitude,coord.radius)

CoordRefSystems.constructor(::Type{<:EarthCenteredInertial2D{Datum}}) where {Datum} = EarthCenteredInertial2D{Datum}

CoordRefSystems.reconstruct(C::Type{EarthCenteredInertial2D{Datum}}, raw) where {Datum} = begin
  (altitude,radius) = raw.*units(C)
  constructor(C)(altitude,radius)
end

lentype(::Type{<:EarthCenteredInertial2D{Datum,D,L}}) where {Datum,D,L}=L

==(coords₁::ECI2D{Datum},coord₂::ECI2D{Datum}) where {Datum}=
coords₁.altitude==coord₂.altitude && coords₁.radius==coord₂.radius

Random.rand(rng::Random.AbstractRNG,::Type{EarthCenteredInertial2D{Datum}}) where {Datum}=let
  majoraxis_earth=majoraxis(ellipsoid(Datum)) |> x-> uconvert(km,x) |> ustrip
  EarthCenteredInertial2D{Datum}(
    360*rand(rng),
    majoraxis_earth+800*abs(rand(rng))
  )
end
Random.rand(rng::Random.AbstractRNG,::Type{EarthCenteredInertial2D})= rand(rng,ECI2D{WGS84Latest})
=#
#-------------
# CONVERSIONS
#-------------

# Adapted from CelesTrack
# reference: https://celestrak.org/columns/v02n01/
#            https://celestrak.com/columns/v02n02/
#            https://celestrak.com/columns/v02n03/


#
#  (1-squared_eccentricity_earth) * tan(ϕ) = (1-squared_eccentricity_earth) * tan(ϕ′)
#
#


function Base.convert(::Type{EarthCenteredEarthFixed{Datum}},coords::LatLonAlt{Datum}) where {Datum}
  🌎  = ellipsoid(Datum)
  ϕ′  = ustrip(coords.lat)
  squared_eccentricity_earth  = oftype(ϕ′, eccentricity²(🌎))
  majoraxis_earth   = majoraxis(🌎) |> x-> uconvert(km,x) |> ustrip
  N   = majoraxis_earth / sqrt(1 - squared_eccentricity_earth * sind(ϕ′)^2)

  h   = coords.alt |> x-> uconvert(km,x) |> ustrip

  EarthCenteredEarthFixed{Datum}(
    (N+h)*cosd(ϕ′)*cosd(coords.lon),
    (N+h)*cosd(ϕ′)*sind(coords.lon),
    (N*(1-squared_eccentricity_earth)+h)*sind(ϕ′))
end

# Use Zhu Algorithm Closed Form
#
function Base.convert(::Type{LatLonAlt{Datum}},coords::EarthCenteredEarthFixed{Datum}) where {Datum}
  🌎  = ellipsoid(Datum)

  X = coords.E |> x-> uconvert(km,x) |> ustrip
  Y = coords.F |> x-> uconvert(km,x) |> ustrip
  Z = coords.G |> x-> uconvert(km,x) |> ustrip
  λ = atand(Y/X)

  squared_eccentricity_earth  = oftype(λ, eccentricity²(🌎))
  majoraxis_earth   = majoraxis(🌎) |> x-> uconvert(km,x) |> ustrip
  P   = hypot(X,Y)

  # initial value
  ϕ′ = atand(P/Z)



  @inline function N(ϕ′)
    majoraxis_earth / sqrt(1-squared_eccentricity_earth*sind(ϕ′)^2)
  end

  N′  = N(ϕ′)

  while true
    ϕ = atand(Z/P/(1-squared_eccentricity_earth*N′*cosd(ϕ′)/P))
    err=abs(ϕ-ϕ′)
    #
    ϕ′ = ϕ
    N′= N(ϕ′)
    if err < eps()
      break
    end
  end
  h′ = P/cosd(ϕ′)-N′

  LatLonAlt{Datum}(Deg(ϕ′),Deg(λ),KM(h′))

end



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
  NE=ellipsoid(NormalizedEarth)



  return _NormalizedEarth🌎[]
end


function setNormalizedEarth(squared_eccentricity_earth::T) where T<:IEEEFloat
  @assert(0<=(squared_eccentricity_earth)<=1, "The eccentricity squared must be between 0 and 1")
  _NormalizedEarth🌎[]=ellipsfrome²(squared_eccentricity_earth)

  getNormalizedEarth()
  return nothing
end
setNormalizedEarth() = setNormalizedEarth(eccentricity²(CoordRefSystems.ellipsoid(WGS84Latest)))


abstract type NormalizedEarth🌎 <: RevolutionEllipsoid end
struct NormalizedEarth<:Datum end
ellipsoidparams(::Type{NormalizedEarth🌎}) = _NormalizedEarth🌎[]
ellipsoid(::Type{NormalizedEarth}) = NormalizedEarth🌎
eccentricity(ellipsoid(NormalizedEarth))
##############
# CONVERSIONS
##############

@inline function _from_lla2d_to_ecef2d(::Type{Datum},θ::T,h::T) where {T,Datum}
  🌎 = ellipsoid(Datum)
  majoraxis_earth = majoraxis(🌎) |> ustrip
  squared_eccentricity_earth= eccentricity²(🌎)
  sinθ=sind(θ)
  cosθ=cosd(θ)
  N=majoraxis_earth/sqrt(1-squared_eccentricity_earth*sinθ*sinθ)
  w=(N+h)*cosθ
  z=(N*(1-squared_eccentricity_earth)+h)*sinθ
  return (w,z)
end

@inline function _from_ecef2d_to_lla2d(::Type{Datum},p::T,z::T) where {T,Datum}
  🌎 = ellipsoid(Datum)
  majoraxis_earth = T(ustrip(majoraxis(🌎)))
  minoraxis_earth = T(ustrip(minoraxis(🌎)))
  squared_eccentricity_earth = T(eccentricity²(🌎))
  e′² = squared_eccentricity_earth / (1 - squared_eccentricity_earth)
  ψ = atand(majoraxis_earth * z, minoraxis_earth * p)
  ϕ = mod1(atand(z + minoraxis_earth * e′² * sind(ψ)^3, p - majoraxis_earth * squared_eccentricity_earth * cosd(ψ)^3),360)

  N = majoraxis_earth / sqrt(1 - squared_eccentricity_earth * sind(ϕ)^2)
  ## Fix for the condition cosθ=0, in that case subtract the
  cosθ=cosd(ϕ)
  @debug "cosθ=$cosθ"
  if abs(cosθ)>atol(T)
    h = p/cosθ - N
  else
    h =abs(z)-minoraxis_earth
  end

  return (ϕ,h)
end


function Base.convert(::Type{ECEF2D{Datum}},coords::LLA2D{Datum,T}) where {T,Datum}
  h=ustrip(coords.h)
  θ=ustrip(coords.θ)
  return ECEF2D{Datum}(_from_lla2d_to_ecef2d(Datum,θ,h)...)
end

function Base.convert(::Type{LLA2D{Datum}},coords::ECEF2D{Datum,T}) where {T,Datum}
  z = ustrip(coords.z)
  p = ustrip(coords.w)
  return LLA2D{Datum}(_from_ecef2d_to_lla2d(Datum,z,p)...)
end
