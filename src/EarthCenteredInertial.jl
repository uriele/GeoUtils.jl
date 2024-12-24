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
  a=majoraxis(ellipsoid(Datum)) |> x-> uconvert(km,x) |> ustrip
  e²=eccentricity²(ellipsoid(Datum))
  h=rand(rng)*100;
  ϕ=-90+180*rand(rng)
  γ=-180+360*rand(rng)
  N=a/sqrt(1-e²*sin(ϕ)^2)
  EarthCenteredEarthFixed{Datum}(
    (N+h)*cosd(ϕ)*cosd(γ),
    (N+h)*cosd(ϕ)*sind(γ),
    (N*(1-e²)+h)*sind(ϕ)
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
  a=majoraxis(ellipsoid(Datum)) |> x-> uconvert(km,x) |> ustrip
  EarthCenteredInertial{Datum}(
    -90+180*rand(rng),
    -180+360*rand(rng),
    a+800*abs(rand(rng))
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
  a=majoraxis(ellipsoid(Datum)) |> x-> uconvert(km,x) |> ustrip
  EarthCenteredInertial2D{Datum}(
    360*rand(rng),
    a+800*abs(rand(rng))
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
#  (1-e²) * tan(ϕ) = (1-e²) * tan(ϕ′)
#
#


function Base.convert(::Type{EarthCenteredEarthFixed{Datum}},coords::LatLonAlt{Datum}) where {Datum}
  🌎  = ellipsoid(Datum)
  ϕ′  = ustrip(coords.lat)
  e²  = oftype(ϕ′, eccentricity²(🌎))
  a   = majoraxis(🌎) |> x-> uconvert(km,x) |> ustrip
  N   = a / sqrt(1 - e² * sind(ϕ′)^2)

  h   = coords.alt |> x-> uconvert(km,x) |> ustrip

  EarthCenteredEarthFixed{Datum}(
    (N+h)*cosd(ϕ′)*cosd(coords.lon),
    (N+h)*cosd(ϕ′)*sind(coords.lon),
    (N*(1-e²)+h)*sind(ϕ′))
end

# Use Zhu Algorithm Closed Form
#
function Base.convert(::Type{LatLonAlt{Datum}},coords::EarthCenteredEarthFixed{Datum}) where {Datum}
  🌎  = ellipsoid(Datum)

  X = coords.E |> x-> uconvert(km,x) |> ustrip
  Y = coords.F |> x-> uconvert(km,x) |> ustrip
  Z = coords.G |> x-> uconvert(km,x) |> ustrip
  λ = atand(Y/X)

  e²  = oftype(λ, eccentricity²(🌎))
  a   = majoraxis(🌎) |> x-> uconvert(km,x) |> ustrip
  P   = hypot(X,Y)

  # initial value
  ϕ′ = atand(P/Z)

  @debug "ϕ′: $(ϕ′)"

  @inline function N(ϕ′)
    a / sqrt(1-e²*sind(ϕ′)^2)
  end

  N′  = N(ϕ′)

  while true
    ϕ = atand(Z/P/(1-e²*N′*cosd(ϕ′)/P))
    err=abs(ϕ-ϕ′)
    #@debug "ϕ: $(rad2deg(ϕ)), N: $N′, err: $(rad2deg(err))"
    ϕ′ = ϕ
    N′= N(ϕ′)
    if err < eps()
      break
    end
  end
  h′ = P/cosd(ϕ′)-N′

  LatLonAlt{Datum}(Deg(ϕ′),Deg(λ),KM(h′))

end
