using CoordRefSystems: Deg,Rad
using CoordRefSystems: Geographic
using CoordRefSystems: Datum
using Unitful
using Unitful: Length
import Unitful.° as deg
import Unitful.km as km
using CoordRefSystems
const KM{T}=Quantity{T,u"𝐋",typeof(km)}
struct EarthCenteredInertial{Datum,D<:Deg,L<:KM} <: Geographic{Datum}
  altitude::D
  azimuth ::D
  radius  ::L
end

isa(1deg,Deg)
EarthCenteredInertial{Datum}(altitude::Deg,azimuth::Deg,radius::KM) where {Datum}=
EarthCenteredInertial{Datum,float(Deg),float(KM)}(altitude,azimuth,radius)


EarthCenteredInertial{Datum}(altitude::Deg,azimuth::Deg,radius::Length) where {Datum}=
  EarthCenteredInertial{Datum}(promote(altitude,azimuth)...,radius)
EarthCenteredInertial{Datum}(altitude::Deg,azimuth::Deg,radius::Number) where {Datum}=
  EarthCenteredInertial{Datum}(promote(altitude,azimuth)...,KM(radius))
EarthCenteredInertial{Datum}(altitude::Rad,azimuth::Rad,radius::Length) where {Datum}=
  EarthCenteredInertial{Datum}(rad2deg(altitude),rad2deg(azimuth),radius)
EarthCenteredInertial{Datum}(altitude::Rad,azimuth::Rad,radius::Number) where {Datum}=
  EarthCenteredInertial{Datum}(rad2deg(altitude),rad2deg(azimuth),KM(radius))
EarthCenteredInertial{Datum}(altitude::Deg,azimuth::Rad,radius::Length) where {Datum}=
  EarthCenteredInertial{Datum}(altitude,rad2deg(azimuth),radius)
EarthCenteredInertial{Datum}(altitude::Deg,azimuth::Rad,radius::Number) where {Datum}=
  EarthCenteredInertial{Datum}(altitude,rad2deg(azimuth),KM(radius))
EarthCenteredInertial{Datum}(altitude::Rad,azimuth::Deg,radius::Length) where {Datum}=
  EarthCenteredInertial{Datum}(rad2deg(altitude),azimuth,radius)
EarthCenteredInertial{Datum}(altitude::Rad,azimuth::Deg,radius::Number) where {Datum}=
  EarthCenteredInertial{Datum}(rad2deg(altitude),azimuth,KM(radius))

EarthCenteredInertial{Datum}(altitude::Number,azimuth::Number,radius::Length) where {Datum}=
  EarthCenteredInertial{Datum}(Deg(altitude),Deg(azimuth),radius)
EarthCenteredInertial{Datum}(altitude::Number,azimuth::Number,radius::Number) where {Datum}=
  EarthCenteredInertial{Datum}(Deg(altitude),Deg(azimuth),KM(radius))

EarthCenteredInertial(altitude,azimuth,radius) where {Datum}=
  EarthCenteredInertial{WGS84Latest}(altitude,azimuth,radius)

const ECI{Datum}=EarthCenteredInertial{Datum}

struct EarthCenteredInertial2D{Datum,D<:Deg,L<:KM} <: Geographic{Datum}
  altitude::D
  h  ::L
end

EarthCenteredInertial2D{Datum}(altitude::Deg,radius::KM) where {Datum}=
  EarthCenteredInertial2D{Datum,float(Deg),float(KM)}(altitude,radius)
EarthCenteredInertial2D{Datum}(altitude::Deg,radius::Number) where {Datum}=
  EarthCenteredInertial2D{Datum}(altitude,KM(radius))
EarthCenteredInertial2D{Datum}(altitude::Rad,radius) where {Datum}=
  EarthCenteredInertial2D{Datum}(deg2rad(altitude),radius)
EarthCenteredInertial2D{Datum}(altitude::Number,radius) where {Datum}=
  EarthCenteredInertial2D{Datum}(Deg(altitude),radius)

EarthCenteredInertial2D(altitude,radius) where {Datum}=
  EarthCenteredInertial2D{WGS84Latest}(altitude,radius)

const ECI2D{Datum}=EarthCenteredInertial2D{Datum}

ECI2D(1,1)
ECI(1,1,1)
Base.convert(::Type{EarthCenteredInertial{Datum,D,L}},coord::EarthCenteredInertial{Datum}) where {Datum,D,L}=
  EarthCenteredInertial{Datum,D,L}(coord.altitude,coord.azimuth,coord.radius)
ECI{WGS84}(ECI(1,1,1).altitude,
ECI(1,1,1).azimuth,
ECI(1,1,1).radius)
convert(EarthCenteredInertial,ECI(1.,1.,1.))

using CoordRefSystems: raw,constructor,reconstruct
raw(coords::EarthCenteredInertial) = ustrip(coords.altitude), ustrip(coords.azimuth), ustrip(coords.radius)

CoordRefSystems.constructor(::Type{<:EarthCenteredInertial{Datum}}) where {Datum} = EarthCenteredInertial{Datum}
CoordRefSystems.reconstruct(C::Type{EarthCenteredInertial{Datum}}, raw) where {Datum} = begin
  (altitude,azimuth,radius) = raw.*units(C)
  constructor(C)(altitude,azimuth,radius)
end
