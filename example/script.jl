using GeoUtils,WGLMakie,Makie
using Unitful
using Unitful:°
using CoordRefSystems
using CoordinateTransformations
using Base: IEEEFloat
const TALT_TROPOSPHERE=6.5     # Lapse rate in K/km
const TALT_STRATOSPHERE1=-1     # Lapse rate in K/km
const TALT_STRATOSPHERE2=-2.8     # Lapse rate in K/km
const TALT_MESOSPHERE1=2.8     # Lapse rate in K/km
const TALT_MESOSPHERE2=2.0     # Lapse rate in K/km
const TDELTAALT= 1  # 1km
const TTHETA= x-> cosd(x)
const TMINC= uconvert(u"°C",0K) |> ustrip
#
# - Implement a function that computes the temperature at a given height using the
# ENVIROMENTAL LAPSE RATE (ELR)
# - I do not use stratopause and mesopause because their lapse rates is zero.
function _lapse_rate(h::T)::T where T
        (h<11)*TALT_TROPOSPHERE+
        (20<h<32)*TALT_STRATOSPHERE1+
        (32<=h<47)*TALT_STRATOSPHERE2+
        (51<h<71)*TALT_MESOSPHERE1+
        (71<=h<=85)*TALT_MESOSPHERE2 |> x-> T(x)
end

h=LinRange(0,120,1000)
dh=diff(h)[1]
T0=40.
T=[T0]
for _h in h
  push!(T,T[end]-_lapse_rate(_h)*dh)
end
fig=Figure(aspect_ratio=1)
ax=Axis(fig[1,1],xlabel="Temperature (°C)",ylabel="Height (km)")
lines!(ax,T,vcat(0,h),color=:blue)

# - Assume the temperature varies linearly with the height with a lapse rate of -6.5K/km
# in the troposphere (up to 11km) and is constant in the stratosphere
# - Assume it increses in the stratosphere with a rate of 1K/km
# - Assume the temperature changes with the latitude as Tlat=Tmin+(Tmax-Tmin)*cos(θ) with
# minimum temperature at the poles and maximum at the equator

_position(1,10)
function _t(h::T,θ::T;Tmin=-10,Tmax=30,hmin=0,hmax=120)::T where T<:IEEEFloat
  Tlat=Tmin+(Tmax-Tmin)*cosd(θ)
  h= min(max(hmin,h),hmax)
  @info h
  Temp=Tlat+h*TALT/TDELTAALT # |> t-> max(TMINC,t)
  return Temp
end
_t(h::N1,θ::N2;kwargs...) where {N1<:Number,N2<:Number} =_t(float.((h,θ))...;kwargs...)

_t(80,0;Tmax=50)
# - Uses the barometric formula to compute the pressure as a function of the height and temperature
function _p(h::T,Temperature::T;Pb=101330,hb=0,hmax=120)::T where T<:IEEEFloat
  Pb=
  Pl=Pmin+(Pmax-Pmin)*cosd(θ)
  return Pl*exp(-h*TALT/TDELTAALT)
end
