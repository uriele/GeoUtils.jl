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
const TMINC= uconvert(u"°C",0K) |> ustrip

const ALTTROPOSPHERE=11
const ALTSTRATOSPHERE1=20
const ALTSTRATOSPHERE2=32
const ALTSTRATOPAUSE=47
const ALTMESOSPHERE1=51
const ALTMESOSPHERE2=71
const ALTMESOPAUSE=85
#
# - Implement a function that computes the temperature at a given height using the
# ENVIROMENTAL LAPSE RATE (ELR)
# - I do not use stratopause and mesopause because their lapse rates is zero.
function _lapse_rate(h::T)::T where T
        (h<ALTTROPOSPHERE)*TALT_TROPOSPHERE+
        (ALTSTRATOSPHERE1<h<ALTSTRATOSPHERE2)*TALT_STRATOSPHERE1+
        (ALTSTRATOSPHERE2<=h<ALTSTRATOPAUSE)*TALT_STRATOSPHERE2+
        (ALTMESOSPHERE1<h<ALTMESOSPHERE2)*TALT_MESOSPHERE1+
        (ALTMESOSPHERE2<=h<=ALTMESOPAUSE)*TALT_MESOSPHERE2 |> x-> T(x)
end

h=LinRange(0,120,1000)
dh=diff(h)[1]
T0=40.
T=[T0]
for _h in h
  push!(T,T[end]-_lapse_rate(_h)*dh)
end

# - Assume the temperature varies linearly with the height with a lapse rate of -6.5K/km
# in the troposphere (up to 11km) and is constant in the stratosphere
# - Assume it increses in the stratosphere with a rate of 1K/km
# - Assume the temperature changes with the latitude as Tlat=Tmin+(Tmax-Tmin)*cos(θ) with
# minimum temperature at the poles and maximum at the equator

function _t(h::T,θ::T;Tmin=-10,Tmax=30,hmin=0,hmax=120,dh=0.01)::T where T<:IEEEFloat
  Tlat=Tmin+(Tmax-Tmin)*cosd(θ)
  h= min(max(hmin,h),hmax)
 # @info h
  hh=divrem(h,dh) |> x-> (x[1]*dh,x[2])
#  @info hh
  Temp=reduce(+,[-_lapse_rate(_h)*dh for _h in 0:dh:hh[1]])+_lapse_rate(h)*hh[2]
#  @info Temp
  return Temp.+Tlat
end
_t(h::N1,θ::N2;kwargs...) where {N1<:Number,N2<:Number} =_t(float.((h,θ))...;kwargs...)

# - Uses the barometric formula to compute the pressure as a function of the height and temperature
function _p(h::T,Temperature::T;Pb=101330,dh=0.01)::T where T<:IEEEFloat
  GRAVITY= 9.80665  #m/s^2
  R= 8.3144598 #J/(mol*K)
  MAIR=0.0289644 #kg/mol
  C=GRAVITY*MAIR/R

  Tb=uconvert(u"K",Temperature*u"°C") |> ustrip
  hb=uconvert(u"m",h*u"km") |> ustrip
#  @info hb
#  @info Tb
  P0=101325 #Pa
  return P0*exp(-C*hb/Tb)
end

_p(h::N1,Temperature::N2;kwargs...) where {N1<:Number,N2<:Number} =_p(float.((h,Temperature))...;kwargs...)

function EarthModel!(ax,h,θ;a=1,eccentricity=0.0818191908426215)
  b=a*sqrt(1-eccentricity^2)
  for _h in h
    lines!(ax,[((a+_h)*cosd(_θ),(b+_h)*sind(_θ)) for _θ in LinRange(0,361,720)],color=:black)
  end
  for _θ in θ
    lines!(ax,[((a+h[1])*cosd(_θ),(b+h[1])*sind(_θ));
    ((a+h[end])*cosd(_θ),(b+h[end])*sind(_θ))
     ],color=:black)
  end
end

@kwdef struct WedgeAtm{T}
  center::Point2f=Point2f(0.,0.)
  h₁::T=0.0
  h₂::T=1.0
  θ₁::T=0.0
  θ₂::T=90.0
  a ::T=1.0
  b ::T=0.9
end

@inline _ellipse_point(a,b,h,θ)=Point2f((a+h)*cosd(θ),(b+h)*sind(θ))
function mycoordinates(w::WedgeAtm{T},nvertex=60) where T
  p=Point2f[]

  for θ in LinRange(w.θ₁,w.θ₂,nvertex)
    push!(p,_ellipse_point(w.a,w.b,w.h₁,θ))
  end
  for θ in LinRange(w.θ₂,w.θ₁,nvertex)
    push!(p,_ellipse_point(w.a,w.b,w.h₂,θ))
  end
  return p
end

a=majoraxis(ellipsoid(WGS84)) |> x-> uconvert(km,x) |> ustrip
b=minoraxis(ellipsoid(WGS84)) |> x-> uconvert(km,x) |> ustrip |> x-> x/a

θ=convert_to_array("data_atm/INP_FILES/in_lat.dat") |> x-> [x;x[1]+360].+90
h1=convert_to_array("data_atm/INP_FILES/in_alt.dat",16)
temp=convert_to_array("data_atm/INP_FILES/in_temp.dat",13) |> x-> reshape(x,(length(h1),length(θ[1:end-1]))) |> x-> x[1:end-1,1:end] |>
x-> uconvert.(u"°C",x.*K) |> x-> ustrip.(x)
press=convert_to_array("data_atm/INP_FILES/in_pres.dat",13) |> x-> reshape(x,(length(h1),length(θ[1:end-1]))) |> x-> x[1:end-1,1:end] |>
x-> uconvert.(u"Pa",x.*u"mbar") |> x-> ustrip.(x)
h=h1./a
pp=Vector{Point2f}[]
[push!(pp,mycoordinates(WedgeAtm(b=b,θ₁=θ₁,θ₂=θ₂,h₁=h₁,h₂=h₂))) for (h₁,h₂) in zip(h[1:end-1],h[2:end]),(θ₁,θ₂) in zip(θ[1:end-1],θ[2:end])];

hh=h1[1:end-1]+diff(h1)/2
tt=θ[1:end-1]+diff(θ)/2
ref=[refractive_index(Mathar4();wavelength=10.0,temperature=temp[i,j],pressure=press[i,j]) for i in 1:size(temp,1), j in 1:size(temp,2)]
ref1=[refractive_index(;wavelength=10.0,temperature=temp[i,j],pressure=press[i,j],CO2ppm=0.0) for i in 1:size(temp,1), j in 1:size(temp,2)]

ref2=copy(ref1)
ref2[hh.>=100,:].=repeat(ref2[findfirst(hh.<=100),:]',20,1)




size(ref[:])
fig=Figure(size=(600,800))
ax,p=poly(fig[1,1][1,1],pp[:],color=ref[:])
Colorbar(fig[1,1][1,2],p)
lines(fig[2,1][1,1], hh, temp[:,div(end,2)],color=:black)
lines(fig[2,1][2,1], hh, press[:,div(end,2)],color=:black)
ax,lin=lines(fig[2,1][3,1], hh, ref[:,div(end,2)],color=:black,label="Mathar")
lines!(ax, hh, ref1[:,div(end,2)],color=:red,label="Ciddor")
lines!(ax, hh, ref2[:,div(end,2)],color=:red,label="Ciddor Tover")
axislegend(ax,title="Refractive index",fontsize=10)

n=SemiCircularMatrix(reshape(ref,size(temp)))


heatmap(hh,tt,n)
heatmap(hh,tt,temp)
heatmap(hh,tt,press)

fig=Figure()
ax,lin=lines(fig[1,1], hh, ref[:,div(end,2)],color=:black,label="Mathar")
lines!(ax, hh, ref1[:,div(end,2)],color=:red,label="Ciddor")
lines!(ax, hh, ref2[:,div(end,2)],color=:blue,label="Ciddor Tover")
axislegend(ax,title="Refractive index",fontsize=10)
xlims!(ax,100,120)


fig=Figure()
ax,pl=lines(fig[1,1], hh, press[:,div(end,2)],color=:red)
xlims!(ax,80,120)
ylims!(ax,-1000,150000)

refractive_index(Mathar4();wavelength=10.0,temperature=1020.0,pressure=10.0,CO2ppm=0.0)
