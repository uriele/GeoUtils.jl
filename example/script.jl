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
fig=Figure(aspect_ratio=1)
ax=Axis(fig[1,1],xlabel="Temperature (°C)",ylabel="Height (km)")
lines!(ax,T,vcat(0,h),color=:blue)

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


N= 20 # Number of elevations
M= 120 # Number of radii

h=LinRange(4,120,N)
θ=LinRange(0,360,M+1)[1:end-1]

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

fig=Figure(aspect_ratio=1)
ax=Axis(fig[1,1],xlabel="x (km)",ylabel="y (km)")
a=majoraxis(ellipsoid(WGS84)) |> x-> uconvert(km,x) |> ustrip
EarthModel!(ax,h,θ;a=50)

a=50
eccentricity= 0.0818191908426215
b=a*sqrt(1-eccentricity^2)
satellite= cosd(45)*(a+140),sind(45)*(b+140)
tang_v(x,y,h;a=a,b=b)=2.0.*(x/(a+h)^2,y/(b+h)^2)

tang_v(satellite...,140)

atand(tang_v(satellite...,140)...)


x0=a*cosd(45)
x1=(a+140)*cosd(45)
y0=b*sind(45)
y1=(b+140)*sind(45)

(y1-y0)/(x1-x0)

limb_angle= 20


scatter!(ax,[satellite],color=:red)

a=50
b=48
N(t)=a^2/sqrt(a^2*cosd(t)^2+b^2*sind(t)^2)
N1(t,h)=(a+h)^2/sqrt((a+h)^2*cosd(t)^2+(b+h)^2*sind(t)^2)
X(t,h=0)=(N(t)+h)*cosd(t)
Y(t,h=0)=(b^2/a^2*N(t)+h)*sind(t)

X1(t,h=0)=N1(t,h)*cosd(t)
Y1(t,h=0)=(b+h)^2/(a+h)^2*N1(t,h)*sind(t)

X2(t,h=0)=(a+h)*cosd(t)
Y2(t,h=0)=(b+h)*sind(t)

[hypot(X(45,h)-X1(45,h),Y(45,h)-Y1(45,h)) for h in LinRange(0,100,1000)] |> x-> extrema(x)
[hypot(X(25)-X(25,h),Y(25)-Y(25,h)) for h in LinRange(0,100,1000)] |> x-> extrema(x)
[hypot(X1(25,0)-X1(25,t),Y1(25,0)-Y1(25,t)) for t in LinRange(0,100,1000)] |> x-> extrema(x)
[hypot(X2(25,0)-X2(25,t),Y2(25,0)-Y2(25,t)) for t in LinRange(0,100,1000)] |> x-> extrema(x)

fig=Figure(aspect_ratio=1)
ax=Axis(fig[1,1],xlabel="x (km)",ylabel="y (km)")
lines!(ax,[(X(t,10),Y(t,10)) for t in LinRange(0,360,720)],color=:black)
lines!(ax,[(X1(t,10),Y1(t,10)) for t in LinRange(0,360,720)],color=:blue,linestyle=:dash)
lines!(ax,[(X2(t,10),Y2(t,10)) for t in LinRange(0,360,720)],color=:red,linestyle=:dot)

xlims!(ax,-60.5,-59.5)
ylims!(ax,-0.5,0.5)

using Symbolics
@variables x y h t a b e²

N(t)=a^2/(a^2*cosd(t)^2-b^2*sind(t)^2)^(1/2)
X(t)=(N(t)+h)*cos(t)
Y(t)=((1-e²)*N(t)+h)*sin(t)

Dt=Differential(t)
Dx=Differential(x)
Dy=Differential(y)

expand_derivatives(Dx(x^2/a^2+y^2/b^2))

simplify(expand_derivatives(Dt(N(t)))/N(t))

simplify(sin(t)*cos(t)/2)

simplify(X(t)^2/a^2+Y(t)^2/b^2)

simplify(N(t)^2)


expand_derivatives(Dt(1/sqrt(cos(t)^2+(1-e²)*sin(t)^2)))


lines([begin
  (M1*[(a+h)*cosd(tt);(b+h)*sind(tt)]) |>
  x-> (x[1]/(1+h/a),x[2]/(1+h/b)) |>
  x-> sum(x-> x^2,x) |> sqrt |> x-> (x[1]*cosd(tt),x[2]*sind(tt))
end
  for tt in 0:360])

a=50
b=48

M=[a 0;0 b]
M1=inv(M)
h=10
lines!([Tuple(M*[(1+h/a)*cosd(tt);(1+h/b)*sind(tt)]) for tt in 0:360],color=:black)
a= 100
b=50
phi(t,h)=atand((N(t)*b^2/a^2+h)/(N(t)+h)*tand(t))
X(t,h=0)=(a+h)*cosd(t)
Y(t,h=0)=(b+h)*sind(t)
N(t)=a^2/sqrt(a^2*cosd(t)^2+b^2*sind(t)^2)
X1(t,h=0)=(N(t)+h)*cosd(t)
Y1(t,h=0)=(b^2/a^2*N(t)+h)*sind(t)


a
b
h
hypot(a*cosd(45)-(a+h)*cosd(45),b*sind(45)-(b+h)*sind(45))


fig=Figure(aspect_ratio=1)
ax=Axis(fig[1,1],xlabel="x (km)",ylabel="y (km)")
lines!(ax,[(X(t),Y(t)) for t in LinRange(0,360,720)],color=:black)
lines!(ax,[(X(t,10),Y(t,10)) for t in LinRange(0,360,720)],color=:black)
lines!(ax,[(X1(t),Y1(t)) for t in LinRange(0,360,720)],color=:red,linestyle=:dash)
lines!(ax,[(X1(t,10),Y1(t,10)) for t in LinRange(0,360,720)],color=:blue,linestyle=:dash)


hypot(X(45,0)-X(45,10),Y(45,0)-Y(45,10))


norm_vec(t,h)= [X(t,h)/(a+h)^2,Y(t,h)/(b+h)^2] |> x-> (x[1],x[2])./hypot(x[2],x[1])

norm_vec1(t,h)= [X1(t,h)/(a+h)^2,Y1(t,h)/(b+h)^2] |> x-> atand(x[2],x[1])

atand(norm_vec(45,0)[2:-1:1]...)
atand(norm_vec(45,100)[2:-1:1]...)
norm_vec1(45,0)
norm_vec1(45,100)

r(0,45).-(X(45,0),Y(45,0))
r(s,t,h=0)=(X(t,h),Y(t,h)) |> x -> x.+s.*norm_vec(t,h)
TT(s,t,h=0)= (X(t,h),Y(t,h)) |> x -> x.+s.*(-norm_vec(t,h)[2],norm_vec(t,h)[1])

fig=Figure(aspect_ratio=1)
ax=Axis(fig[1,1],xlabel="x (km)",ylabel="y (km)")
lines!(ax,[(X(t),Y(t)) for t in LinRange(0,361,720)],color=:black)
lines!(ax,[r(s,45,0) for s in LinRange(-10,10,100)],color=:black)
lines!(ax,[TT(s,45,0) for s in LinRange(0,100,1000)],color=:black)
lines!(ax,[TT(s,180-45,0) for s in LinRange(-100,0,1000)],color=:black)
lines!(ax,[TT(s,5,0) for s in LinRange(-20,100,1000)],color=:black)
lines!(ax,[TT(s,180-5,0) for s in LinRange(-100,20,1000)],color=:black)
lines!(ax,[r(s,5,0) for s in LinRange(-20,100,1000)],color=:black)
lines!(ax,[r(s,180-5,0) for s in LinRange(-20,100,1000)],color=:black)


using LinearAlgebra

dot([TT(0,45)...],[r(0,45)...])

TT(1,45)
r(1,45)

lines([dot([r(1,t).-r(0,t)],[TT(1,t).-TT(0,t)]) for t in LinRange(0,360,720)])


h=0.0

expand(TT(t,45,h)[1]^2*(b+h)^2+TT(t,45,h)[2]^2*(a+h)^2-(a+h)^2*(b+h)^2)

(A,B,C)=expand(TT(t,45,h)[1]^2*(b+h)^2+TT(t,45,h)[2]^2*(a+h)^2-(a+h)^2*(b+h)^2) |> x-> Symbolics.coeff.(x,[t^2,t,1])

B
C
A

((B/A)^2-4*C/A)
