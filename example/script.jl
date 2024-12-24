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
 # @ h
  hh=divrem(h,dh) |> x-> (x[1]*dh,x[2])
#  @ hh
  Temp=reduce(+,[-_lapse_rate(_h)*dh for _h in 0:dh:hh[1]])+_lapse_rate(h)*hh[2]
#  @ Temp
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
#  @ hb
#  @ Tb
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


using Makie,WGLMakie

fig=Figure(aspect_ratio=1)
θ=LinRange(0,10,10)

h=LinRange(0,50,50)


dh=diff(h)[1]./68765
dθ=diff(θ)[1]
ax=Axis(fig[1,1],xlabel="x (km)",ylabel="y (km)")
pp=Vector{Point2f}[]
[push!(pp,mycoordinates(WedgeAtm(b=0.5,θ₁=θ,θ₂=θ+dθ,h₁=h,h₂=h+dh))) for h in h./68765,θ in θ];


temp=[_t(h,θ) for h in h,θ in θ ]
pres=[_p(h,_t(h,θ)) for h in h,θ in θ ]
ref=[refractive_index(temperature=_t(h,θ),pressure=_p(h,_t(h,θ))) for h in h,θ in θ ]

poly(pp,color=pres[:])
poly(pp,color=temp[:])
fig=Figure(aspect_ratio=1)
ax,p=poly(fig[1,1],pp,color=ref[:])
Colorbar(fig[1,2],p)

fig=Figure(aspect_ratio=1)
ax,p=poly(fig[1,1],pp,color=temp[:])
Colorbar(fig[1,2],p)
lines(fig[2,1], h, temp[:,div(end,2)],color=:black)
lines(fig[3,1], h, pres[:,div(end,2)],color=:black)
lines(fig[4,1], h, ref[:,div(end,2)],color=:black)


PolarAxis(fig[1,1],h, temp[:,div(end,2)],color=:black)


datum= WGS84Latest
🌎 = ellipsoid(datum)
ϕ = ustrip.(deg2rad(float(10deg)))
e = oftype(ϕ, eccentricity(🌎))
e² = oftype(ϕ, eccentricity²(🌎))

ϕ′ = oftype(ϕ, deg2rad(10.0))

ϕa=atan(1/(1-e²)*tan(ϕ′))

ϕb=ϕ′

a= majoraxis(🌎) |> x-> uconvert(km,x) |> ustrip

_C(x)=1/sqrt(1-e²*sin(x)^2)


R=hypot(a*cos(ϕ′),b*sin(ϕ′))

h=1000
ϕb=ϕ′


fromGtoE(ϕb,h)=(1-e²*_NN(ϕb)/(_NN(ϕb)+h))
_NN=ϕ->a/sqrt(1-e²*sin(ϕ)^2)

fromEtoG(ϕb,h)=fromGtoE(ϕb,h)^-1

GtoE(ϕb,h)=atan(fromGtoE(ϕb,h)*tan(ϕb))

EtoG(ϕ,h)= begin
  ϕ′=ϕ
  diff=1
  count=0
  while diff>eps()
    ϕ″=atan(fromEtoG(ϕ′,h)*tan(ϕ))
    diff=abs(ϕ″-ϕ′)
    ϕ′=ϕ″
    count+=1
  end
   count
  return ϕ′
end

lines!(LinRange(0,2*pi,360),[rad2deg(GtoE(EtoG((x),1000),1000)) for x in LinRange(0,2*pi,360)])

using Makie,WGLMakie

lines(LinRange(0,2*pi,360),[atand(tan(x)) for x in LinRange(0,2*pi,360)])


#=
let
  figure=Figure()
  @ h,b

  ax=Axis(figure[1,1])
  lines!(ax,tt,[(X(0,t,b))^2+(Y(0,t,b))^2/b^2 for t in tt])
  lines!(ax,tt,[(X(h,t,b))^2/(1+h)^2+(Y(h,t,b))^2/(b+h)^2-(Y(h,t,b))*cosd(t)/(b+h)-(X(h,t,b))*sind(t)/(1+h)
  for t in tt])
  figure
end


=#

#=
x0,y0=(-1,-1)
dx0,dy0=(-y0/b^2,x0) |> x-> x./hypot(x...)
R(t)=[cosd(t) -sind(t);sind(t) cosd(t)]
dx1,dy1=R(45)*[dx0,dy0]
r(t1)=[x0+dx1*t1,y0+dy1*t1]
r1(t1)=[x0+x0*t1,y0+y0/b^2*t1]
r2(t1)=[x0+dx0*t1,y0+dy0/b^2*t1]
h=120000/a+0.3

AA=dx1^2/(1+h)^2+dy1^2/(b+h)^2
BB=2*((x0*dx1)/(1+h)^2+(y0*dy1)/(b+h)^2)
CC=(x0)^2/(1+h)^2+(y0)^2/(b+h)^2-1

hd=hypot(dx1/(1+h),dy1/(b+h))
ddx1=dx1/(1+h)
ddy1=dy1/(b+h)
(ddx1,ddy1)=(ddx1,ddy1)./hd
xx0=x0/(1+h)
yy0=y0/(b+h)
hypot(ddx1,ddy1)
AA1=ddx1^2+ddy1^2
BB1=2*((xx0*ddx1)+(yy0*ddy1))
CC1=(xx0)^2+(yy0)^2-1


dd=(BB)^2-4*AA*CC
t1=(-BB-sqrt(dd))/(2*AA)



dd1=(BB1/2)^2-CC1
tt1=(-BB1/2-sqrt(dd1))
tt1=tt1/hd

function better_interception(x0,y0,)

figure=Figure()
ax=Axis(figure[1,1])
tt=LinRange(-1,2,100)
lines!(ax,[Point2(r(t)) for t in tt])
lines!(ax,[Point2(r1(t)) for t in tt])
lines!(ax,[Point2(r2(t)) for t in tt])

lines!(ax,[(X(h,t,b), Y(h,t,b)) for t in LinRange(0,361,721)])
lines!(ax,[((1+h)*cosd(t), (b+h)*sind(t)) for t in LinRange(0,361,721)])
scatter!(ax,(r(t1)[1],r(t1)[2]))
scatter!(ax,(r(tt1)[1],r(tt1)[2]))

scatter!(ax,[(r(t1)[1],-sqrt(1-r(t1)[1]^2/(1+h)^2)*(b+h))])
r(t1)[2]+sqrt(1-r(t1)[1]^2/(1+h)^2)*(b+h)
r(tt1)[2]+sqrt(1-r(tt1)[1]^2/(1+h)^2)*(b+h)

e22=1-b^2
X(h,t,b)^2/(1)^2+Y(h,t,b)^2/(b)^2

2*h/sqrt(b^2*X(h,t,b)^2+Y(h,t,b)^2)
h



let
bb=a*0.5
ee2=1-b^2
t=0

function conic_coeff(points)
  A=zeros(points,5)
  [A[i,:]=[x^2 y^2 x y 1] for (i,(x,y)) in enumerate(points)]

  _,_,V=svd(A)
  return V[:,end] |> x -> x./x[end] |> x -> x[abs.(x).<1e-10]=0
end
h=a
Na(t)=a^2/sqrt(a^2*cosd(t)^2+bb^2*sind(t)^2)
Nb(t)=bb^2/sqrt(a^2*cosd(t)^2+bb^2*sind(t)^2)

@ conic_coeff([(Na(t),Nb(t)) for t in LinRange(0,90,10)])

for t in LinRange(0,90,10)

  xx,yy=((Na(t)+h)*cosd(t),(Nb(t)+h)*sind(t))
  @ t,xx^2/(a+h)^2+yy^2/(bb+h)^2,
  xx^2/(a+h)^2+yy^2/(bb+h)^2-
  2*h/a^2*yy-2*h/bb^2*xx
end

for t in LinRange(0,90,10)
  @ t,(a^2)*cosd(t)^2/a^2+(bb^2)*sind(t)^2/bb^2
end
end



function conic_coeff(points)
  A=zeros(length(points),6)
  [A[i,:]=[x^2 y^2 x y x*y 1] for (i,(x,y)) in enumerate(points)]
  @ "A $A"
  _,_,V=svd(A)
  @ "V $(V[:,5])"
  return V[:,5] |> x -> x./x[end]
end
let
  h=800000
  bb=b*a

  Na(t)=a^2/sqrt(a^2*cosd(t)^2+bb^2*sind(t)^2)
  Nb(t)=bb^2/sqrt(a^2*cosd(t)^2+bb^2*sind(t)^2)
  Nah(t,h)=(a+h)^2/sqrt((a+h)^2*cosd(t)^2+(bb+h)^2*sind(t)^2)
  Nbh(t,h)=(bb+h)^2/sqrt((a+h)^2*cosd(t)^2+(bb+h)^2*sind(t)^2)

  cc=conic_coeff([((Na(t)+h)*cosd(t),(Nb(t)+h)*sind(t)) for t in LinRange(0.1,89.9,6)])
  (x,y)=((Na(45)+h)*cosd(45),(Nb(45)+h)*sind(45))
  (x1,y1)=((Na(45))*cosd(45),(Nb(45))*sind(45))
  (x2,y2)=((Nah(45,h))*cosd(45),(Nbh(45,h))*sind(45))
  @ cc
  @ atand(y/(bb+h)^2,x/(a+h)^2)
  @ atand(y2/(bb+h)^2,x2/(a+h)^2)
   @ extrema([let
   ((Na(t)+h)*cosd(t),(Nb(t)+h)*sind(t),(Nah(t,h))*cosd(t),(Nbh(t,h))*sind(t)) |>
   x-> hypot(x[1]-x[3],x[2]-x[4])
  end
  for t in LinRange(0,90,900)])


  @ atand(y1/(bb)^2,x1/(a)^2)
  @ atand(y2/(bb+h)^2,x2/(a+h)^2)
  @ extrema([let
  ((Na(t))*cosd(t),(Nb(t))*sind(t),(Nah(t,h))*cosd(t),(Nbh(t,h))*sind(t)) |>
  x-> hypot(x[1]-x[3],x[2]-x[4])
 end
 for t in LinRange(0,90,900)])


end

=#

#read_local_atmosphere("data_atm/INP_FILES";normalize=true)

#orbit=read_orbit("data_atm/INP_FILES/orbit.dat")
#orbit.w
#_,_,n,da=discretize_atmosphere1(atmosphere,100,100)
#da
#n

figure=Figure()
Ax=Axis(figure[1,1])
for j in axes(da,2)
  lines!(Ax,[(w,z) for (w,z) in da[:,j]])
end
for i in axes(da,1)
  lines!(Ax,[(w,z) for (w,z) in da[i,1:end+1]])
end

figure
scatter([(w,z) for (w,z) in zip(orbit.w[:],orbit.z[:])])

let

  h=800000/a
  tt=LinRange(0,360,100)
  figure=Figure(size=(1000,1000))
  ax=Axis(figure[1:2,1:2][1,1],title="b/a: $(round(b,digits=4)) h: $(a*h*10^-3)km",
  xlabel="km",ylabel="km")

  lines!(ax,[(X(a*h,t,a*b;a=a),Y(a*h,t,b*a;a=a)) for t in tt],color=:red,linestyle=:dash)
  lines!(ax,[(X1(a*h,t,a*b;a=a),Y1(a*h,t,b*a;a=a)) for t in tt],color=:blue,linestyle=:dot)
  bb=b

  let

    ax1=Axis(figure[3,1][1,1:2],title="x-x̂",
    xlabel="θ°",ylabel="m")
    ax2=Axis(figure[3,2][1,1:2],title="1-(h-ĥ)/h",
    xlabel="θ°",ylabel="°",ytickformat="{:.3e}")
    ax3=Axis(figure[4,1][1,1:2],title="n(x)-n̂(x̂)",
    xlabel="θ°",ylabel="mdeg")
    ax4=Axis(figure[4,2][1,1:2],title="n(x)-n(x̂)",
    xlabel="θ°",ylabel="mdeg")
    for b in LinRange(0.995,bb,10)


      lines!(ax1,tt,[
          hypot(X(a*h,t,a*b;a=a)-X1(a*h,t,a*b;a=a),
          Y(a*h,t,a*b;a=a)-Y1(a*h,t,a*b;a=a),
          )
          for t in tt],label="b/a: $(round(b,digits=4))")

        lines!(ax2,tt,[1-
        hypot(X1(a*h,t,a*b;a=a)-X(0,t,a*b;a=a),
        Y1(a*h,t,a*b;a=a)-Y(0,t,a*b;a=a),
        )./(h*a)
        for t in tt],label="b/a: $(round(b,digits=4))")
        lines!(ax3,tt,[
        (
        atand(Y1(a*h,t,a*b;a=a)/(a*b+a*h)^2,X1(a*h,t,a*b;a=a)/(a*(1+h))^2)-
        atand(Y(a*h,t,a*b;a=a)/(a*b)^2,X(a*h,t,a*b;a=a)/a^2)
        )*1e3
        for t in tt],label="b/a: $(round(b,digits=4))")

        lines!(ax4,tt,[
        (
        atand(Y1(a*h,t,a*b;a=a)/(a*b)^2,X1(a*h,t,a*b;a=a)/(a*(1))^2)-
        atand(Y(a*h,t,a*b;a=a)/(a*b)^2,X(a*h,t,a*b;a=a)/a^2)
        )*1e3
        for t in tt],label="b/a: $(round(b,digits=4))")
    end
    Legend(figure[3,1][1,3],ax1)
    Legend(figure[3,2][1,3],ax2)
    Legend(figure[4,1][1,3],ax3)
    Legend(figure[4,2][1,3],ax4)
  end
    figure
  save("$(pwd())/approximation_ellipse.png",figure)
end
