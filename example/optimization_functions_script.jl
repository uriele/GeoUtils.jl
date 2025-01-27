## The new optimization functions
# to be used and tested before the next release

using LinearAlgebra
minoraxis_earth = 0.5 #minor axis of the ellipse <1
const TOLF32=1e-7
const TOLF64=1e-14
@inline atol(::Type{T}) where T<:Float32 =TOLF32
@inline atol(::Type{T}) where T<:Float64 =TOLF64
h=0.1
h2=h*h
#function to minimize

function fg!(g::V1,
  t::T,θ::T,origin::V2,
  direction::V3,b::T,h²::T)::Tuple{T,T}  where {
    V1<:AbstractVector{T},
    V2<:AbstractVector{T},
    V3<:AbstractVector{T}
  } where T<:Real

  cosθ=cos(θ)
  sinθ=sin(θ)
  bsinθ=b*sinθ
  bcosθ=b*cosθ
  b²=b*b
  X=origin[1] + t*direction[1]-cosθ
  Y=origin[2] + t*direction[2]-bsinθ

  AB² = X*X + Y*Y
  distance=AB²-h²
  distance² =distance*distance
  #=
  if distance²<atol(T)
    g.=T(0)
    return h²
  end
  =#
  twice_AB²_minus_h² = 2*(AB²-h²)
  AB_Dᵣ = X*direction[1]+Y*direction[2]
  AB_Tₑ = X*sinθ-Y*bcosθ
  Dᵣ² = direction[1]*direction[1]+direction[2]*direction[2]
  Tₑ² = bcosθ*bcosθ+sinθ*sinθ
  Dᵣ_Tₑ = direction[1]*sinθ-direction[2]*bcosθ
  AB_A=X*cosθ+Y*bsinθ
  eight_AB_Dᵣ_AB_Tₑ=8*AB_Dᵣ*AB_Tₑ

  g1=2*AB_Dᵣ*(twice_AB²_minus_h²)
  g2=2*AB_Tₑ*(1+twice_AB²_minus_h²)

  ∇²f_11= 2*(Dᵣ²)*twice_AB²_minus_h² + 8*(AB_Dᵣ*AB_Dᵣ)
  ∇²f_22= 2*(AB_A+Tₑ²)*twice_AB²_minus_h² + 8*(AB_Tₑ*AB_Tₑ)+2*(AB_A+Tₑ²)
  ∇²f_12= Dᵣ_Tₑ*twice_AB²_minus_h² + eight_AB_Dᵣ_AB_Tₑ
  Denom=∇²f_11*∇²f_22-∇²f_12*∇²f_12
  Denom+= (Denom<atol(T))*atol(T)


  g[1]=( ∇²f_22*g1-∇²f_12*g2)./Denom
  g[2]=(-∇²f_12*g1+∇²f_11*g2)./Denom

  return distance²+AB²,h²+(distance²>atol(T))*distance²

end



function myoptimize!(g,x,origin,direc,minoraxis,h2)
    _old=0.0
    _h2=999.0
    while true
      _fun,_h2=fg!(g,x[1],x[2],origin,direc,minoraxis,h2)
      x.-=0.5.*g
      #@info "x: $x"
      err=abs(_fun-_old)
      _old=_fun
      if err<10^-10
        #@info "err: $err h: $(sqrt(_h2))"
        #@info "b: $(minoraxis), h2: $h2 x: $x"
        break
      end
    end
    return _h2
end





#=
  x
  _h2
x
#using Makie,WGLMakie
x
fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
NN(t)=1/sqrt(cosd(t)^2+minoraxis^2*sind(t)^2)
lines!(ax, [Point2((NN(a)+h)*cosd(a),(minoraxis^2*NN(a)+h)*sind(a)) for a in LinRange(0,361,36001)], color = :blue)
lines!(ax, [Point2((NN(a))*cosd(a),(minoraxis^2*NN(a))*sind(a)) for a in LinRange(0,361,36001)], color = :blue)
scatter!(ax, [Point2(_origin+x[1]*direc) for (x,direc) in zip(eachcol(x),eachcol(direc))], color = :blue)

[
  (scatterlines!(ax, [Point2(_origin+x[1]*direc), Point2(cos(x[2]),minoraxis*sin(x[2]))], color = :black),
lines!(ax, [Point2(_origin), Point2(_origin+1.1*x[1]*direc)] , color = :blue))
  for (x,direc) in zip(eachcol(x),eachcol(direc))]

using Symbolics
@variables t θ r₀[1:2] r₁[1:2] h² b
@variables AB² AB_x AB_y
@variables AB²_minus_h²
@variables AB_Dᵣ AB_Tₑ
@variables tᵢ
@variables AB_A Tₑ²  Dᵣ²  Dᵣ_Tₑ

f1(t,θ)= (r₀[1] + t*r₁[1]-cos(θ))^2 + (r₀[2] + t*r₁[2]-b*sin(θ))^2 |>
  _AB²-> (_AB²-h²)^2
f2(θ)= (r₀[1] +tᵢ *r₁[1] -cos(θ))^2 + (r₀[2] +tᵢ *r₁[2] -b*sin(θ))^2
f(t,θ)= f1(t,θ)+f2(θ)

∇_θ =Differential(θ)
∇_t =Differential(t)

_∇f1= [∇_t(f1(t,θ)); ∇_θ(f1(t,θ))]|> ∂f -> expand_derivatives.(∂f) |> ∂f-> simplify.(∂f)
_∇f2= [∇_t(f2(θ)); ∇_θ(f2(θ))]|> ∂f -> expand_derivatives.(∂f) |> ∂f-> simplify.(∂f)

∇f1= substitute.(_∇f1,(t*r₁[1])=>AB_x-(r₀[1]-cos(θ))) |>
_∂f-> substitute.(_∂f,(t*r₁[2])=>AB_y-(r₀[2]-b*sin(θ))) |>
_∂f-> substitute.(_∂f,(tᵢ*r₁[1])=>AB_x-(r₀[1]-cos(θ))) |>
_∂f-> substitute.(_∂f,(tᵢ*r₁[2])=>AB_y-(r₀[2]-b*sin(θ))) |>
_∂f-> substitute.(_∂f,(-h²)=>AB²_minus_h²-AB_x^2-AB_y^2) |>
_∂f-> substitute.(_∂f,(2*AB_x*r₁[1]+2*AB_y*r₁[2])=>2*AB_Dᵣ) |>
_∂f-> substitute.(_∂f,(2*AB_x*sin(θ)-2*AB_y*b*cos(θ))=> 2*AB_Tₑ)

∇f2= substitute.(_∇f2,(t*r₁[1])=>AB_x-(r₀[1]-cos(θ))) |>
_∂f-> substitute.(_∂f,(t*r₁[2])=>AB_y-(r₀[2]-b*sin(θ))) |>
_∂f-> substitute.(_∂f,(tᵢ*r₁[1])=>AB_x-(r₀[1]-cos(θ))) |>
_∂f-> substitute.(_∂f,(tᵢ*r₁[2])=>AB_y-(r₀[2]-b*sin(θ))) |>
_∂f-> substitute.(_∂f,(-h²)=>AB²_minus_h²-AB_x^2-AB_y^2) |>
_∂f-> substitute.(_∂f,(2*AB_x*r₁[1]+2*AB_y*r₁[2])=>2*AB_Dᵣ) |>
_∂f-> substitute.(_∂f,(2*AB_x*sin(θ)-2*AB_y*b*cos(θ))=> 2*AB_Tₑ)


∇f=∇f1+∇f2

_∇²f1[1]


_∇²f1= [∇_t(_∇f1[1]) ∇_θ(_∇f1[1]); ∇_t(_∇f1[2]) ∇_θ(_∇f1[2])] |> ∂f -> expand_derivatives.(∂f) |> ∂f-> simplify.(∂f)


∇²f1=_∇²f1 |>
_∂f-> substitute.(_∂f,(t*r₁[1])=>AB_x-(r₀[1]-cos(θ))) |>
_∂f-> substitute.(_∂f,(t*r₁[2])=>AB_y-(r₀[2]-b*sin(θ))) |>
_∂f-> substitute.(_∂f,(-h²)=>AB²_minus_h²-AB_x^2-AB_y^2) |>
_∂f-> substitute.(_∂f,(2*AB_x*r₁[1]+2*AB_y*r₁[2])=>2*AB_Dᵣ) |>
_∂f-> substitute.(_∂f,(2*AB_x*sin(θ)-2*AB_y*b*cos(θ))=> 2*AB_Tₑ) |>
_∂f-> substitute.(_∂f,(t*r₁[2])=>AB_y-(r₀[2]-b*sin(θ))) |>
_∂f-> substitute.(_∂f,(tᵢ*r₁[1])=>AB_x-(r₀[1]-cos(θ))) |>
_∂f-> substitute.(_∂f,(tᵢ*r₁[2])=>AB_y-(r₀[2]-b*sin(θ))) |>
_∂f-> substitute.(_∂f,(-h²)=>AB²_minus_h²-AB_x^2-AB_y^2) |>
_∂f-> substitute.(_∂f,(2*AB_x*r₁[1]+2*AB_y*r₁[2])=>2*AB_Dᵣ) |>
_∂f-> substitute.(_∂f,(2*AB_x*cos(θ))=>-2*AB_y*b*sin(θ)+ 2*AB_A) |>
_∂f-> substitute.(_∂f,(2*sin(θ)^2)=> 2*Tₑ²-2*b^2*cos(θ)^2) |>
_∂f-> substitute.(_∂f,(r₁[1]^2+r₁[2]^2)=>Dᵣ²) |>
_∂f-> substitute.(_∂f,(2*r₁[1]*sin(θ)-2*b*r₁[2]*cos(θ))=>Dᵣ_Tₑ)

∇f1


@variables ∂f1_t ∂f1_θ ∂f2_t ∂f2_θ
_∇f2[2]
_∇²f2= [∇_t(_∇f2[1]) ∇_θ(_∇f2[1]); ∇_t(_∇f2[2]) ∇_θ(_∇f2[2])] |> ∂f -> expand_derivatives.(∂f) |> ∂f-> simplify.(∂f)

∇²f2=_∇²f2  |>
_∂f-> substitute.(_∂f,∇f2[2]=>∂f2_θ) |>
_∂f-> substitute.(_∂f,(t*r₁[2])=>AB_y-(r₀[2]-b*sin(θ))) |>
_∂f-> substitute.(_∂f,(tᵢ*r₁[1])=>AB_x-(r₀[1]-cos(θ))) |>
_∂f-> substitute.(_∂f,(tᵢ*r₁[2])=>AB_y-(r₀[2]-b*sin(θ))) |>
_∂f-> substitute.(_∂f,(-h²)=>AB²_minus_h²-AB_x^2-AB_y^2) |>
_∂f-> substitute.(_∂f,(2*AB_x*r₁[1]+2*AB_y*r₁[2])=>2*AB_Dᵣ) |>
_∂f-> substitute.(_∂f,(2*AB_x*cos(θ))=>-2*AB_y*b*sin(θ)+ 2*AB_A) |>
_∂f-> substitute.(_∂f,(2*sin(θ)^2)=> 2*Tₑ²-2*b^2*cos(θ)^2)

∇f2


∇²f= expand.(∇²f1+∇²f2) |> ∂f -> simplify.(∂f)


_c1= Symbolics.coeff.(∇²f,AB²_minus_h²)

_c0= substitute.(∇²f, AB²_minus_h² => 0)

∇f


@variables  c1[1:2,1:2] c0[1:2,1:2]


inv([c1[1,1] c1[1,2];
       c1[2,1] c1[2,2]])


H_=inv([c1[1,1]*AB²_minus_h²+c0[1,1] c1[1,2]*AB²_minus_h²+c0[1,2];
       c1[2,1]*AB²_minus_h²+c0[2,1] c1[2,2]*AB²_minus_h²+c0[2,2]])

       H_[1,1]
       H_[2,2]
       H_[1,1]
       H_[1,1]

H=inv(_c1*AB²_minus_h² + _c0)
D=((c0[1, 1] + AB²_minus_h²*c1[1, 1])*(c0[2, 2] + AB²_minus_h²*c1[2, 2]) - (c0[1, 2] + AB²_minus_h²*c1[1, 2])*(c0[2, 1] + AB²_minus_h²*c1[2, 1]))
@variables dummy
d0=substitute(D,AB²_minus_h²=>0)
d1= D-d0 |> expand |> ff-> Symbolics.coeff.(ff,AB²_minus_h²^2)
d2= D-d0 |> expand |> ff-> substitute(ff,AB²_minus_h²^2=>0) |> ff-> Symbolics.coeff.(ff,AB²_minus_h²) |> simplify


D=((2(AB_A + Tₑ²) + ((4//1)*AB_A + (4//1)*Tₑ²)*AB²_minus_h² + (8//1)*(AB_Tₑ^2))*((8//1)*(AB_Dᵣ^2) + 4AB²_minus_h²*Dᵣ²) - ((8AB_Dᵣ*AB_Tₑ + 2AB²_minus_h²*Dᵣ_Tₑ)^2))


∇²f- ( _c1*AB²_minus_h² + _c0 ) |> ∂f -> expand.(∂f) |> ∂f -> simplify.(∂f)

=#
#=
nu=10000
NUM=nu
#for NUM in nu
bb1=[]
bb2=[]
bb3=[]
using BenchmarkTools
  for NUM in [10,100,1000,10000, 100000, 1000000]
    direc=zeros(2,NUM)
    rand(2,NUM)|> x-> [direc[1:2,i].=-x/hypot(x...) for (i,x) in enumerate(eachcol(x))]
    #direc.= -[1,1]/sqrt(2)
    x=zeros(2,NUM)
    x[1:2,:].=[0.1,pi/3]
    g=zeros(2,NUM)

  @info "eachcol $NUM serial"

  push!(bb1, @benchmark let
    for i in axes(x,2) # zip(eachslice(x,dims=2),eachslice(g,dims=2),eachslice(direc,dims=2))
        @inbounds @fastmath myoptimize!(view($g,1:2,i),view($x,1:2,i),view($_origin,1:2),view($direc,1:2,i),$minoraxis,$h2)
    end
  end
    )


    @info "eachcol $NUM @threads"
    push!(bb2, @benchmark let
      @Threads.threads for i in axes(x,2) # zip(eachslice(x,dims=2),eachslice(g,dims=2),eachslice(direc,dims=2))
          @inbounds @fastmath myoptimize!(view($g,1:2,i),view($x,1:2,i),view($_origin,1:2),view($direc,1:2,i),$minoraxis,$h2)
      end
    end
    )

  @info "eachcol $NUM @batch"
  push!(bb3, @benchmark let
    @batch for i in axes(x,2) # zip(eachslice(x,dims=2),eachslice(g,dims=2),eachslice(direc,dims=2))
        @inbounds @fastmath myoptimize!(view($g,1:2,i),view($x,1:2,i),view($_origin,1:2),view($direc,1:2,i),$minoraxis,$h2)
    end
  end
  )
end
=#
#=
df[!,:speedup_serial].= df[!,:time_old_serial]./df[!,:time_new_serial]
df[!,:speedup_threads].= df[!,:time_old_threads]./df[!,:time_new_threads]
df[!,:speedup_batch].= df[!,:time_old_batch]./df[!,:time_new_batch]


df[!,:gcspeedup_serial].= df[!,:gctime_old_serial]./df[!,:gctime_new_serial]
df[!,:gcspeedup_threads].= df[!,:gctime_old_threads]./max.(1e-15,df[!,:gctime_new_threads])
df[!,:gcspeedup_batch].= df[!,:gctime_old_batch]./max.(1e-15,df[!,:gctime_new_batch])


df[!,:reduced_alloc_serial].= df[!,:alloc_old_serial]./df[!,:alloc_new_serial]
df[!,:reduced_alloc_threads].= df[!,:alloc_old_threads]./max.(1e-15,df[!,:alloc_new_threads])
df[!,:reduced_alloc_batch].= df[!,:alloc_old_batch]./max.(1e-15,df[!,:alloc_new_batch])


df[!,:reduced_memory_serial].=testave.memory[:,1]./_testave.memory[:,1]
df[!,:reduced_memory_threads].=testave.memory[:,2]./_testave.memory[:,2]
df[!,:reduced_memory_batch].=testave.memory[:,3]./_testave.memory[:,3]


df[!,[:speedup_serial,:speedup_threads,:speedup_batch]]
df[!,[:gcspeedup_serial,:gcspeedup_threads,:gcspeedup_batch]]
df[!,[:reduced_alloc_serial,:reduced_alloc_threads,:reduced_alloc_batch]]
df[!,[:reduced_memory_serial,:reduced_memory_threads,:reduced_memory_batch]]
=#

using WGLMakie,Makie



fig=Figure()
ax=Axis(fig[1,1])
lines!(ax,[Point(cos(t),minoraxis_earth*sin(t)) for t in LinRange(0,2π,1000)],color=:blue)
_NN(t)=1/sqrt(cos(t)^2+sin(t)^2*minoraxis_earth^2)

h=sqrt(h2)
lines!(ax,[Point((_NN(t)+h).*cos(t),(minoraxis_earth^2*_NN(t)+h)*sin(t)) for t in LinRange(0,2π,1000)],color=:blue)

for i in 1:10:1000 #axes(x,2)
  lines!(ax,[Point(_origin[1,i]+_direc[1,i]*t,
  _origin[2,i]+_direc[2,i]*t) for t in [0,4*x[1,i]]],color=:red);
  scatter!(ax,[Point(_origin[1,i]+_direc[1,i]*x[1,i],
  _origin[2,i]+_direc[2,i]*x[1,i]),
  Point(cos(x[2,i]),minoraxis_earth*sin(x[2,i]))])


end
fig

T=Float64
N=100*10*100
_origin=Array{T,2}(undef,2,100*10)
_direc =Array{T,2}(undef,2,100*10)
tpos=LinRange(0,2π,100)
x=Array{T,2}(undef,2,100*10)
count=0
for i in 1:100
  for j in 1:10
    count+=1
    _origin[1,count]=1.2*cos(tpos[i])
    _origin[2,count]=1.2*sin(tpos[i])
    tx=-_origin[2,count]
    ty=_origin[1,count]
    dx=tx*cosd(j/10-27)+ty*sind(j/10-27)
    dy=-tx*sind(j/10-27)+ty*cosd(j/10-27)
    magn=1/hypot(dx,dy)
    _direc[1,count]=dx*magn
    _direc[2,count]=dy*magn
    x[1,count]=0.0
    x[2,count]=tpos[i]
  end
end
hh=Array{Float64,1}(undef,100*10)
g1=similar(x)
using BenchmarkTools
@benchmark for i in axes($x,2) # zip(eachslice(x,dims=2),eachslice(g,dims=2),eachslice(direc,dims=2))
  @inbounds @fastmath $hh[i]=myoptimize!(view($g1,:,i),
    view($x,:,i),
    view($_origin,:,i),
    view($_direc,:,i),$minoraxis_earth,$h2)
end
-



hum=LinRange(0,0.9,1000)

lines(hum,[refractive_index(;humidity=h,wavelength=10.0,CO2ppm=0.0) for h in hum],color=:blue)
