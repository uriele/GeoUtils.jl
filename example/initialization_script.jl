
using GeoUtils
using StructArrays
using CoordRefSystems
using StaticArrays
using Base:IEEEFloat
using LinearAlgebra: ⋅
using Unitful: ustrip
#using Makie,WGLMakie
const MODEL=Ref{AirModel}(Carlotti())
setModel(model::AM=Carlotti()) where AM<:AirModel= MODEL[]=model
setModel(::Type{AM}) where AM<:AirModel= MODEL[]=AM()
getModel()=MODEL[]
setModel()
getModel()
const INTERPOLATION=Ref{AbstractPressureInterpolation}(LinearPressure())
setPressureInterpolation(interpolation::PI=LinearPressure()) where PI<:AbstractPressureInterpolation= INTERPOLATION[]=interpolation
setPressureInterpolation(::Type{PI}) where PI<:AbstractPressureInterpolation= INTERPOLATION[]=PI()
getPressureInterpolation()=INTERPOLATION[]
setPressureInterpolation(LogarithmicPressure)
setPressureInterpolation()


function initialize_raytracing_plot(h_levels,θ_radii)
  fig=Figure(size=(600,800))
  ax=Axis(fig[1:2,1],
  xlabel="w",
  ylabel="z"
  )

    ax2a=Axis(fig[3,1][1,1],
    xlabel="Altitude (km)",
    ylabel="Refractive Index"
    )
    ax2b=Axis(fig[3,1][1,2],
    xlabel="Iteraction",
    ylabel="Altitude (km)",
    )
    #scatter!(ax,hhh.*a)
    #scatter!(ax,allintersections)
    for h in h_levels
      lines!(ax,[let
        convert(ECEF2D{NormalizedEarth},LLA2D{NormalizedEarth}(h,tt)) |>
        x-> (x.w,x.z)
        end
        for tt in LinRange(0,361,1000)])
    end
    for tt in θ_radii
      lines!(ax,[let
        convert(ECEF2D{NormalizedEarth},LLA2D{NormalizedEarth}(h,tt)) |>
        x-> (x.w,x.z)
        end
        for h in h_levels])
    end
  return (fig,ax,ax2a,ax2b)
end

########
# Read atmosphere
using BenchmarkTools

atmosphere=read_local_atmosphere("./data_atm/INP_FILES");


 θᵢ=" 0.00      1.00      2.00      3.00      4.00      5.00      6.00
7.00      8.00      9.00     10.00     11.00     12.00     13.00
14.00     15.00     16.00     17.00     18.00     19.00     20.00
21.00     22.00     23.00     24.00     25.00     26.00     27.00
28.00     29.00     30.00     31.00     32.00     33.00     34.00
35.00     36.00     37.00     38.00     39.00     40.00     41.00
42.00     42.50     43.00     43.50     44.00     44.50     45.00
45.25     45.50     45.75     46.00     46.10     46.20     46.30
46.40     46.50     46.60     46.70     46.80     46.90     47.00
47.10     47.20     47.30     47.40     47.50     47.60     47.70
47.80     47.90     48.00     48.10     48.20     48.30     48.40
48.50     48.60     48.70     48.80     48.90     49.00     49.10
49.20     49.30     49.40     49.50     49.60     49.70     49.80
49.90     50.00     50.10     50.20     50.30     50.40     50.50
50.60     50.70     50.80     50.90     51.00     51.10     51.20
51.30     51.40     51.50     51.60     51.70     51.80     51.90
52.00     52.10     52.20     52.30     52.40     52.50     52.60
52.70     52.80     52.90     53.00     53.10     53.20     53.30
53.40     53.50     53.60     53.70     53.80     53.90     54.00
54.25     54.50     54.75     55.00     55.25     55.50     56.00
56.50     57.00     57.50     58.00     59.00     60.00     61.00
62.00     63.00     64.00     65.00     66.00     67.00     68.00
69.00     70.00     71.00     72.00     73.00     74.00     75.00
76.00     77.00     78.00     79.00     80.00     81.00     82.00
83.00     84.00     85.00     86.00     87.00     88.00     89.00
90.00     91.00     92.00     93.00     94.00     95.00     96.00
97.00     98.00     99.00    100.00    101.00    102.00    103.00
104.00    105.00    106.00    107.00    108.00    109.00    110.00
111.00    112.00    113.00    114.00    115.00    116.00    117.00
118.00    119.00    120.00    121.00    122.00    123.00    124.00
125.00    126.00    127.00    128.00    129.00    130.00    131.00
132.00    133.00    134.00    135.00    136.00    137.00    138.00
139.00    140.00    141.00    142.00    143.00    144.00    145.00
146.00    147.00    148.00    149.00    150.00    151.00    152.00
153.00    154.00    155.00    156.00    157.00    158.00    159.00
160.00    161.00    162.00    163.00    164.00    165.00    166.00
167.00    168.00    169.00    170.00    171.00    172.00    173.00
174.00    175.00    176.00    177.00    178.00    179.00    180.00
181.00    182.00    183.00    184.00    185.00    186.00    187.00
188.00    189.00    190.00    191.00    192.00    193.00    194.00
195.00    196.00    197.00    198.00    199.00    200.00    201.00
202.00    203.00    204.00    205.00    206.00    207.00    208.00
209.00    210.00    211.00    212.00    213.00    214.00    215.00
216.00    217.00    218.00    219.00    220.00    221.00    222.00
223.00    224.00    225.00    226.00    227.00    228.00    229.00
230.00    231.00    232.00    233.00    234.00    235.00    236.00
237.00    238.00    239.00    240.00    241.00    242.00    243.00
244.00    245.00    246.00    247.00    248.00    249.00    250.00
251.00    252.00    253.00    254.00    255.00    256.00    257.00
258.00    259.00    260.00    261.00    262.00    263.00    264.00
265.00    266.00    267.00    268.00    269.00    270.00    271.00
272.00    273.00    274.00    275.00    276.00    277.00    278.00
279.00    280.00    281.00    282.00    283.00    284.00    285.00
286.00    287.00    288.00    289.00    290.00    291.00    292.00
293.00    294.00    295.00    296.00    297.00    298.00    299.00
300.00    301.00    302.00    303.00    304.00    305.00    306.00
307.00    308.00    309.00    310.00    311.00    312.00    313.00
314.00    315.00    316.00    317.00    318.00    319.00    320.00
321.00    322.00    323.00    324.00    325.00    326.00    327.00
328.00    329.00    330.00    331.00    332.00    333.00    334.00
335.00    336.00    337.00    338.00    339.00    340.00    341.00
342.00    343.00    344.00    345.00    346.00    347.00    348.00
349.00    350.00    351.00    352.00    353.00    354.00    355.00
356.00    357.00    358.00    359.00" |> x-> split(x,"\n") |> x-> join(x," ") |>
x-> split(x," ") |> x-> filter(y-> y≠"",x) |>
x-> [parse(Float64,y) for y in x]; #.+90.0



# normalizetion factor for h_knots

majoraxis_earth=majoraxis(ellipsoid(WGS84Latest)) |> x-> uconvert(km,x) |> ustrip;

squared_eccentricity_earth=eccentricity²(ellipsoid(NormalizedEarth));
reduced_minoraxis_earth=minoraxis(ellipsoid(NormalizedEarth));


hᵢ="0.000000000000000E+000   3.60000000000000        3.80000000000000
4.00000000000000        4.20000000000000        4.40000000000000
4.60000000000000        4.80000000000000        5.00000000000000
5.20000000000000        5.40000000000000        5.60000000000000
5.80000000000000        6.00000000000000        6.20000000000000
6.40000000000000        6.60000000000000        6.80000000000000
7.00000000000000        7.20000000000000        7.40000000000000
7.60000000000000        7.80000000000000        8.00000000000000
8.20000000000000        8.40000000000000        8.60000000000000
8.80000000000000        9.00000000000000        9.20000000000000
9.40000000000000        9.60000000000000        9.80000000000000
10.0000000000000        10.2000000000000        10.4000000000000
10.6000000000000        10.8000000000000        11.0000000000000
11.2000000000000        11.4000000000000        11.6000000000000
11.8000000000000        12.0000000000000        12.2000000000000
12.4000000000000        12.6000000000000        12.8000000000000
13.0000000000000        13.2000000000000        13.4000000000000
13.6000000000000        13.8000000000000        14.0000000000000
14.2000000000000        14.4000000000000        14.6000000000000
14.8000000000000        15.0000000000000        15.2000000000000
15.4000000000000        15.6000000000000        15.8000000000000
16.0000000000000        16.2000000000000        16.4000000000000
16.6000000000000        16.8000000000000        17.0000000000000
17.2000000000000        17.4000000000000        17.6000000000000
17.8000000000000        18.0000000000000        18.2000000000000
18.4000000000000        18.6000000000000        18.8000000000000
19.0000000000000        19.2000000000000        19.4000000000000
19.6000000000000        19.8000000000000        20.0000000000000
21.0000000000000        22.0000000000000        23.0000000000000
24.0000000000000        25.0000000000000        29.0000000000000
37.0000000000000        44.0000000000000        68.0000000000000
68.1500000000000    "  |> x-> split(x,"\n") |> x-> join(x," ") |>
x-> split(x," ") |> x-> filter(y-> y≠"",x) |>
x-> [parse(Float64,y) for y in x]./majoraxis_earth;
# Generate discretized atmosphere using the knots
pressure,temperature,refractive=discretize_atmosphere(atmosphere,hᵢ,θᵢ; model=MODEL[],interpolation_pressure=INTERPOLATION[]);


############################################
# get coordinates in LLA and ECEF
#########################################
begin
  idx270=findfirst(θᵢ.==90.0+180);
  θᵢ=SemiCircularVector(θᵢ) |> x->  x[idx270:end+idx270-1];

  θᵢ.+=90.0;
  θᵢ.=mod.(θᵢ,360.0);
  pressure.=pressure[:,idx270:end+idx270-1];
  temperature.=temperature[:,idx270:end+idx270-1];
  refractive.=refractive[:,idx270:end+idx270-1];
end


##########################################
# Used reversed in discretize_atmosphere
sort!(hᵢ;rev=true);
##########################################

# I will not need it in the new code, only the angle and the height is required
(
  scale_levels2,
  line_radii2
)=getIntersectionObjects(NormalizedEarth(),hᵢ,θᵢ);

# although I will need a new structure for the angles
θᵢ_effective= GeoUtils._surface_from_geodesic_θrad_to_geocentric_θrad.(deg2rad.(θᵢ));


# create rays
#@benchmark orbit=read_orbit("./data_atm/INP_FILES/orbit.dat")
#orbit1=GeoUtils.__read_orbit("./data_atm/INP_FILES/orbit.dat");
#normalize_datum!(orbit1)
#|> # read the orbit file
transposed_orbit=read_orbit("./data_atm/INP_FILES/orbit.dat")

orbit=StructArray{SatOrbit{Float64},1}(undef,size(transposed_orbit,2),size(transposed_orbit,1))

for (o1,to1) in zip(StructArrays.components(orbit),StructArrays.components(transposed_orbit))
  @inbounds o1.=transpose(to1);
end

normalize_datum!(orbit);


rr2=StructArray{Ray2D{Float64}}(undef,length(orbit));
rays_from_orbit!(rr2,orbit)


x0=fill(0.0,(2,length(rr2)));
g= similar(x0)
slack=fill(0.0,(3,length(rr2)));
ii=1

hᵢ_max=maximum(hᵢ)
hᵢ_max²=hᵢ_max*hᵢ_max
θmin=atan(rr2.origin[ii][2],rr2.origin[ii][1])

#for ii in 1:1 #eachindex(rr2)


#  θmin=x0[2,ii]+0.2;#θminatan(rr2[ii].origin[2],rr2[ii].origin[1])
  θmax=θmin+10/180*pi
  #@info "i: $ii"
  x0[2,ii]


  GeoUtils.optimization_augmented_lagrangian!(view(g,:,ii),view(x0,:,ii),rr2.origin[ii],rr2.direction[ii],reduced_minoraxis_earth,hᵢ_max²,θmin,θmax,view(slack,:,ii));
  slack
  θmin
  #@info "=============="
  rr2.origin[ii]=rr2.origin[ii]+rr2.direction[ii]*x0[1,ii]
end

using Makie,WGLMakie

perp=StructArray{Ray2D{Float64}}(undef,length(rr2))
for ii in eachindex(rr2)
  perp.origin[ii]=Vec2(cos(x0[2,ii]),reduced_minoraxis_earth*sin(x0[2,ii]))
end

v=rr2.origin-perp.origin
hs=[hypot(x[1],x[2]) for x in v]
hs[1]-hᵢ_max
diff([extrema(hs)...])



fig=Figure()
ax=Axis(fig[1,1])
scatter!(ax,rr2.origin)
scatter!(ax,perp.origin)
scatter!(ax,Point2(orbit.w[1],orbit.z[1]))


[]

(v[1]⋅[sin(x0[2,1])/reduced_minoraxis_earth;-cos(x0[2,1])])

x0


# Find intersection with outmost level
function find_input_interception_new_lagrangian2!(index,x,rays::Ar, # ray
  n₀,refractive,h_levels,θ_radii) where Ar<:AbstractArray{Ray2D{T}} where T
  h_levels_max=maximum(h_levels)
  @info "h_levels_max=$h_levels_max"
  h_levels_max²=h_levels_max*h_levels_max
  @info "h_levels_max²=$h_levels_max²"
  ii=findfirst(h_levels.==h_levels_max)
  minoraxis_earth=minoraxis(ellipsoid(NormalizedEarth))
  invb2=1/(minoraxis_earth*minoraxis_earth)
  θ_max=maximum(θ_radii)

   @inbounds for i in eachindex(rays)

    x[2,i]=atan(rays[i].origin[2],rays[i].origin[1])
    h²=GeoUtils.optimize_distance_h²_from_surface_lagrangian2!(view(x,:,i),
      rays.origin[i],rays.direction[i],minoraxis_earth,h_levels_max²,
      T(0),T(2pi))

    rays.origin[i]=rays[i](x[1,i])
    x[2,i]=mod(x[2,i],2pi)
    jj=ifelse(x[2,1]>θ_max,length(θ_radii), findlast(θ_radii.<= x[2,i]))


    index[1,i]=ii
    index[2,i]=jj
    n₁=refractive[ii,jj]

    ############
    n₀₁=n₀/n₁
    n₀₁²=n₀₁*n₀₁
    ############
    _hypot=hypot(rays[i].origin[1],rays[i].origin[2]*invb2)
    Normal_x=rays[i].origin[1]/_hypot
    Normal_y=rays[i].origin[2]*invb2/_hypot
    cosθₙ=-(Normal_x*rays[i].direction[1]+Normal_y*rays[i].direction[2])
    sinθₙ² =n₀₁²*(1-cosθₙ*cosθₙ)
    if sinθₙ²≤1
      tmp=(n₀₁*cosθₙ-sqrt(1-sinθₙ²))
      rays.direction[i]=Vec2(n₀₁*rays[i].direction[1]+tmp*Normal_x,n₀₁*rays[i].direction[1]+tmp*Normal_y)
    else
      rays.direction[i]=Vec2(rays[i].direction[1]-2*cosθₙ*Normal_x,rays[i].direction[2]-2*cosθₙ*Normal_y)
    end
  end
end

function _first_input(index::A1,x::A2,_origin::B1,_direction::B2,n₀::T,b::T,invb2::T,hlevel²::T,refractive,θ_radii,θ_max) where {A1<:AbstractArray,A2<:AbstractArray,B1<:AbstractArray{V},B2<:AbstractArray{V}} where V<:Vec2{T} where T
    x[2]=atan(_origin[1][2],_origin[1][1])
    h²=GeoUtils.optimize_distance_h²_from_surface_lagrangian2!(x,
      _origin[1],_direction[1],b,h_levels_max²,
      x[2],T(4pi))

    _origin[1]=_origin[1]+_direction[1]*x[1]
    x[2]=mod(x[2],2pi)
    jj=ifelse(x[2]>θ_max,length(θ_radii), searchsortedlast(θ_radii, x[2]))

    index[1]=1
    index[2]=jj
    n₁=refractive[jj]


    ############
    n₀₁=n₀/n₁
    n₀₁²=n₀₁*n₀₁
    ############
    _hypot=hypot(_origin[1][1],_origin[1][2]*invb2)
      Normal_x=_origin[1][1]/_hypot
      Normal_y=_origin[1][2]*invb2/_hypot
    # It's a radius
    if (h²-hlevel²)> 10^-10
      (Normal_x,Normal_y)=(Normal_y,-Normal_x)
    end
    _bendme!(_direction,Normal_x,Normal_y,n₀₁,n₀₁²)
    return h²
end

const ACCEPTED_ERROR=10^-10

function _next_input(index::A1,x::A2,_origin::B1,_direction::B2,n₀::T,hlevel²::T,refractive_radius::T,
  refractive_level::T,θ_prev::T,θ_next::T,b::T,invb2::T,isdescending::Bool)::Tuple{T,Bool} where {A1<:AbstractArray,A2<:AbstractArray,
  B1<:AbstractArray{V},B2<:AbstractArray{V}} where V<:Vec2{T} where T<:IEEEFloat

  # invert theta if looping
  (θ_min,θ_max)=(θ_prev>θ_next) ? (θ_next,θ_prev) : (θ_prev,θ_next)

  h²=GeoUtils.optimize_distance_h²_from_surface_lagrangian2!(x,
    _origin[1],_direction[1],b,hlevel²,
    θ_min,θ_max)

  _origin[1]=_origin[1]+_direction[1]*x[1]
  Δh²=abs(h²-hlevel²)

  n₁=n₀
  # if true becomes false
  isdescending=!isdescending #flip
#  @info "Δh² = $Δh² h²=$h² hlevel²=$hlevel², x=$x"
  if Δh²<=GeoUtils.atol(T)  # if it hits the level
    @info "get level"
    isdescending=!isdescending #flop
    index[1]+=(isdescending)-!(isdescending)
    n₁=refractive_level;
  elseif abs(x[2]-θ_next)<=GeoUtils.atol(T) # if it hits the next radius
    @info "get radius"
    isdescending=!isdescending #flop
    index[2]+=index[2]+1
    n₁=refractive_radius;
  end

  ############
  n₀₁=n₀/n₁
  n₀₁²=n₀₁*n₀₁
  ############
  _hypot=hypot(_origin[1][1],_origin[1][2]*invb2)
    Normal_x=_origin[1][1]/_hypot
    Normal_y=_origin[1][2]*invb2/_hypot
  # It's a radius
  if (h²-hlevel²)> 10^-10
    (Normal_x,Normal_y)=(Normal_y,-Normal_x)
  end
  _bendme!(_direction,Normal_x,Normal_y,n₀₁,n₀₁²)
  return h²,isdescending
end

function _bendme!(_direction::A,Normal_x::T,Normal_y::T,n₀₁::T,n₀₁²::T) where A<:AbstractArray{V} where V<:Vec2{T} where T
  cosθₙ=-(Normal_x*_direction[1][1]+Normal_y*_direction[1][2])
  sinθₙ² =n₀₁²*(1-cosθₙ*cosθₙ)
  if sinθₙ²≤1
    tmp=(n₀₁*cosθₙ-sqrt(1-sinθₙ²))
    dx=n₀₁*_direction[1][1]+tmp*Normal_x
    dy=n₀₁*_direction[1][2]+tmp*Normal_y
    invmagnitude=1/hypot(dx,dy)
    _direction[1]=Vec2(dx*invmagnitude,dy*invmagnitude)
  else
    dx=_direction[1][1]-2*cosθₙ*Normal_x
    dy=_direction[1][2]-2*cosθₙ*Normal_y
    invmagnitude=1/hypot(dx,dy)
    _direction[1]=Vec2(dx,dy)
  end
end
rays=rays[:]
rays2=deepcopy(rays)
index=zeros(Int,2,length(rays))
x0=zeros(Float64,2,length(rays))
b=minoraxis(ellipsoid(NormalizedEarth))
invb2=1/(b*b)
h_levels_max=maximum(hᵢ)
h_levels_max²=h_levels_max*h_levels_max
θ_max=maximum(θᵢ_effective)
#ProfileCanvas.@profview_allocs

hᵢ²=hᵢ.*hᵢ
qt=ones(Float64,length(rays2[1:3]))
nlevels=size(refractive,1)
for i in 1 #eachindex(rays2[1])
  isfirst=true
  isdescending=true
  isinside=true
  count=0
  while isinside
    if isfirst
      h²=_first_input(view(index,:,i),view(x0,:,i),view(rays2.origin,i),
      view(rays2.direction,i),1.0,b,invb2,h_levels_max²,view(refractive,1,:),θᵢ_effective,θ_max)
      isfirst=false
      qt[i]=min(h²,qt[i])
      θ_next=θᵢ_effective[index[2,i]+1]
      θ_prev=θᵢ_effective[index[2,i]]

      #@info "θ_next=$θ_next θ_prev=$θ_prev θ:$(x0[2,i])"
    else
     # @info "itr: $count x0: $(x0[:,i])"
      next_level=isdescending-!isdescending
      #precompute current n
      n₀=refractive[index[:,i]...]
      #precompute n of next radius
      n_next_radius=refractive[index[1,i],index[2,i]+1]
      #precompute n of next level
      #@info " count:$count index=$(index[:,i]) next_level=$next_level x0=$(x0[:,i])"
      n_next_level=refractive[index[1,i]+next_level,index[2,i]]
      #precompute h² of next level
      hlevel²=hᵢ²[index[1,i]+isdescending]
      #precompute θ of next level
      θ_next=θᵢ_effective[index[2,i]+1]
      θ_prev=θᵢ_effective[index[2,i]]
      #@info "θ_next=$θ_next θ_prev=$θ_prev hlevel²=$hlevel² isdescending:$isdescending"
      #@info " t:$(x0[:,i]) θ:$(x0[2,i])"
      #@info "hold=$(hᵢ²[index[1,i]]) hnew=$(hᵢ²[index[1,i]+isdescending])"

      (h²,isdescending)=_next_input(view(index,:,i),view(x0,:,i),view(rays2.origin,i),
      view(rays2.direction,i),n₀,hlevel²,
      n_next_radius,n_next_level,
      θ_prev,θ_next,
      b,invb2,
      isdescending)

    #  @info "θ_next=$θ_next θ_prev=$θ_prev θ:$(x0[2,i])"

      qt[i]=min(h²,qt[i])
    end
    @info "itr: $count x0: $(x0[:,i]), qt: $(sqrt(qt[i])*majoraxis_earth)km"
    isinside= 1<=index[1,i]<nlevels
    count+=1
    isinside=count<10
  end
end

using Makie,WGLMakie
fig=Figure()
ax=Axis(fig[1,1])
_NN(t,b)=1/sqrt(cos(t)^2+b^2*sin(t)^2)
tt=LinRange(0,2pi,1000)
h1=hᵢ[1]
h2=hᵢ[2]
b
lines!(ax,[Point2f((_NN(t,b)+h1)*cos(t),(b^2*_NN(t,b)+h1)*sin(t)) for t in tt])
lines!(ax,[Point2f((_NN(t,b)+h2)*cos(t),(b^2*_NN(t,b)+h2)*sin(t)) for t in tt])
direc1=θ_next |> x-> [cos(x),sin(x)/b] |> x-> x./hypot(x[1],x[2])
lines!(ax,[Point2f(cos(θᵢ_effective[index[2,1]])+c*direc1[1],b*sin(θᵢ_effective[index[2,1]])+c*direc1[2]) for c in [0,h1]])

direc2=θ_prev |> x-> [cos(x),sin(x)/b] |> x-> x./hypot(x[1],x[2])
lines!(ax,[Point2f(cos(θᵢ_effective[index[2,1]+1])+c*direc2[1],b*sin(θᵢ_effective[index[2,1]+1])+c*direc2[2]) for c in [0,h1]])
scatter!(ax,[Point2f(rays2.origin[1][1],rays2.origin[1][2])])
scatter!(ax,[Point2f(rays3.origin[1][1],rays3.origin[1][2])])

rays3=deepcopy(rays2)
x1=deepcopy(x0)
index1=deepcopy(index)
isdescending=true
next_level=isdescending-!isdescending

n₀=refractive[index[:,1]...]

n_next_radius=refractive[index[1,1],index[2,1]+1]

n_next_level=refractive[index[1,1]+next_level,index[2,1]]
      #precompute h² of next level
hcurrent2=hᵢ²[index[1,1]]

hlevel²=hᵢ²[index[1,1]+isdescending]
      #precompute θ of next level
θ_next=θᵢ_effective[index[2,1]+1]
θ_prev=θᵢ_effective[index[2,1]]
x1[2]

(h²,isdescending)=_next_input(view(index1,:,1),view(x1,:,1),view(rays3.origin,1),
      view(rays3.direction,1),n₀,hlevel²,
      n_next_radius,n_next_level,
      θ_prev,θ_next,
      b,invb2,
      isdescending)


# Find intersection with outmost level
  n₀;  #refractive index ray
  refractive_index_map=refractive, # map of refractive index
  h_levels=h_levels, # levels for i
  θ_radii=θ_radii) where Ar<:AbstractArray{Ray2D{T}} where T
  h_levels_max=maximum(h_levels)
  h_levels_max²=h_levels_max*h_levels_max
  ii=findfirst(h_levels.==h_levels_max)
  minoraxis_earth=minoraxis(ellipsoid(NormalizedEarth))
  θ_max=maximum(θ_radii)

  Normal=[0.0,0.0]
  for i in eachindex(rays)

    h²=GeoUtils.optimize_distance_h²_from_surface_lagrangian!(view(x,:,i),
      rays.origin[i],rays.direction[i],minoraxis_earth,h_levels_max²,
      T(0),T(2pi))

    rays.origin[i]=rays[i](x[1,i])
    θ=mod(x[2,i],2pi)
    jj=ifelse(θ>θ_max,length(θ_radii), findlast(θ_radii.<= θ))


    index[1,i]=ii
    index[2,i]=jj
    n₁=refractive_index_map[ii,jj]

    ############
    n₀₁=n₀/n₁
    n₀₁²=n₀₁*n₀₁
    ############
    Normal.=GeoUtils._normal_vector_new(rays.origin[i],minoraxis_earth)
    cosθₙ=-Normal⋅rays.direction[i]
    sinθₙ² =n₀₁²*(1-cosθₙ*cosθₙ)
    if sinθₙ²≤1
      rays.direction[i]=n₀₁*rays.direction[i]+(n₀₁*cosθₙ-sqrt(1-sinθₙ²))*Normal
    else
      rays.direction[i]=rays.direction[i]-2*cosθₙ*Normal
    end
  end
end



# Find intersection with outmost level
function find_first_intersection_ellipse(ray::Ray2D, # ray
  n₀;  #refractive index ray
  refractive_index_map=refractive, # map of refractive index
  h_levels=h_levels, # levels for i
  θ_radii=θ_radii, # radii angles for j
  scale_levels=scale_levels) # scaling mapping to unit circle
  h_levels_max=maximum(h_levels)
  # index of outmost level (to be sure it is sorted in the right order)
  ii=findfirst(h_levels.==h_levels_max)
   # find intersection point
  t_intersection=advance(BottomIntersection(),ray,scale_levels[ii])

  # find j on ellipse
  θ,h=geocentric_xy_to_geodesic_θ(ray(t_intersection+GeoUtils.SHIFTORIGIN)...)

  #θ=ray(t_intersection) |> x-> atand(x[2],x[1]) |> x-> mod(x,360)
  #@info "θ=$θ h=$h h=$(h_levels_max)"
  θ_max=maximum(θ_radii)

  # Take advantage of the semicircular matrix, if the theta is greater
  # that 359.0, we loop around to find the leftmost value
  jj=ifelse(θ>θ_max,length(θ_radii), findlast(θ_radii.<= θ))
  #@info "jj=$jj, θ_max=$θ_max,θ=$θ"

  n₁=refractive_index_map[ii,jj]
  n₀= ifelse(n₀>0,n₀,n₁)
  #@info "n₀=$n₀, n₁=$n₁"
  # compute Interface
  return let # return the initial ray
    ray_out,target=bend(BottomIntersection(),ray,t_intersection;
      n₀=n₀, # refractive index space
      n₁=n₁, # refractive index atmosphere
      h=h_levels[ii]) # level
    (ray_out,target,h_levels[ii], # return the level indexndex
    (ii,jj),t_intersection) # return the intersection i
    end
end
#=
@kwdef mutable struct Register{T}
  flag::Bool=false
  iteration::Int=0
  h::T=T(NaN)
  θ::T=T(NaN)
  w::T=T(NaN)
  z::T=T(NaN)
end
Register(flag,w::T,z::T) where T= let
  lla=convert(LLA2D{NormalizedEarth},ECEF2D{NormalizedEarth}(w,z)) |> x-> (x.h,x.θ)
  Register{T}(flag,0,w,z,lla...)
end
next!(register::Register)= register.iteration+=1
setTangent!(register::Register,w,z)= begin
  register.w=w
  register.z=z
  register.h=convert(LLA2D{NormalizedEarth},ECEF2D{NormalizedEarth}(w,z)) |> x-> x.h
  register.θ=convert(LLA2D{NormalizedEarth},ECEF2D{NormalizedEarth}(w,z)) |> x-> x.θ
  register.flag=true
end
setTangent(register::Register,v::AbstractVector{T}) where T = setTangent!(register,v...)


settangentquote!(register::Register,w,z,count)= begin
  register.w=w
  register.z=z
  register.h=convert(LLA2D{NormalizedEarth},ECEF2D{NormalizedEarth}(w,z)) |> x-> x.h
  register.θ=atand(z/w) |> x-> mod(x-90,360)
  register.iteration=count
  register.flag=true
end

Base.show(io::IO,register::Register)= begin
  if register.flag
    println(io,"Tangent Point Found at $(register.iteration) with h: $(register.h), θ: $(register.θ)")
    println(io,"w: $(register.w), z: $(register.z)")
  else
    println(io,"Tangent Point Not Found  ($(register.iteration))")
  end
end

#==========#
_angleme(r::Ray2D) = r(0)-r(1) |> x->(x[2],x[1])  |> x-> atand(x...)
_angleme(r1::Ray2D,r::Ray2D) = r1(0)-r(0) |> x->(x[2],x[1])  |> x-> atand(x...)
##############


Register()=Register{Float64}()
=#
