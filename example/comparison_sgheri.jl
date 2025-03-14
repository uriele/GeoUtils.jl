using NCDatasets
using DataFrames
using Makie,WGLMakie
# Open the file
nc=Dataset("Output_CAIRT/cairt_trace_out.nc")
orbit_nc=nc.group["orbit"]
atm_nc=nc.group["atm"]
scans_nc=nc.group["scans"]
integrals_nc=nc.group["integrals"]

using GeoUtils

using Unitful:°,hPa
const MODEL=Ref{AirModel}(Carlotti())
const INTERPOLATION=Ref{AbstractPressureInterpolation}(LinearPressure())

setModel(A::AirModel)= MODEL[]=A
getModel()=MODEL[]
setModel()=MODEL[]=(Ciddor())
setPressureInterpolation(A::AbstractPressureInterpolation)= INTERPOLATION[]=A
setPressureInterpolation()= INTERPOLATION[]=(LinearPressure())
getPressureInterpolation()=INTERPOLATION[]
getPressureInterpolation()
setModel()

# only radii
#setDebugIntersection(0)
include("initialization_script.jl")

#get discretized atmosphere
z,phi= atm_nc["z"][:],atm_nc["phi"][:]
pressure,temperature= atm_nc["pressure"][:,:],atm_nc["temperature"][:,:]

N,M=length(z),length(phi)

n_index=Matrix{Union{Missing,Float64}}(undef,N-1,M-1)
pressure_radial=(pressure[:,2:end]+pressure[:,1:end-1])./2
temperature_radial=(temperature[:,2:end]+temperature[:,1:end-1])./2

temperature_ave=(temperature_radial[1:end-1,:]+temperature_radial[2:end,:])./2
# Compute average pressure longitudinally
diff_press=(pressure_radial[1:end-1,:]-pressure_radial[2:end,:])
diff_log_press=log.(pressure_radial[1:end-1,:])-log.(pressure_radial[2:end,:])

pressure_ave=diff_press./diff_log_press

# get the first scans
scans_nc["engineering_quotes"][:]

orbit_nc
scans_nc["scans_phi"][:]
integrals_nc["phi_idx"][:]==integrals_nc["z_idx"][:]

scans_nc
[(phi,z) for (phi,z) in zip(integrals_nc["phi_idx"][:,1,1],integrals_nc["z_idx"][:,1,1])]

pressure_linear,temperature_,refractive_linear_carlotti,_,_=discretize_atmosphere(atmosphere,h_knots,θ_knots; model=Carlotti(),interpolation_pressure=INTERPOLATION[]);
_,_,refractive_linear_ciddor,_,_=discretize_atmosphere(atmosphere,h_knots,θ_knots; model=Ciddor(),interpolation_pressure=INTERPOLATION[]);
_,_,refractive_linear_mathar,_,_=discretize_atmosphere(atmosphere,h_knots,θ_knots; model=Mathar4(),interpolation_pressure=INTERPOLATION[]);
pressure_log,_,refractive_log_carlotti,_,_=discretize_atmosphere(atmosphere,h_knots,θ_knots; model=Carlotti(),interpolation_pressure=LogarithmicPressure());
_,_,refractive_log_ciddor,_,_=discretize_atmosphere(atmosphere,h_knots,θ_knots; model=Ciddor(),interpolation_pressure=LogarithmicPressure());
_,_,refractive_log_mathar,_,_=discretize_atmosphere(atmosphere,h_knots,θ_knots; model=Mathar4(),interpolation_pressure=LogarithmicPressure());

figure=Figure(size=(1300,800))

ax_11=Axis(figure[1,1][1,1],yticks=(0:90:360, ["0°","90°","180°","270°","360°"]))
ax_12=Axis(figure[1,2][1,1],yticks=(0:90:360, ["0°","90°","180°","270°","360°"]))
ax_13=Axis(figure[1,3][1,1],yticks=(0:90:360, ["0°","90°","180°","270°","360°"]))
ax_21=Axis(figure[2,1][1,1],yticks=(0:90:360, ["0°","90°","180°","270°","360°"]))
ax_22=Axis(figure[2,2][1,1],yticks=(0:90:360, ["0°","90°","180°","270°","360°"]))
ax_23=Axis(figure[2,3][1,1],yticks=(0:90:360, ["0°","90°","180°","270°","360°"]))
[xlims!(ax,0,68) for ax in [ax_11,ax_12,ax_13,ax_21,ax_22,ax_23]]
[ylims!(ax,0,360) for ax in [ax_11,ax_12,ax_13,ax_21,ax_22,ax_23]]
sideinfo1=Label(figure[1,0],"Linear Pressure Interpolation",rotation=pi/2,tellheight=false,fontsize=20)
sideinfo2=Label(figure[2,0],"Logarithmic Pressure Interpolation",rotation=pi/2,tellheight=false,fontsize=20)
topinfo1=Label(figure[0,1],"Carlotti",tellwidth=false,fontsize=30)
topinfo2=Label(figure[0,2],"Ciddor",tellwidth=false,fontsize=30)
topinfo3=Label(figure[0,3],"Mathar",tellwidth=false,fontsize=30)

h1=heatmap!(ax_11,h_knots.*majoraxis_earth,θ_knots,refractive_linear_carlotti)
Colorbar(figure[1,1][1,2],h1)
h2=heatmap!(ax_12,h_knots.*majoraxis_earth,θ_knots,refractive_linear_ciddor)
Colorbar(figure[1,2][1,2],h2)
h3=heatmap!(ax_13,h_knots.*majoraxis_earth,θ_knots,refractive_linear_mathar)
Colorbar(figure[1,3][1,2],h3)
h4=heatmap!(ax_21,h_knots.*majoraxis_earth,θ_knots,refractive_log_carlotti)
Colorbar(figure[2,1][1,2],h4)
h5=heatmap!(ax_22,h_knots.*majoraxis_earth,θ_knots,refractive_log_ciddor)
Colorbar(figure[2,2][1,2],h5)
h6=heatmap!(ax_23,h_knots.*majoraxis_earth,θ_knots,refractive_log_mathar)
Colorbar(figure[2,3][1,2],h6)

save("$(pwd())/output_refractive_index.png",figure)
figure=Figure(size=(1200,1200))
ax_11=Axis(figure[1,1][1,1],yticks=(0:90:360, ["0°","90°","180°","270°","360°"]))
ylims!(ax_11,0,360)
xlims!(ax_11,0,68)
ax_12=Axis(figure[1,2][1,1],yticks=(0:90:360, ["0°","90°","180°","270°","360°"]))
ylims!(ax_11,0,360)
xlims!(ax_11,0,68)


leftinfo1=Label(figure[0,1],"Linear Pressure Interpolation",tellwidth=false,fontsize=30)
rightinfo1=Label(figure[0,2],"Logarithmic Pressure Interpolation",tellwidth=false,fontsize=30)
h1=heatmap!(ax_11,h_knots.*majoraxis_earth.*km,θ_knots,pressure_linear.*10^-2.0.*hPa)
Colorbar(figure[1,1][1,2],h1)
h2=heatmap!(ax_12,h_knots.*majoraxis_earth.*km,θ_knots,pressure_log.*10^-2)
Colorbar(figure[1,2][1,2],h2)
Label(figure[2,:],L"\left|P_{linear}-P_{log}\right|_2",fontsize=30)
ax3=Axis(figure[3:4,1:2][1,1],yticks=(0:90:360, ["0°","90°","180°","270°","360°"]))
ylims!(ax3,0,360)
xlims!(ax3,0,68)
h3=heatmap!(ax3,h_knots.*majoraxis_earth.*km,θ_knots,10^-2.0.*sqrt.((pressure_log.-pressure_linear).^2))
Colorbar(figure[3:4,1:2][1,2],h3)
save("$(pwd())/output_pressure.png",figure)


# generate output vectors
T=Float64
# output rays to be used and avoid overwriting
#rays1=deepcopy(rays[:])
rays=rays[:]
rays1=similar(rays[:])
# position of the tangent point only updated if there is an intersection
tangent_quote=fill(T(Inf),length(rays1))
register=Vector{Register}(undef,length(rays1)) # index of the current position of the ray
# index of the current position of the ray
wedge_index=fill((1,1),length(rays1))
# find the first int

n₀=1.0


refractive=refractive_linear_carlotti

altitudes=Vector{Matrix{T}}(undef,6)


function getTangentquote(ray,local_index,local_register;hmax=hmax,levels=size(h_levels,1)-1)
  hmin=hmax ;              #every time we start a new ray, we initialize the minimum height to be the maximum height
  ############################
  inside_atmosphere=true;  # Initialize the search for the intersection
  tangent_found = false;   # Initialize the tangent point
  ############################
  previous_intersection=BottomIntersection();
  while inside_atmosphere
    ray,_,h,local_index,previous_intersection=code_intersection(ray,local_index;
      previous_intersection=previous_intersection);

    isouterspace=local_index[1]<1;
    isearth=local_index[1]>levels;

    inside_atmosphere= !(isouterspace || isearth);

    # update the position of the minimum
      if hmin>=h
        hmin=h;
        next!(local_register)
      elseif !tangent_found   # if the tangent point is not found
        setTangent!(local_register,ray.origin...)
        tangent_found=true;
      end
  end
  return local_register
end

function print_sequence(altitude_quote,orbit,majoraxis_earth=majoraxis_earth;error_limit=1)
  diff_altitude=altitude_quote-orbit.h.*majoraxis_earth;
  open("$(pwd())/output_sequence.dat","w") do io
    (_,ngeom)=size(altitude_quote)
    for seq in axes(altitude_quote,1)
      println(io,"SEQUENCE")
      println(io,"#")
      println(io,"           $(seq)")
      println(io,"  No of geometries")
      println(io,"#")
      println(io,"           $(ngeom)")
      println(io,"# z[km]  w[km]   angle[deg]  altitude_bianca[km]  altitude_marco[km]  diff[km]")
      for geom in axes(altitude_quote,2)
        info="$(orbit.z[seq,geom]*majoraxis_earth)  "
        info*="$(orbit.w[seq,geom]*majoraxis_earth)  "
        info*="$(orbit.ang[seq,geom])  "
        info*="$(orbit.h[seq,geom]*majoraxis_earth)  "
        info*="$(altitude_quote[seq,geom])  "
        info*="$(diff_altitude[seq,geom])"
        println(io,info)
      end
      println(io,"#")
    end
  end
  high_error_cartesian=findall(abs.(diff_altitude).>error_limit)
  return (diff_altitude,extrema(diff_altitude),high_error_cartesian)
end

 code_intersection(ray::Ray2D,index,refractive;kwargs...)=
new_intersection(ray, # ray
  index; # wedge index
  refractive_map=refractive, # map of refractive index
  h_levels=h_levels, # levels for i
  θ_radii=θ_radii, # radii angles for j
  scale_levels=scale_levels, # scaling mapping to unit circle
  line_radii=line_radii, # line radii
  tangent_quote=tangent_quote, # tangent quote
  register=register,kwargs...
)

for (iii,refr) in enumerate([refractive_linear_carlotti,refractive_linear_ciddor,refractive_linear_mathar,
            refractive_log_carlotti,refractive_log_ciddor,refractive_log_mathar])

  refractive=refr
  _find_first_intersection_ellipse(ray,n₀;kwargs...)=find_first_intersection_ellipse(ray,n₀;refractive_index_map=refractive,kwargs...)


  begin
    [
      begin
        #@info i
        rays1[i],flag,tangent_quote[i],wedge_index[i],_=
        _find_first_intersection_ellipse(ray,n₀)
        register[i]= Register(flag,rays1[i](0)...)
      end
      for (i,ray) in enumerate(rays)
    ];

    rays=reshape(rays,800,19)
    index=deepcopy(wedge_index)



    hmax=maximum(h_levels)  #initialize the maximum height



    NOPLOT=true
    if !NOPLOT
      fig,ax,ax2a,ax2b=initialize_raytracing_plot(h_levels,θ_radii)
      fig
      ax3=Axis(fig[4,1],
      xlabel="intersection point",
      ylabel="difference in altitude from line [km]"
      )
    end

    rays=reshape(rays,800,19)
    rays1=reshape(rays1,800,19)
    index=reshape(index,800,19)
    register=reshape(register,800,19)


    for i in eachindex(rays1)  #@info i
    r=rays1[i]

      ############################
      hmin=hmax               #every time we start a new ray, we initialize the minimum height to be the maximum height
      local_index=index[i]    #initialize the wedge index
      ############################

      ############################
      # TO PLOT THE RAYS
      ############################

      if !NOPLOT
        refraction_in_node=Float32[]
        intersections=Point2f[]
        directions=Point2f[]
        altitudes=Float32[]
        push!(intersections,Point2f(rays[i].origin...))

        push!(intersections,Point2f(r.origin...))
        push!(directions,Point2f(r.direction...))
        push!(refraction_in_node,refractive[local_index...])
        push!(altitudes,hmin)
      end

      ############################
      inside_atmosphere=true  # Initialize the search for the intersection
      tangent_found = false   # Initialize the tangent point
      ############################
      previous_intersection=BottomIntersection()

      #_count=0
      #while _count<20
      count=1
      (tangent_quote_w,tangent_quote_z)=get_origin(r)
      tangent_quote_number=0
      while inside_atmosphere
        count+=1
        r1,target,h,local_index,previous_intersection=code_intersection(r,local_index,refractive;
          previous_intersection=previous_intersection)
          # If the ray is in the atmosphere

        isouterspace=local_index[1]<1;
        isearth=local_index[1]>size(refractive,1);


        inside_atmosphere= !(isouterspace || isearth);

        #Refractive Index in Current Node
        n₀= if local_index[1]<1
          1.0;
        elseif local_index[1]>size(refractive,1)
          0.01;
        else
          refractive[local_index...];
        end;
        hh,tt=convert(LLA2D{NormalizedEarth},ECEF2D{NormalizedEarth}(r1(0)...)) |> x-> (x.h,x.θ)
        #@info "hextr: $(hh.*majoraxis_earth),hout: $(h.*majoraxis_earth), hmin: $(hmin.*majoraxis_earth),previous_intersection: $(previous_intersection)"
        #@info "hextr: $(hh.*majoraxis_earth),hout: $(h.*majoraxis_earth), θ_ext: $(tt), hmin: $(hmin.*majoraxis_earth)"
          if hh<hmin
            hmin=hh
            tangent_quote_number=count
            (tangent_quote_w,tangent_quote_z)=get_origin(r1)
            flag=true
          end

        if !NOPLOT
          push!(refraction_in_node,n₀);
          push!(intersections,Point2f(r1.origin...));
          push!(directions,Point2f(r1.direction...));
          push!(altitudes,h);
        end
          r=r1
      end


      settangentquote!(register[i],tangent_quote_w,tangent_quote_z,tangent_quote_number)


      #==================================#
      if !NOPLOT
        scatter!(ax2b,Point2(tangent_quote_number,hmin*majoraxis_earth),marker='x',markersize=20,color=:black)

        scatterlines!(ax,intersections[1:end])
        _A=intersections[2]
        _B=intersections[end]
        _AB=_B-_A
        _f(x)=_AB[2]/_AB[1]*(x-_A[1])+_A[2]
        scatterlines!(ax3,[let
          abs(_f(intersection[1])-intersection[2]).*majoraxis_earth
        end
        for intersection in intersections[2:end]
        ])


        scatterlines!(ax2b,altitudes.*majoraxis_earth)
        scatterlines!(ax2a,altitudes.*majoraxis_earth,refraction_in_node)
      end
    end
    if !NOPLOT
      save("$(pwd())/output_highest_error.png",fig)
    end


    #(diff_altitude,extdiff,high_error_cartesian)=print_sequence(StructArray(register).h.*majoraxis_earth,orbit;error_limit=0.2)
  end

  altitudes[iii]=StructArray(register).h.*majoraxis_earth
end



(diff_altitude,extdiff,high_error_cartesian)=print_sequence(altitudes[1],orbit;error_limit=0.2)

size(high_error_cartesian)

extrema(diff_altitude./(orbit.h.*majoraxis_earth).*100)

df=DataFrame(
  seq=repeat(1:800,outer=19),
  geom=repeat(1:19,inner=800),
  bianca=(orbit.h[:]).*majoraxis_earth,
  carlotti_linear=altitudes[1][:],
  absdiff_bianca_marco_abs=abs.(diff_altitude[:]),
  absdiff_bianca_marco_perc=abs.(diff_altitude[:])./(orbit.h[:].*majoraxis_earth).*100,
  ciddor_linear=altitudes[2][:],
  mathar_linear=altitudes[3][:],
  carlotti_log=altitudes[4][:],
  ciddor_log=altitudes[5][:],
  mathar_log=altitudes[6][:]
)

df[!,:linear_vs_log_ciddor_perc]=df[!,:ciddor_linear].-df[!,:ciddor_log] |> x-> @. x/df[!,:ciddor_linear]*100 |> x-> @. abs(x) |> x-> @. ustrip(x)
df[!,:linear_vs_log_mathar_perc]=df[!,:mathar_linear].-df[!,:mathar_log] |> x-> @. x/df[!,:mathar_linear]*100 |> x-> @. abs(x) |> x-> @. ustrip(x)
df[!,:linear_vs_log_carlotti_perc]=df[!,:carlotti_linear].-df[!,:carlotti_log] |> x-> @. x/df[!,:carlotti_linear]*100 |> x-> @. abs(x) |> x-> @. ustrip(x)
df[!,:carlotti_vs_ciddor_linear_perc]=df[!,:carlotti_linear].-df[!,:ciddor_linear] |> x-> @. x/df[!,:carlotti_linear]*100 |> x-> @. abs(x)  |> x-> @. ustrip(x)
df[!,:carlotti_vs_mathar_linear_perc]=df[!,:carlotti_linear].-df[!,:mathar_linear] |> x-> @. x/df[!,:carlotti_linear]*100 |> x-> @. abs(x) |> x-> @. ustrip(x)
df[!,:carlotti_vs_ciddor_log_perc]=df[!,:carlotti_log].-df[!,:ciddor_log] |> x-> @. x/df[!,:carlotti_log]*100 |> x-> @. abs(x) |> x-> @. ustrip(x)
df[!,:carlotti_vs_mathar_log_perc]=df[!,:carlotti_log].-df[!,:mathar_log] |> x-> @. x/df[!,:carlotti_log]*100 |> x-> @. abs(x) |> x-> @. ustrip(x)

summdf=describe(df[!,[:absdiff_bianca_marco_abs,:absdiff_bianca_marco_perc,
:linear_vs_log_ciddor_perc,:linear_vs_log_mathar_perc,:linear_vs_log_carlotti_perc,:carlotti_vs_ciddor_linear_perc,:carlotti_vs_mathar_linear_perc,:carlotti_vs_ciddor_log_perc,:carlotti_vs_mathar_log_perc]],
:mean, :min,:max,:median,:std)

using XLSX

dfc=df[!,[:seq,:geom,:bianca,:carlotti_linear]]
dfc[!,:diff_bianca_marco]=df[!,:bianca].-df[!,:carlotti_linear]
dfc[!,:diff_bianca_marco_perc]=dfc[!,:diff_bianca_marco]./dfc[!,:bianca].*100
extrema(dfc[!,:diff_bianca_marco_perc])
linindex=[(seq_geom[1]-1).*19+seq_geom[2] for seq_geom in sort(high_error_cartesian)]

XLSX.writetable("$(pwd())/250113_output_ray_tracing.xlsx","raw_data"=>df,"summary"=>summdf,"error above 200m"=>dfc1;overwrite=true)

dfc1=subset(dfc,:diff_bianca_marco_perc => x->@. abs(x)>0.)
sort(dfc,:diff_bianca_marco_perc)

summdf[!,:abs_diff]=["bianca vs marco km","bianca vs marco %",
"linear vs log ciddor %","linear vs log mathar %","linear vs log carlotti %",
"carlotti vs ciddor linear %","carlotti vs mathar linear %","carlotti vs ciddor log %","carlotti vs mathar log %"]

summdf=summdf[!,[:abs_diff,:mean,:std,:min,:max,:median]]

sort(dfc,:diff_bianca_marco_perc)

figure=Figure()

axis=Axis3(figure[1,1][1,1],xlabel="Sequence",
ylabel="Limb Angle [deg]",
zlabel="Altitude [km]")

d1=T[]
ang1=T[]
[
  begin
    @info orb
    push!(d1,sqrt(((orb.h).*majoraxis_earth-alt)^2))
    push!(ang1,orb.ang)
  end
  for (orb,alt) in zip(orbit,altitudes[1])
]
d1=reshape(d1,800,19)
ang1=reshape(ang1,800,19)


pp=[
  Point(i,ang1[i,j],d1[i,j])  for i in axes(d1,1), j in axes(d1,2)
]

cc=[p[3] for p in pp]
scatter!(axis,pp[:],color=cc[:])
save("$(pwd())/output_error_carlotti_vs_seq_limb_top.png",figure)


figure=Figure(size=(800,1200))
ax1=Axis(figure[1,1:2][1,1],xlabel="error",ylabel="frequency",title="Error Frequency")
ax2=Axis(figure[2,1][1,1],xlabel="angle",ylabel="frequency",title="Latitude with Error>0.19km",
xticks=([   0,
45,
90,
135,
180,
225,
270,
315,
360],["0°","45°","90°","135°","180°","225°","270°","315°","360°"]))
xlims!(ax2,0,360)
ax3=Axis(figure[2,2][1,1],xlabel="angle",ylabel="frequency",title="Latitude with Error<0.01km",
xticks=([   0,
45,
90,
135,
180,
225,
270,
315,
360],["0°","45°","90°","135°","180°","225°","270°","315°","360°"]))
xlims!(ax3,0,360)



hb2=hist!(ax1,((orbit.h[:]).*majoraxis_earth-altitudes[1][:]),bins=50,color=:values)


orb_bigdiff=orbit[findall(abs.(diff_altitude).>0.19)]
orb_smalldiff=orbit[findall(abs.(diff_altitude).<0.01)]


hist!(ax2,[let
  convert(LLA2D{NormalizedEarth},ECEF2D{NormalizedEarth}(o.w,o.z)) |>
  o-> o.θ #-90 |> o-> mod(o,360)
end
for o in orb_bigdiff],bins=100,color=:values)

hist!(ax3,[let
  convert(LLA2D{NormalizedEarth},ECEF2D{NormalizedEarth}(o.w,o.z)) |>
  o-> o.θ #-90 |> o-> mod(o,360)
end
for o in orb_smalldiff],bins=100,color=:values)

save("$(pwd())/output_hist_error_bianca_vs_mine.png",figure)


subset(df,:seq=> s-> s.==800)


LinRange(0,360,1000) |> x-> lines(x, hypot.(cosd.(x),sind.(x)./0.988))


LinRange(-1,1,1000) |> x-> lines!(acosd.(x),hypot.(x,_y.(x)))

sqrt(cosd(90)^2+(sind(90)*0.988)^2)


_y(x)= cos(asin(x))./0.988
b=0.4
h=0.5
x=LinRange(0,360,1000)

  x1=@. cosd(x)
  y1= @. b*sind(x)

  N=@. 1/sqrt(x1^2+y1^2)
  x2=@. N*cosd(x)
  y2=@. b^2*N*sind(x)



  hx=@. h*cosd(x)
  hy=@. h*sind(x)

  n=@. hypot(cosd(x),sind(x)/b)

  hx2=@. h/n*(x2-h)
  hy2=@. h/n*sind(x)/b


  figure=Figure()
  ax=Axis(figure[1,1][1,1],xlabel="x",ylabel="y")
  lines!(ax,x1+hx2,y1+hy2)
  lines!(ax,x1,y1)
  lines!(ax,x2+hx,y2+hy)


extrema(@. x1^2+y1^2)
