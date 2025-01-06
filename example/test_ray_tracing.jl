using GeoUtils
const MODEL=Ref{AirModel}(Carlotti())
const INTERPOLATION=Ref{AbstractPressureInterpolation}(LinearPressure())

setModel(A::AirModel)= MODEL[]=A
getModel()=MODEL[]
setModel()=MODEL[]=(Carlotti())
setPressureInterpolation(A::AbstractPressureInterpolation)= INTERPOLATION[]=A
setPressureInterpolation()= INTERPOLATION[]=(LinearPressure())
getPressureInterpolation()=INTERPOLATION[]
getPressureInterpolation()
getModel()
#setModel(NoAtmosphere())
setModel(Carlotti())
getDebugIntersection()
# only radii
#setDebugIntersection(0)
include("initialization_script.jl")

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

[
  begin
    #@info i
    rays1[i],flag,tangent_quote[i],wedge_index[i],_=
    find_first_intersection_ellipse(ray,n₀)
    register[i]= Register(flag,rays1[i](0)...)
  end
  for (i,ray) in enumerate(rays)
];

index=deepcopy(wedge_index)

@inline code_intersection(ray::Ray2D,index;kwargs...)=
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

hmax=maximum(h_levels)  #initialize the maximum height

@inline function getTangentquote(ray,local_index,local_register;hmax=hmax,levels=size(h_levels,1)-1)
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


using Profile,PProf
using BenchmarkTools
  #
  @profile begin
  for (ray,local_register,local_index) in zip(rays1,register,index)
    ############################
    local_register=getTangentquote(ray,local_index,local_register;hmax=hmax)
  end
end
@profile  begin
  Threads.@threads for i in eachindex(rays1) #,register,index)
    register[i]=getTangentquote(rays1[i],index[i],register[i];hmax=hmax)
  end
end
pprof(;webport=57777)
nt=Threads.nthreads()



figure=Figure()
ax=Axis(figure[1,1])
scatterlines!(ax,altitude_quote[1,:])
scatterlines!(ax,orbit.h[1,:].*majoraxis_earth)


high_error=findall(abs.(diff_altitude[:]).>1)
high_error_cartesian=findall(abs.(diff_altitude).>1)

diff_altitude1=copy(diff_altitude)
diff_altitude1[abs.(diff_altitude1).<1e-4].=0.0



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

altitude_quote= round.(altitude_quote,digits=2)
diff_altitude,max_error,position=print_sequence(altitude_quote,orbit,majoraxis_earth;error_limit=0.5)
psort=sort([(c[1],c[2]) for c in position])

max_error
 mat1=zeros(5,size(psort,1))
[mat1[:,i]=[c[1] c[2] diff_altitude[c[1],c[2]] orbit.h[c[1],c[2]]*majoraxis_earth altitude_quote[c[1],c[2]] ] for (i,c) in enumerate(psort)]
mat1
mat1'
figure=Figure()
ax=Axis(figure[1,1],
xlabel="Difference in Altitude km",
ylabel="Frequency")
hist!(ax,abs.(diff_altitude[:]),normalization=:probability,bins=100,color=:values,bar_labels = :values,
label_formatter=x-> round(x, digits=4))
xlims!(ax,0,0.5)
ylims!(ax,0,0.8)
save("$(pwd())/output_histogram.png",figure)


figure=Figure()
ax=Axis(figure[1,1],
xlabel="Difference in Altitude km",
ylabel="Number")
his=hist!(ax,abs.(diff_altitude[abs.(diff_altitude).>0.5]),color=:values)
Colorbar(figure[1,2],his)
save("$(pwd())/output_histogram_over0.5.png",figure)



fig,ax,ax2a,ax2b=initialize_raytracing_plot(h_levels,θ_radii)
fig
ax3=Axis(fig[4,1],
xlabel="intersection point",
ylabel="difference in altitude from line [km]"
)


for (ii,i) in enumerate(high_error)
  #@info i
  r=rays1[i]

  ############################
  hmin=hmax               #every time we start a new ray, we initialize the minimum height to be the maximum height
  local_index=index[i]    #initialize the wedge index
  ############################

  ############################
  # TO PLOT THE RAYS
  ############################
  refraction_in_node=Float32[]
  intersections=Point2f[]
  directions=Point2f[]
  altitudes=Float32[]
  push!(intersections,Point2f(rays[i].origin...))

  push!(intersections,Point2f(r.origin...))
  push!(directions,Point2f(r.direction...))
  push!(refraction_in_node,refractive[local_index...])
  push!(altitudes,hmin)

  ############################
  inside_atmosphere=true  # Initialize the search for the intersection
  tangent_found = false   # Initialize the tangent point
  ############################
  previous_intersection=BottomIntersection()

  #_count=0
  #while _count<20
  while inside_atmosphere

    r1,target,h,local_index,previous_intersection=code_intersection(r,local_index;
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

      push!(refraction_in_node,n₀);
      push!(intersections,Point2f(r.origin...));
      push!(directions,Point2f(r.direction...));
      push!(altitudes,h);

      r=r1
  end
  #==================================#
  scatterlines!(ax,intersections[1:end])
  _A=intersections[2]
  _B=intersections[end]
  _AB=_B-_A
  _f(x)=_AB[2]/_AB[1]*(x-_A[1])+_A[2]
  scatterlines!(ax3,[let
    abs(_f(intersection[1])-intersection[2]).*majoraxis_earth
  end
  for intersection in intersections[2:end]
  ],label="Ray ($(high_error_cartesian[ii][1]),$(high_error_cartesian[ii][2]))")


  scatterlines!(ax2b,altitudes.*majoraxis_earth;
  label="Ray ($(high_error_cartesian[ii][1]),$(high_error_cartesian[ii][2]))")
  scatterlines!(ax2a,altitudes.*majoraxis_earth,refraction_in_node;
  label="Ray ($(high_error_cartesian[ii][1]),$(high_error_cartesian[ii][2]))")
end

save("$(pwd())/output_highest_error.png",fig)
