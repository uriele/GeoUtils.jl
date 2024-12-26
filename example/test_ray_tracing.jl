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
setModel(NoAtmosphere())
#setModel(Carlotti())
getDebugIntersection()
# only radii
setDebugIntersection(1)
include("initialization_script.jl")


# generate output vectors
T=Float64
# output rays to be used and avoid overwriting
#rays1=deepcopy(rays[:])
rays1=similar(rays[1,:])
# position of the tangent point only updated if there is an intersection
tangent_quote=fill(T(Inf),length(rays1))
register=Vector{Register}(undef,length(rays1)) # index of the current position of the ray
# index of the current position of the ray
wedge_index=fill((1,1),length(rays1))
# find the first int

#n₀=refractive[1,1]
n₀=1.0

[
  begin
    rays1[i],flag,tangent_quote[i],wedge_index[i],_=
    find_first_intersection_ellipse(ray,n₀)
    register[i]= Register(flag,rays1[i](0)...)
  end
  for (i,ray) in enumerate(rays[1,:])
];



# find all the other Intersection
index=deepcopy(wedge_index)

code_intersection(ray::Ray2D,index;kwargs...)=
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

#==========#
_angleme(r::Ray2D) = r(0)-r(1) |> x->(x[2],x[1])  |> x-> atand(x...)
_angleme(r1::Ray2D,r::Ray2D) = r1(0)-r(0) |> x->(x[2],x[1])  |> x-> atand(x...)
##############

fig,ax,ax2a,ax2b=initialize_raytracing_plot(h_levels,θ_radii)
fig
ax3=Axis(fig[4,1])
for i in  1:19
  r=rays1[i]

  ############################
  hmin=maximum(h_levels)  #initialize the maximum height
  index=wedge_index[i]    #initialize the wedge index
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
  push!(refraction_in_node,refractive[index...])
  push!(altitudes,hmin)

  ############################
  inside_atmosphere=true  # Initialize the search for the intersection
  tangent_found = false   # Initialize the tangent point
  ############################
  previous_intersection=BottomIntersection()
while inside_atmosphere
  r,target,h,index,previous_intersection=code_intersection(r,index;
      previous_intersection=previous_intersection)
      # If the ray is in the atmosphere


      isouterspace=index[1]<1;
    isearth=index[1]>size(refractive,1);

    @debug "isouterspace: $isouterspace, isearth: $isearth"
    inside_atmosphere= !(isouterspace || isearth);

    # Refractive Index in Current Node
    n₀= if index[1]<1
      1.0;
    elseif index[1]>size(refractive,1)
      3.0;
    else
      refractive[index...];
    end;

      push!(refraction_in_node,n₀);
      push!(intersections,Point2f(r.origin...));
      push!(directions,Point2f(r.direction...));
      push!(altitudes,h);

      # update the position of the minimum
      if hmin>=h
         hmin=h;
         next!(register[i])
      elseif !tangent_found   # if the tangent point is not found
         setTangent!(register[i],r(0)...)
         tangent_found=true
         register[i]
      end;

      @info "h:$(h*majoraxis_earth)km  $(_angleme(r))"
      @info "previous_intersection"
      @info register[i]

  end
  scatterlines!(ax,intersections[1:end])

  scatterlines!(ax3,[intersections[1],intersections[end]])
  scatterlines!(ax3,intersections[1:end])
  scatterlines!(ax2b,altitudes.*majoraxis_earth;label="Ray $i")
  scatterlines!(ax2a,altitudes.*majoraxis_earth,refraction_in_node;label="Ray $i")

end
xlims!(ax,(0.,0.5))
ylims!(ax,(0.5,1))

iijj=wedge_index[i]

rl=line_radii[iijj[2]]
rr=line_radii[iijj[2]+1]
_fl(x)=-(rl.coeff_a*x+rl.coeff_c)/rl.coeff_b
_fr(x)=-(rr.coeff_a*x+rr.coeff_c)/rr.coeff_b

fig=Figure()
ax=Axis(fig[1,1])
index=(index[1],index[2]-1)
lines!(ax,[line_radii[index[2]].PointA,line_radii[index[2]].PointB],linestyle=:dash)
lines!(ax,[line_radii[index[2]+1].PointA,line_radii[index[2]+1].PointB],linestyle=:dash)
lines!(ax,[r1(0),r(1)])
scatter!(ax,r1(0))

convert(LLA2D{NormalizeEarth},ECEF2D{NormalizeEarth}(...))
