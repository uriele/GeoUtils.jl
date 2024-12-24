using GeoUtils
const MODEL=Ref{AirModel}(Carlotti())
const INTERPOLATION=Ref{AbstractPressureInterpolation}(LinearPressure())

setModel(A::AirModel)= MODEL[]=A
getModel()=MODEL[]
setModel()=MODEL[]=(Carlotti())
setPressureInterpolation(A::AbstractPressureInterpolation)= INTERPOLATION[]=A
setPressureInterpolation()= INTERPOLATION[]=(LinearPressure())
getPressureInterpolation()=INTERPOLATION[]


include("initialization_script.jl")


# generate output vectors
T=Float64
# output rays to be used and avoid overwriting
rays1=deepcopy(rays[:])
# position of the tangent point only updated if there is an intersection
tangent_quote=fill(T(Inf),length(rays1))
register=Vector{Register}(undef,length(rays1)) # index of the current position of the ray
# index of the current position of the ray
wedge_index=fill((1,1),length(rays1))
# find the first int
[
  begin
    rays1[i],flag,tangent_quote[i],wedge_index[i]=find_first_intersection_ellipse(ray,1.0)
    register[i]= Register(flag,rays1[i](0)...)
  end
  for (i,ray) in enumerate(rays[:])
];
# find all the other Intersection
index=deepcopy(wedge_index)

code_intersection(ray::Ray2D,index;kwargs...)=
new_intersection(ray, # ray
  index; # wedge index
  refractive_index_map=refractive, # map of refractive index
  h_levels=h_levels, # levels for i
  θ_radii=θ_radii, # radii angles for j
  scale_levels=scale_levels, # scaling mapping to unit circle
  line_radii=line_radii, # line radii
  tangent_quote=tangent_quote, # tangent quote
  register=register,kwargs...
)



fig,ax,ax2a,ax2b=initialize_raytracing_plot(h_levels,θ_radii)
fig


for i in  101:500
#i=10
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
    inside_atmosphere= !(isouterspace || isearth)

    # Refractive Index in Current Node
    n₀= if index[1]<1
      1.0
orbi    elseif index[1]>size(refractive,1)
      NaN
    else
      refractive[index...]
    end

      push!(refraction_in_node,n₀)
      push!(intersections,Point2f(r.origin...))
      push!(directions,Point2f(r.direction...))
      push!(altitudes,h)

      # update the position of the minimum
      if hmin>=h
         hmin=h;
         next!(register[i])
      elseif !tangent_found   # if the tangent point is not found
         setTangent!(register[i],r(0)...)
         tangent_found=true
         register[i]
      end
      @debug "h:$(h*a)km  $wm $(origin(r)) $(direction(r))"
      @debug register[i]
      @debug hypot(direction(r)...)
      @debug hypot(r(0)...)
  end
  scatterlines!(ax,intersections)
  scatterlines!(ax2b,altitudes.*a;label="Ray $i")
  scatterlines!(ax2a,altitudes.*a,refraction_in_node;label="Ray $i")
end

save("$(pwd())/ray_tracing_test.png",fig)
