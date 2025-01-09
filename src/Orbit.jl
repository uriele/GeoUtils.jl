
#TO DO if necessary refactor and add a 3D ellipsoid and a Radius plane
# For now, it is not needed

"""
  Ellipsoid(majoraxis,minoraxis)

Generate a structure with scale and unscale functions based on the major and minor axis.
If majoraxis<minoraxis, they are inverted.

Ellipsoid is used to compute the intersection of the Ray2D with the atmosphere.
"""
struct Ellipsoid{T<:IEEEFloat}
  # The map is assumed to be centered at the origin so I do not need to store the center
  # and the affine translation
  scale ::SVec2{T}
  unscale::SVec2{T}

  @inline function Ellipsoid(majoraxis_earth::T,minoraxis_earth::T) where T
    (majoraxis_earth,minoraxis_earth)= (majoraxis_earth>=minoraxis_earth) ? (majoraxis_earth,minoraxis_earth) : (minoraxis_earth,majoraxis_earth)
    unscale=SVec2{T}(majoraxis_earth,minoraxis_earth)
    scale=one(T)./unscale
    return new{T}(scale,unscale)
  end
end
Ellipsoid(majoraxis_earth::T) where T = Ellipsoid(majoraxis_earth,majoraxis_earth)


"""
    Radius([datum=NormalizedEarth],θ)

Generate a radius for the earth Ellipse that is used to compute the intersection of the Ray2D with the atmosphere using the angle θ and the datum.
If the datum is not provided, it will use the NormalizedEarth datum.

"""
struct Radius{T<:IEEEFloat}
  PointA::SVec2{T}
  PointB::SVec2{T}
end


@inline function _create_radii_from_θ(datum::Datum, θ::T) where {Datum,T<:IEEEFloat}
  majoraxis_earth= majoraxis(ellipsoid(datum)) |> ustrip
  minoraxis_earth= minoraxis(ellipsoid(datum)) |> ustrip
  minoraxis_earth/=majoraxis_earth
  origin=convert(ECEF2D{datum},LLA2D{datum}(0.,θ)) |> x-> SVec2([x.w,x.z])
  direction=MVec2(origin...)
  #@info "origin: $(origin) direction: $(direction)"
  _normal_vector!(direction,zero(T),minoraxis_earth)
  #@info "n: $(direction)"
  return (origin...,(origin+direction)...)
end

@inline _create_radii_from_θ(datum::Datum, θ) where {Datum}= _create_radii_from_θ(datum, float(θ))
@inline _create_radii_from_θ(θ)=_create_radii_from_θ(NormalizedEarth,θ)

Radius(a::T,b::T,c::T,d::T) where T = Radius(SVec2(a,b),SVec2(c,d))
Radius(datum::Datum,θ::T) where {Datum,T} = Radius(_create_radii_from_θ(datum,θ)...)
Radius(θ::T) where T = _create_radii_from_θ(θ)

"""
    distance_from_unit_circle(ray::Ray2D{T})::T where T

Calculate the distance from the unit circle to a ray or a point
"""
distance_from_unit_circle(ray::R2) where {R2<:AbstractRay2D{T}} where  T = _distance_from_unit_circle(ray.origin,ray.direction)::T

@inline function _distance_from_unit_circle(origin::SVec2{T},direction::SVec2{T})::T where T

  _inv_hyp_direct= one(T)/hypot(direction...)
  direction=direction*_inv_hyp_direct
  C=dot(origin,origin)-one(T)
  halfB=-dot(origin,direction)


  μ =_fastquadratic(halfB,C)

  t1=μ+halfB |> x-> x>0 ? x.*_inv_hyp_direct : Inf
  t2=-μ+halfB |> x-> x>0 ? x.*_inv_hyp_direct : Inf
  return min(t1,t2)

end


"""
  distance_from_radii(ray::Ray2D{T},r::Radius{T})::T where T

Calculate the distance from a ray to one of Earth Radii or from a line with coefficients coeff_a,coeff_b,coeff_c

```math
coeff_a*x+coeff_b*y+coeff_c=0
```
"""
distance_from_radii(ray::Ray2D{T},radius::Radius{T}) where T=
  distance_from_segment(ray::Ray2D{T},radius.PointA,radius.PointB)


@inline function _distance_from_segment(
  A::SV1,B::SV2,
  C::SV3,D::SV4
  )::T where {SV1<:StaticVector{2,T},
  SV2<:StaticVector{2,T},
  SV3<:StaticVector{2,T},
  SV4<:StaticVector{2,T}} where T
  # find intersection between segments AB,CD
  # Assume unit ray
  # A+t*(B-A) = p = C+s*(D-C)
  # [B-A  C-D]*[t;s]=[C-A]
  #
  M=[(B-A) (C-D)];
  y=(C-A);
  #
  _,R=qr([M y]);
  τ = atol(T);
  # Check L1 norm of the residual
  rₐ =sum(>(τ),sum(abs,R, dims=2))
  # Calculate rank of A
  r = sum(>(τ),sum(abs,view(R,:,1:2), dims=2))
  r= (r==rₐ)*r+(r!=rₐ)*-1
  # note: if it's collinear we will need to match the
  # fitting with the ellipse, so for our purposes collinear
  # is not intersecting
  ((r>=2) && (M\y))[1] |> x-> (x>0) ? x : T(Inf)
end

"""
    distance_from_segment(ray::Ray2D{T},C::SVec2{T},D::SVec2{T})::T where T

Calculate the distance from a segment AB to a segment CD
"""
distance_from_segment(ray::Ray2D{T},C::SVec2{T},D::SVec2{T}) where T = _distance_from_segment(ray.origin,ray(1),C,D)



"""
    advance(intersection::I,ray::Ray2D{T},obj)::T where {I<:IntersectionLocation,T}

Advance the ray to the intersection point with the intersection object I.
"""
advance(::I,ray::Ray2D{T},obj::X) where {I<:IntersectionLocation,T,X} = throw(ArgumentError("advance function not implemented for $(I) and $(X)"))

function advance(::L, r::Ray2D{T},scale::SVec2{T})::T where {L<:LevelIntersection,T}
  _orig=scale.*r.origin
  _direc=scale.*r.direction
  return _distance_from_unit_circle(_orig,_direc)
end
function advance(::R, ray::Ray2D{T},radius::Radius)::T where {R<:RadiusIntersection,T}
  return _distance_from_segment(ray.origin,ray(1),radius.PointA,radius.PointB)
end

# Constant of the shift inside the medium
const SHIFTORIGIN=1e-12

# direction of the outward and inward normal, used to correct the normal to the point of intersection
const OUTWARD_NORMAL=1.0
const INWARD_NORMAL=-1.0


"""
  bend(::IntersectionLocation,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) where T<:IEEEFloat

Bend the ray at the intersection point I with the object x. The intersection point is at distance t from the origin of the ray.
The bend point is shifted by a small amount (1e-12) inside the interface to avoid numerical issues.

# Input
- `::IntersectionLocation`: The intersection object, can be a
  1. TopIntersection()
  2. BottomIntersection()
  3. LeftIntersection()
  4. RightIntersection()
  5. TopLeftIntersection()
  6. TopRightIntersection()
  7. BottomLeftIntersection()
  8. BottomRightIntersection()
- `ray::Ray2D`: The ray to be bent
- `t::T`: The distance from the origin of the ray to the intersection point
- `n₀::T=1.0`: The refractive index of the medium where the ray is coming from
- `n₁::T=1.0`: The refractive index of the first medium where the ray is going to
- `n₂::T=1.0`: The refractive index of the second medium where the ray is going to [only for corners, otherwise it is ignored]
- `h::T=0.0`: The height of the atmosphere [only for intersections with the ellipse, otherwise it is ignored]

# Output
- `Ray2D{T}`: The new ray after the bend
- `Bool`: If the ray is reflected or not
"""
bend(::IsIntersection,x::T,ray::Ray2D,t;kwargs...) where T= error("bend function not implemented for $(T)")
bend(::IsIntersection,x::Type{T},ray::Ray2D,t;kwargs...) where T= error("bend function not implemented for $(T)")
bend(x,ray::Ray2D,t;kwargs...) = bend(IntersectionStyle(x),x,ray,t;kwargs...)
#TopLevelIntersection -1
#--------------------
#         |  n
#         V
####################
bend(::TopIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)    = _bend_ellipse(ray,t,n₀,n₁,h,INWARD_NORMAL)
#RightRadiusIntersection -1
#
#      n    |
# <---------|
#           |
####################
bend(::RightIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)  = _bend_radii(ray,t,n₀,n₁,0.0,INWARD_NORMAL)


#BottomLevelIntersection 1
#         ^
#         |  n
#--------------------
####################

bend(::BottomIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) = _bend_ellipse(ray,t,n₀,n₁,h,OUTWARD_NORMAL)
#LeftRadiusIntersection 1
# |
# |----------> n
# |
#####################
bend(::LeftIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)   = _bend_radii(ray,t,n₀,n₁,0.0,OUTWARD_NORMAL)

#LeftTop 1 -1
bend(::LeftTopIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)=
  bend(LeftIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(TopIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)
#LeftBottom 1 1
bend(::LeftBottomIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(LeftIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(BottomIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)

#RightTop -1 -1
bend(::RightTopIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(RightIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(TopIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)
#RightBottom -1 1
bend(::RightBottomIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)=
  bend(RightIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(BottomIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)

#TopLeft -1 1
bend(::TopLeftIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(TopIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(LeftIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)
#TopRight -1 1
bend(::TopRightIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(TopIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(RightIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)

#BottomLeft 1 1
bend(::BottomLeftIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(BottomIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(LeftIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)

#BottomRight 1 -1
bend(::BottomRightIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(BottomIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(RightIntersection(),x,0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)


@inline function _bend_initialize(n₀::T,n₁::T) where T
  n₀₁=n₀/n₁
  return (n₀₁,n₀₁*n₀/n₁)
end


# bend the ellipse and return the new ray AND if it is reflected
@inline function _bend_ellipse!(N::MVec2{T},neworigin::SV,direction::SV, n₀::T,n₁::T, h::T,
  normalized_minoraxis_earth::T, outward::T) where {SV<:StaticVector{2,T}} where T
  ## shift the origin to avoidnumerical issues
  (n₀₁,n₀₁²)=_bend_initialize(n₀,n₁)
  N.=neworigin
  _normal_vector!(N,h,normalized_minoraxis_earth).*outward |> x-> x/hypot(x...)
  ray_out,isReflected=_bend_common(neworigin,direction,N,n₀₁,n₀₁²)
  return (ray_out,isReflected)
end

# bend the ellipse by changing the ray INPLACE AND return if it is reflected
@inline function _bend_ellipse!!(N::MVec2{T},ray::MRay2D{T}, n₀::T,n₁::T,
   h::T,normalized_minoraxis_earth::T, outward::T) where T
  ## shift the origin to avoidnumerical issues
  (n₀₁,n₀₁²)=_bend_initialize(n₀,n₁)
  N.=ray.origin
  _normal_vector!(N,h,normalized_minoraxis_earth).*outward |> x-> x/hypot(x...)
  return _bend_common!(ray,N,n₀₁,n₀₁²)
end


# bend the radii and return the new ray AND if it is reflected
@inline function _bend_radii!(N::MVec2{T},neworigin::SV,direction::SV,
  n₀::T,n₁::T,normalized_minoraxis_earth::T, outward::T) where {SV<:StaticVector{2,T}} where T
  (n₀₁,n₀₁²)=_bend_initialize(n₀,n₁)
  N.=neworigin
  _tangent_vector!(N,zero(T),normalized_minoraxis_earth::T).*outward |> x-> x/hypot(x...)
  ray_out,isReflected=_bend_common(neworigin,direction,N,n₀₁,n₀₁²)
  return (ray_out,isReflected)
end


# bend the radii by changing the ray INPLACE AND return if it is reflected
@inline function _bend_radii!!(N::MVec2{T},ray::MRay2D{T},n₀::T,n₁::T,normalized_minoraxis_earth::T, outward::T) where T

  @info "before: $(ray.origin) $(N)"
  (n₀₁,n₀₁²)=_bend_initialize(n₀,n₁)
  N.=ray.origin
  @info "after: $(ray.origin) $(N)"
  _tangent_vector!(N,zero(T),normalized_minoraxis_earth::T).*outward |> x-> x/hypot(x...)
  return _bend_common!(ray,N,n₀₁,n₀₁²)
end



# return both the new ray and if it is reflected
@inline function _bend_common(origin::SVec2{T},_direction::SVec2{T},
  N::SVec2{T},n₀₁::T,n₀₁²::T) where T
  # Base condition
  isReflected=false

  cosθₙ=-N⋅_direction

  sinθₙ² =n₀₁²*(1-cosθₙ^2)

  if sinθₙ²≤1
    newdir=n₀₁*_direction+(n₀₁*cosθₙ-sqrt(1-sinθₙ²))*N
  else
    newdir=_direction-2*cosθₙ*N
    isReflected=true
  end
  # SHIFT ORIGIN TO AVOID NUMERICAL ISSUES
  return Ray2D(origin+newdir*SHIFTORIGIN,newdir),isReflected
end

# modifies the ray in place and return if it is reflected
@inline function _bend_common!(ray::MRay2D{T},
  N::SV2,n₀₁::T,n₀₁²::T) where {SV2<:StaticVector{2,T}} where T
  # Base condition
  isReflected=false

  cosθₙ=-N⋅ray.direction

  sinθₙ² =n₀₁²*(1-cosθₙ^2)

  if sinθₙ²≤1
    ray.direction.=n₀₁*ray.direction+(n₀₁*cosθₙ-sqrt(1-sinθₙ²))*N
  else
    ray.direction.=ray.direction-2*cosθₙ*N
    isReflected=true
  end
  ray.origin.=ray.origin+ray.direction*SHIFTORIGIN
  # SHIFT ORIGIN TO AVOID NUMERICAL ISSUES
  return isReflected
end


###################################################################################################
# VISUAL REPRESENTATION OF TANGENT AND NORMAL VECTORS
###################################################################################################
#          n=(-nx,ny), t=(ny,nx)
#           ^  n                                                                         ^ t(0,1)
#           | (0,1)                                n \    / t                    n(-1,0) |
#           |                              (-0.5,1)   \  /   (1,0.5)               <-----|
#           -------------> t (1,0)                     \/
###################################################################################################


@inline function _normal_vector!(v,h::T,normalized_minoraxis::T) where T
  v.=[v[1]/(one(T)+h)^2,v[2]/(normalized_minoraxis+h)^2]
  norm_v=hypot(v...)
  v./=norm_v
end

@inline function _tangent_vector!(v,h::T,normalized_minoraxis::T) where T
  _normal_vector!(v,h,normalized_minoraxis)
  (v[1],v[2])=(v[2],-v[1])
end

@inline _rotation_matrix(θ)= SMatrix{2,2}([cosd(θ) sind(θ);-sind(θ) cosd(θ)] )

@inline function _generate_rays_information_from_orbit(orb::Orbit{T}) where T
  (w,z)=(orb.w,orb.z)
  # LIMB ANGLE WITH RESPECT TO CENTER TO THE EARTH
   θ = atand(z/w)
  (tx,ty)=(z,-w)|> x-> x./hypot(x...) .*INWARD_NORMAL
  #################################

  angle= orb.ang*INWARD_NORMAL
  return (w,z,(_rotation_matrix(angle)*[tx,ty])...)
end

"""
    create_static_rays(orb::Orbit)::SRay2D{T}

Create a static ray from the datum and the orbit. If the datum is not provided, it will use the NormalizedEarth datum
"""
function create_static_rays(orb::Orbit{T})::SRay2D{T} where T
  return SRay2D(_generate_rays_information_from_orbit(orb)...)
end
"""
    create_mutable_rays(orb::Orbit)::MRay2D{T}
Create a mutable ray from the datum and the orbit.
"""
function create_mutable_rays(orb::Orbit{T})::MRay2D{T} where T
  return MRay2D(_generate_rays_information_from_orbit(orb)...)
end




@inline function _scale_earth_by_h(datum::Datum,h::T) where {Datum,T<:Real}
Ellipsoid(
  majoraxis(ellipsoid(datum))+h,
  minoraxis(ellipsoid(datum))+h
)
end

@inline _scale_earth_from_h(h::T) where T<:Real=_scale_earth_by_h(NormalizedEarth,h)


"""
    create_radii_from_θ([datum::Datum=NormalizedEarth],θ::T)::Radius{T} where {Datum,T}

Create a radius for the earth Ellipse that is used to compute the intersection of the Ray2D with the atmosphere using the angle θ and the datum.
If the datum is not provided, it will use the NormalizedEarth datum.

"""
create_radii_from_θ(datum::Datum,θ) where {Datum}=Radius(datum,θ)
create_radii_from_θ(θ) = Radius(NormalizedEarth,θ)

"""
    scale_earth_by_h([datum::Datum=NormalizedEarth],h::T)::Ellipsoid{T} where {Datum,T}

Scale the earth Ellipsoid by h. If the datum is not provided, it will use the NormalizedEarth datum.

"""
scale_earth_by_h(datum::Datum,h::T) where {Datum,T<:Real}=_scale_earth_by_h(datum,h)
scale_earth_by_h(h::T) where T<:Real=_scale_earth_by_h(NormalizedEarth,h)



#########################
"""
    getIntersectionObjects(lla::AbstractArray{L}) where L<:LLA2D{Datum} where Datum

Get the intersection objects from the LLA2D array. The intersection objects are the levels and radii of the earth.
Returns the h_levels, scale_levels, θ_radii, and line_radii
"""
function getIntersectionObjects(lla::AbstractArray{L}) where L<:LLA2D{Datum} where Datum
  (nlevels,nradii)=size(lla)
  h_levels = Vector{Float64}(undef,nlevels)
  scale_levels = Vector{Ellipsoid{Float64}}(undef,nlevels)
  θ_radii  = Vector{Float64}(undef,nradii)
  line_radii = Vector{Radius{Float64}}(undef,nradii)

  for i in axes(lla,1)
    h_levels[i]=lla[i,1].h
    scale_levels[i]=_scale_earth_by_h(Datum,lla[i,1].h)
  end

  for j in axes(lla,2)
    θ_radii[j]=lla[1,j].θ
    line_radii[j]=Radius(Datum,lla[1,j].θ)
  end

  return h_levels,StructArray(scale_levels).scale,SemiCircularArray(θ_radii),SemiCircularArray(line_radii)
end

const NEARBY_LEFT_INDEX=+1
const NEARBY_RIGHT_INDEX=-1
const NEARBY_TOP_INDEX=-1
const NEARBY_BOTTOM_INDEX=+1


## USEFUL CONSTANTS FOR READABILITY AND AVOIDING MISTAKES
const LEFT_RADIUS_INDEX=1
const RIGHT_RADIUS_INDEX=0
const BOTTOM_LEVEL_INDEX=1
const TOP_LEVEL_INDEX=0


@enum Intersection begin
  IntersectionWithLevel
  IntersectionWithRadius
end
######################################################


@inline function _compute_distance_from_radii(ray::Ray2D{T},line_left::Radius{T}, line_right::Radius{T})::T where T
##############################################
  # FIND LEFT AND RIGHT INTERSECTION WITH RADIUS
  ##############################################
  ray_at_1=ray.origin+ray.direction
  t_radius_l=_distance_from_segment(ray.origin,ray_at_1,line_left.PointA,line_left.PointB)
  t_radius_r=_distance_from_segment(ray.origin,ray_at_1,line_right.PointA,line_right.PointB)
  return (t_radius_l,t_radius_r)
end
@inline function _compute_distance_from_ellipses(ray::Ray2D{T},scale_top::ST,scale_bottom::ST)::T where {T,ST}
  ##############################################
  # FIND TOP AND BOTTOM INTERSECTION WITH LEVELS
  ##############################################
  (origin_top,origin_bottom)=(scale_top,scale_bottom).*ray.origin
  (direction_top,direction_bottom)=(scale_top,scale_bottom).*ray.direction
  t_level_b=_distance_from_unit_circle(origin_bottom,direction_bottom)
  t_level_t=_distance_from_unit_circle(origin_top,direction_top)
  return (t_level_b,t_level_t)
end


@inline function _find_real_intersection(ray::Ray2D{T},  #could be both a mutable or a static ray
  input_index::TI,
  refractive_map::M,
  len_h::Int,
  scale_levels::V1,
  line_radii::V2) where {T,TI,M,V1,V2}

  (line_left, line_right) = line_radii[input_index[2]+LEFT_RADIUS_INDEX],line_radii[input_index[2]+RIGHT_RADIUS_INDEX]
  (t_radius_l,t_radius_r)=_compute_distance_from_radii(ray,line_left,line_right)

  ##############################################
  # FIND TOP AND BOTTOM INTERSECTION WITH LEVELS
  ##############################################
  (scale_bottom,scale_top)=(scale_levels[input_index[1]+BOTTOM_LEVEL_INDEX],scale_levels[input_index[1]+TOP_LEVEL_INDEX])
  (t_level_b,t_level_t)=_compute_distance_from_ellipses(ray,scale_top,scale_bottom)
  ##############################################
  # FIND THE MINIMUM DISTANCE
  ###############################################
  t_radius= min(t_radius_l,t_radius_r)
  leftright= t_radius_l<t_radius_r ? OUTWARD_NORMAL : INWARD_NORMAL # left or right
  t_level= min(t_level_b,t_level_t)
  (bottomtop,bottomtop_index)= t_level_b<t_level_t ? (OUTWARD_NORMAL,BOTTOM_LEVEL_INDEX) : (INWARD_NORMAL,TOP_LEVEL_INDEX) # bottom or top

  (t_real, radius_or_level,idx_neighbor)=t_radius<t_level ? (t_radius,IntersectionWithRadius,(0,leftright)) : (t_level,IntersectionWithLevel,(bottomtop,0))

  @debug "idx neighbor: $(idx_neighbor), idx: $(input_index)"
  neighbor_index=input_index.+idx_neighbor

  ### NOTE: Might need to check in test what happens when the index is out of bounds
  # on the opposite side of the earth
  # one idea could be adding 2 additinal indexes to the refractive map 0 and n+1
  # let check current performance first
  ##################################################################################
  n₁=T(1.0);
  n₁=ifelse(0<neighbor_index<len_h,refractive_map[neighbor_index...],n₁);

  return (t_real,n₁,neighbor_index,radius_or_level)
end


"""
    ray_intersection!!(N::MVec2{T},ray::Ray2D{T},n₀::T,input_index::TI,refractive_map::M,h_levels::V1,scale_levels::V2,line_radii::V3,majoraxis_earth::T,minoraxis_earth::T,squared_eccentricity_earth::T) where {T,TI,M,V1,V2,V3}

Compute the intersection of the ray with the atmosphere and modifies in place both the normal vector and the ray.
It returns the effective elevation h, the index of the intersection, and if it is a level or a radius intersection.

# Input
- `N::MVec2{T}`: The normal vector to the intersection point
- `ray::Ray2D{T}`: The ray to be bent
- `n₀::T`: The refractive index of the medium where the ray is coming from
- `input_index::TI`: The wedge index
- `refractive_map::M`: The map of refractive index
- `h_levels::V1`: The levels for i
- `scale_levels::V2`: The scaling map to the unit circle
- `line_radii::V3`: The line Radii
- `majoraxis_earth::T`: The major axis of the earth
- `minoraxis_earth::T`: The minor axis of the earth
- `squared_eccentricity_earth::T`: The squared eccentricity of the earth

# Output
- `h_effective::T`: The effective elevation
- `index::TI`: The index of the intersection
- `radius_or_level::Intersection`: The type of intersection

#In-place
- `N::MVec2{T}`: The normal vector to the intersection point
- `ray::Ray2D{T}`: The ray to be bent
"""
function ray_intersection!!(N::MVec2{T},ray::MRay2D{T}, # ray
  n₀::T, # refractive index
  input_index::TI, # wedge index
  refractive_map::M, # map of refractive index
  h_levels::V1, # levels for i
  scale_levels::V2, # scaling mapping to unit circle
  line_radii::V3, # line Radii
  majoraxis_earth::T,minoraxis_earth::T,squared_eccentricity_earth::T)  where {T<:IEEEFloat,TI,M,V1,V2,V3}

  # the length of h_index is used to check if the index is out of bounds
  len_h= length(h_levels)

  (t_real,n₁,index,radius_or_level)=_find_real_intersection(ray,input_index,refractive_map,len_h,scale_levels,line_radii)

  # mutable ray does not need to be reallocated
  ray.origin.=ray(t_real)

  h_effective,_=_extract_h_from_position(ray.origin...,majoraxis_earth,minoraxis_earth,squared_eccentricity_earth);

  isReflected=false;

  if (radius_or_level==IntersectionWithLevel)
    h_level=h_levels[input_index[1]+bottomtop_index]
    isReflected = _bend_ellipse!!(N,ray,n₀,n₁,h_level,minoraxis_earth,bottomtop)
  end

  if (radius_or_level==IntersectionWithRadius)
    isReflected = _bend_radii!!(N,ray,n₀,n₁,minoraxis_earth,leftright)
  end

  index= isReflected ? input_index : index
  return (h_effective,index,radius_or_level)
end


"""
    ray_intersection!(N::MVec2{T},ray::Ray2D{T},n₀::T,input_index::TI,refractive_map::M,h_levels::V1,scale_levels::V2,line_radii::V3,majoraxis_earth::T,minoraxis_earth::T,squared_eccentricity_earth::T) where {T,TI,M,V1,V2,V3}

Compute the intersection of the ray with the atmosphere and modifies in place the normal vector.
It returns the new ray, the effective elevation h, the index of the intersection, and if it is a level or a radius intersection.

# Input
- `N::MVec2{T}`: The normal vector to the intersection point
- `ray::Ray2D{T}`: The ray to be bent
- `n₀::T`: The refractive index of the medium where the ray is coming from
- `input_index::TI`: The wedge index
- `refractive_map::M`: The map of refractive index
- `h_levels::V1`: The levels for i
- `scale_levels::V2`: The scaling map to the unit circle
- `line_radii::V3`: The line Radii
- `majoraxis_earth::T`: The major axis of the earth
- `minoraxis_earth::T`: The minor axis of the earth
- `squared_eccentricity_earth::T`: The squared eccentricity of the earth

# Output
- `ray_out::Ray2D{T}`: The new ray after the bend
- `h_effective::T`: The effective elevation
- `index::TI`: The index of the intersection
- `radius_or_level::Intersection`: The type of intersection

#In-place
- `N::MVec2{T}`: The normal vector to the intersection point
"""
function ray_intersection!(N::MVec2{T},ray::Ray2D{T}, # ray
  n₀::T, # refractive index
  input_index::TI, # wedge index
  refractive_map::M, # map of refractive index
  h_levels::V1, # levels for i
  scale_levels::V2, # scaling mapping to unit circle
  line_radii::V3 # line Radii
  ) where {T<:IEEEFloat,TI,M,V1,V2,V3}
##

  # the length of h_index is used to check if the index is out of bounds
  len_h= length(h_levels)

  (t_real,n₁,index,radius_or_level)=_find_real_intersection(ray,input_index,refractive_map,len_h,scale_levels,line_radii)


  position=ray(t_real)

  h_effective,_=_extract_h_from_position(position...,majoraxis_earth,minoraxis_earth,squared_eccentricity_earth);

  if isa(radius_or_level,LevelIntersection)
    h_level=h_levels[input_index[1]+bottomtop_index]
    ray_out,isReflected = _bend_ellipse!(N,position,ray.direction,n₀,n₁,h_level,minoraxis_earth,normal_direction)
  else
    ray_out,isReflected = _bend_radii!(N,position,ray.direction,n₀,n₁,minoraxis_earth,normal_direction)
  end
  index= isReflected ? input_index : index
  return (ray_out,h_effective,index,radius_or_level)
end



"""
    ray_intersection(ray::Ray2D{T},n₀::T,input_index::TI,refractive_map::M,h_levels::V1,scale_levels::V2,line_radii::V3,majoraxis_earth::T,minoraxis_earth::T,squared_eccentricity_earth::T) where {T,TI,M,V1,V2,V3}

Compute the intersection of the ray with the atmosphere and returns the new ray, the effective elevation h, the index of the intersection, and if it is a level or a radius intersection.

# Input
- `ray::Ray2D{T}`: The ray to be bent
- `n₀::T`: The refractive index of the medium where the ray is coming from
- `input_index::TI`: The wedge index
- `refractive_map::M`: The map of refractive index
- `h_levels::V1`: The levels for i
- `scale_levels::V2`: The scaling map to the unit circle
- `line_radii::V3`: The line Radii
- `majoraxis_earth::T`: The major axis of the earth
- `minoraxis_earth::T`: The minor axis of the earth
- `squared_eccentricity_earth::T`: The squared eccentricity of the earth

# Output
- `ray_out::Ray2D{T}`: The new ray after the bend
- `h_effective::T`: The effective elevation
- `index::TI`: The index of the intersection
- `radius_or_level::Intersection`: The type of intersection
"""

ray_intersection(ray::Ray2D{T}, # ray
  n₀::T, # refractive index
  input_index::TI, # wedge index
  refractive_map::M, # map of refractive index
  h_levels::V1, # levels for i
  scale_levels::V2, # scaling mapping to unit circle
  line_radii::V3 # line Radii
  ) where {T<:IEEEFloat,TI,M,V1,V2,V3}=ray_intersection!(fill(one(T),MVec2),ray,n₀,input_index,refractive_map,h_levels,scale_levels,line_radii)
