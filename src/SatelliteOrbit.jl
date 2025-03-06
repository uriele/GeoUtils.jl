
const ATOL64 = ScopedValue(1.0e-10)
const ATOL32 = ScopedValue(1.0f-5)
abstract type AbstractRay{T<:IEEEFloat} end
"""
    atol(T)
    atol(x::T)

Absolute tolerance used in algorithms for approximate
comparison with numbers of type `T`. It is used in the
source code in calls to the [`isapprox`](@ref) function:

```julia
isapprox(num1::T, num2::T, atol=atol(T))
```
"""
atol(x) = atol(typeof(x))
atol(::Type{Float64}) = ATOL64[]
atol(::Type{Float32}) = ATOL32[]


const Vec3{T} = SVector{3,T}

const Vec2{T} = SVector{2,T}

const Vec2(x::T,y::T=T(0)) where T<:Number=Vec2{T}(x,y)
const Vec3(x::T,y::T=T(0),z::T=T(0)) where T<:Number=Vec3{T}(x,y,z)

const Vec3F64 = Vec3{Float64}

const Vec2F64 = Vec2{Float64}

const Vec3F32 = Vec3{Float32}

const Vec2F32 = Vec2{Float32}

SDiagonal(v::Vec2{T}) where T = Diagonal(v)
SDiagonal(v::Vec3{T}) where T = Diagonal(v)

struct Ray2D{T}<: AbstractRay{T}
  origin::Vec2{T}
  direction::Vec2{T}
  _normalized::Bool
  @inline function Ray2D{T}(origin::Vec2{T},direction::Vec2{T},_normalized=false) where T
    _hypothenuse(direction)==0 && throw(ArgumentError("Direction cannot be zero"))
    return new{T}(origin,_normalize_direction(direction),_normalized)
  end
end


Ray2D(origin::Vec2{T},direction::Vec2{T}) where T<:IEEEFloat = Ray2D{T}(origin,direction)
Ray2D(x::T,y::T,dx::T,dy::T) where T<:IEEEFloat = Ray2D(Vec2(x,y),Vec2(dx,dy))
Ray2D(x::R1,y::R2,dx::R3,dy::R4) where {R1<:Real,R2<:Real,R3<:Real,R4<:Real} = Ray2D(promote(x,y,dx,dy)...)
Ray2D(x::I,y::I,dx::I,dy::I) where I<:Int = Ray2D(float(x),y,dx,dy)
Ray2D(T::Type=Float64)=Ray2D(T.(0.0),T.(0.0),T.(1.0),T.(1.0))

NormalizedRay2D(origin::Vec2{T},direction::Vec2{T}) where T<:IEEEFloat = Ray2D{T}(origin,direction,true)
NormalizedRay2D(x::T,y::T,dx::T,dy::T) where T<:IEEEFloat = NormalizedRay2D(Vec2(x,y),Vec2(dx,dy))
NormalizedRay2D(x::R1,y::R2,dx::R3,dy::R4) where {R1<:Real,R2<:Real,R3<:Real,R4<:Real} = NormalizedRay2D(promote(x,y,dx,dy)...)
NormalizedRay2D(x::I,y::I,dx::I,dy::I) where I<:Int = NormalizedRay2D(float(x),y,dx,dy)


# Interfaces for Ray implementation
origin(r::AbstractRay)=r.origin
direction(r::AbstractRay)=r.direction
isnormalized(r::AbstractRay)=r._normalized
(r::AbstractRay)(t::R) where R<:Real =r.origin+t.*r.direction

"""
  normalize_datum!([datum=WGS84Latest],sorbit::StructArray{R}) where R<:AbstractRays

Normalize in-place the Ray with respect to the Earth's major axis.

Note: sorbit is a StructArray of SatOrbit (SatOrbit is an immutable struct).
"""
function normalize_datum!(datum::Datum,ray::Ar) where {Datum,Ar<:StructArray{R}} where R<:AbstractRay{T} where T
  majoraxis_earth= majoraxis(ellipsoid(datum)) |> ma-> uconvert(km,ma) |> ustrip
  for i in eachindex(ray)
    ray._normalized[i]==true && continue
    #Staticvector changes debug
    ray.origin[i]=ray.origin[i]./majoraxis_earth
    ray._normalized[i]=true
  end
  return nothing
end
normalize_datum!(ray::Ar)  where {Ar<:StructArray{R}} where R<:AbstractRay{T} where T = normalize_datum!(WGS84Latest,ray)


"""
  normalize_datum([datum=WGS84Latest],ray::R)::R where R<:AbstractRay

Normalize the Ray with respect to the Earth's major axis and return a new Ray.
"""
function normalize_datum(datum::Datum,ray::R)::R where {Datum,R<:AbstractRay{T}} where T
  ray._normalized==true && return ray
  majoraxis_earth= majoraxis(ellipsoid(datum)) |> ma-> uconvert(km,ma) |> ustrip
  _origin=ray.origin/majoraxis_earth
  return R(_origin,ray.direction,true)
end
normalize_datum(ray::R) where R<:AbstractRay{T} where T = normalize_datum(WGS84Latest,ray)


function (lm::LinearMap)(r::Ray2D{T}) where T
  return Ray2D(Vec2(lm(r.origin)...),_normalize_direction(Vec2(lm(r.direction)...)))
end

(r::Ray2D{T})(t::R) where {R<:Real,T<:IEEEFloat} = r.origin+T(t)*r.direction

struct Ellipsoid{T<:IEEEFloat}
  # The map is assumed to be centered at the origin so I do not need to store the center
  # and the affine translation
  radius::Vec2{T}
  scale ::LinearMap
  unscale::LinearMap

  @inline function Ellipsoid(majoraxis_earth::T,minoraxis_earth::T) where T
    radius=Vec2{T}(majoraxis_earth,minoraxis_earth)
    unscale=LinearMap(SDiagonal(radius))
    scale=inv(unscale)
    return new{T}(radius,scale,unscale)
  end
end
Ellipsoid(majoraxis_earth::T) where T = Ellipsoid(majoraxis_earth,majoraxis_earth)

struct Radius{T<:IEEEFloat}
  PointA::Vec2{T}
  PointB::Vec2{T}
end
Radius(ray::Ray2D{T}) where T = Radius{T}(ray(0),ray(1))



get_direction(r::Ray2D)=r.direction
get_origin(r::Ray2D)=r.origin



@inline _dot(v::Vec2{T}) where T = dot(v,v)
@inline _dot(v1::Vec2{T},v2::Vec2{T}) where T = dot(v1,v2)
@inline _dot(v::Vec3{T}) where T = dot(v,v)
@inline _dot(v1::Vec3{T},v2::Vec3{T}) where T = dot(v1,v2)

@inline _fastquadratic(halfB::T,C::T) where T = sqrt_llvm(halfB*halfB-C)
@inline _hypothenuse(v::Vec2{T}) where T = hypot(v.x,v.y)
@inline _normalize_direction(v::Vec2{T}) where T = v/_hypothenuse(v)


"""
    distance_from_unit_circle(origin::Vec2{T},direction::Vec2{T})::T where T
    distance_from_unit_circle(ray::Ray2D{T})::T where T

Calculate the distance from the unit circle to a ray or a point
"""
function distance_from_unit_circle(origin::Vec2{T},direction::Vec2{T})::T where T
  # it is supposed to be divided by the square of the direction, but assuming normalized
  #
  # -1+x^2/majoraxis_earth^2+y^2/minoraxis_earth^2=0
  #
  # C= x0^2/majoraxis_earth^2+y0^2/minoraxis_earth^2-1
  # B= x0*dx/majoraxis_earth^2+y0*dy/minoraxis_earth^2
  # A= dx^2/majoraxis_earth^2+dy^2/minoraxis_earth^2
  #t^2+(mu-B/2)(-mu-B/2)*t+(mu^2-C)=0
  #t^2-(B/2)^2*t+(mu^2-C)=0
  C=origin.x*origin.x+origin.y*origin.y-1
  halfB=-(origin.x*direction.x+origin.y*direction.y)
  A= max(hypot(direction.x,direction.y),1e-10)
  C/=A
  halfB/=A

  μ² =halfB*halfB-C

  if μ²<=0
    return T(Inf)
  end

  μ =sqrt(μ²)
  return filter(x-> x>=0, [μ,-μ].+halfB) |> x-> reduce(min,x;init=Inf)

end

"""
    distance_from_unit_ellipse(origin::Vec2{T},direction::Vec2{T})::T where T
    distance_from_unit_elllipse(ray::Ray2D{T})::T where T

Calculate the distance from the unit circle to a ray or a point
"""
function distance_from_ellipse(origin::Vec2{T},direction::Vec2{T};minoraxis_earth::T=1.0,h::T=0.0)::T where T
  # it is supposed to be divided by the square of the direction, but assuming normalized
  #
  # -1+x^2/majoraxis_earth^2+y^2/minoraxis_earth^2=0
  #
  # C= x0^2/majoraxis_earth^2+y0^2/minoraxis_earth^2-1
  # B= x0*dx/majoraxis_earth^2+y0*dy/minoraxis_earth^2
  # A= dx^2/majoraxis_earth^2+dy^2/minoraxis_earth^2
  #t^2+(mu-B/2)(-mu-B/2)*t+(mu^2-C)=0
  #t^2-(B/2)^2*t+(mu^2-C)=0
  invb²=1/(minoraxis_earth+h) |> x-> x*x
  inva²=1/(1+h) |> x-> x*x

  scale=[inva²; invb²]

  origin=scale.*origin
  direction=scale.*direction

  C=origin.x*origin.x+origin.y*origin.y-1
  halfB=-(origin.x*direction.x+origin.y*direction.y)
  A= max(hypot(direction.x,direction.y),1e-10)
  C/=A
  halfB/=A

  μ² =halfB*halfB-C

  if μ²<=0
    return Inf
  end

  μ =sqrt(μ²)
  return filter(x-> x>=0, [μ,-μ].+halfB) |> x-> reduce(min,x;init=Inf)

end


 # direction it is not necessary and save a division
 # C=_dot(origin)-T(1)  # distance of the ray source from the origin
 # halfB=-_dot(origin,direction) # projection of the origin on the direction
 # μ²=halfB*halfB-C
 # if μ²<=0
 #   return Inf
 # end
 # sqrt_llvm(μ²) |>             # use llvm intrinsic for sqrt for fast math,if negative it will return NaN
 # x -> (halfB*ones(2)+[x,-x]) |>
 # x -> reduce(min,filter(x->x>=0,x);init=Inf)  # filter takes advantage of the fact that NaN return false for >0 <0 or ==0

distance_from_unit_circle(ray::Ray2D{T}) where T = distance_from_unit_circle(ray.origin,_normalize_direction(ray.direction))


"""
  distance_from_radii(ray::Ray2D{T},r::Radius{T})::T where T

Calculate the distance from a ray to one of Earth Radii or from a line with coefficients coeff_a,coeff_b,coeff_c

```math
coeff_a*x+coeff_b*y+coeff_c=0
```
"""
distance_from_radii(ray::Ray2D{T},radius::Radius{T}) where T=
  distance_from_segment(ray::Ray2D{T},radius.PointA,radius.PointB)


"""
    distance_from_segment(A::Vec2{T},B::Vec2{T},C::Vec2{T},D::Vec2{T})::T where T
    distance_from_segment(ray::Ray2D{T},C::Vec2{T},D::Vec2{T})::T where T

Calculate the distance from a segment AB to a segment CD
"""
@inline function distance_from_segment(A::Vec2{T},B::Vec2{T},C::Vec2{T},D::Vec2{T})::T where T
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

@inline distance_from_segment(ray::Ray2D{T},C::Vec2{T},D::Vec2{T}) where T = distance_from_segment(ray(0),ray(1),C,D)



#### Trait to define bend type
abstract type IntersectionStyle end
struct IsIntersection <: IntersectionStyle end
struct NoIntersection <: IntersectionStyle end

IntersectionStyle(x)= IntersectionStyle(typeof(x))
IntersectionStyle(::Type) = NoIntersection()

abstract type IntersectionLocation end
abstract type LevelIntersection <: IntersectionLocation end
abstract type RadiusIntersection  <: IntersectionLocation end
abstract type RadiusLevelIntersection <: IntersectionLocation end
abstract type LevelRadiusIntersection <: IntersectionLocation end

struct TopIntersection <: LevelIntersection end
struct BottomIntersection <: LevelIntersection end
struct LeftIntersection <: RadiusIntersection end
struct RightIntersection <: RadiusIntersection end
struct LeftTopIntersection <: RadiusLevelIntersection end
struct LeftBottomIntersection <: RadiusLevelIntersection end
struct RightTopIntersection <: RadiusLevelIntersection end
struct RightBottomIntersection <: RadiusLevelIntersection end
struct TopLeftIntersection <: LevelRadiusIntersection end
struct TopRightIntersection <: LevelRadiusIntersection end
struct BottomLeftIntersection <: LevelRadiusIntersection end
struct BottomRightIntersection <: LevelRadiusIntersection end

IntersectionStyle(::Type{<:LevelIntersection}) = IsIntersection()
IntersectionStyle(::Type{<:RadiusIntersection}) = IsIntersection()
IntersectionStyle(::Type{<:RadiusLevelIntersection}) = IsIntersection()
IntersectionStyle(::Type{<:LevelRadiusIntersection}) = IsIntersection()


@inline function advance(::L, r::Ray2D{T},scale)::T where {L<:LevelIntersection,T}
  _orig=scale(r.origin)
  _direc=scale(r.direction)
  _hyp_direct= hypot(_direc...)
  return distance_from_unit_circle(_orig,_direc./_hyp_direct) |> x-> x/_hyp_direct
end


@inline function advance(::R, ray::Ray2D{T},radius)::T where {R<:RadiusIntersection,T}
  return distance_from_radii(ray,radius)
end

abstract type AbstractInterface{T<:IEEEFloat} end

struct Interface{T}<:AbstractInterface{T}
h::T
n²::T
n::T
end

Interface(;n::T=1.0,h::T=T(0.0)) where {T<:IEEEFloat} = Interface{T}(h,n*n,n)
#####

@inline bend(::IsIntersection,x::T,ray::Ray2D,t;kwargs...) where T= error("bend function not implemented for $(T)")
@inline bend(::IsIntersection,x::Type{T},ray::Ray2D,t;kwargs...) where T= error("bend function not implemented for $(T)")
@inline bend(x,ray::Ray2D,t;kwargs...) = bend(IntersectionStyle(x),x,ray,t;kwargs...)


const SHIFTORIGIN=1e-12

@inline function _bend_initialize(ray::Ray2D,t::T, n₀::T,n₁::T) where T
  _neworigin=ray(t)
  _direction=ray.direction
  n₀₁=n₀/n₁
  return (_neworigin,_direction,n₀₁,n₀₁*n₀/n₁)
end

@inline function _bend_ellipse(ray::Ray2D,t::T, n₀::T,n₁::T, h::T=T(0), outward=1.0) where T
  squared_eccentricity_earth=eccentricity²(ellipsoid(NormalizedEarth))
  ## shift the origin to avoidnumerical issues

  (_neworigin,_direction,n₀₁,n₀₁²)=_bend_initialize(ray,t,n₀,n₁)
  N=_normal_vector(_neworigin...,h,squared_eccentricity_earth).*outward |> x-> x/hypot(x...)
  ray_out,isReflected=_bend_common(_neworigin,_direction,N,n₀₁,n₀₁²)
  return (ray_out,isReflected)
end


@inline function _bend_radii(ray::Ray2D,t::T, n₀::T,n₁::T, h::T=T(0), outward=1.0) where T
  squared_eccentricity_earth=eccentricity²(ellipsoid(NormalizedEarth))
  (_neworigin,_direction,n₀₁,n₀₁²)=_bend_initialize(ray,t,n₀,n₁)
  N=_tangent_vector(_neworigin...,h,squared_eccentricity_earth).*outward |> x-> x/hypot(x...)
  ray_out,isReflected=_bend_common(_neworigin,_direction,N,n₀₁,n₀₁²)

  return (ray_out,isReflected)
end

const OUTWARD_NORMAL=1.0
const INWARD_NORMAL=-1.0


#TopLevelIntersection -1
#--------------------
#         |  n
#         V
####################
@inline bend(::TopIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)    = _bend_ellipse(ray,t,n₀,n₁,h,INWARD_NORMAL)
#RightRadiusIntersection -1
#
#      n    |
# <---------|
#           |
####################
@inline bend(::RightIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)  = _bend_radii(ray,t,n₀,n₁,0.0,INWARD_NORMAL)


#BottomLevelIntersection 1
#         ^
#         |  n
#--------------------
####################

@inline bend(::BottomIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) = _bend_ellipse(ray,t,n₀,n₁,h,OUTWARD_NORMAL)
#LeftRadiusIntersection 1
# |
# |----------> n
# |
#####################
@inline bend(::LeftIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)   = _bend_radii(ray,t,n₀,n₁,0.0,OUTWARD_NORMAL)

#LeftTop 1 -1
@inline bend(::LeftTopIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)=
  bend(LeftIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(TopIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)
#LeftBottom 1 1
@inline bend(::LeftBottomIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(LeftIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(BottomIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)

#RightTop -1 -1
@inline bend(::RightTopIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(RightIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(TopIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)
#RightBottom -1 1
@inline bend(::RightBottomIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)=
  bend(RightIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(BottomIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)

#TopLeft -1 1
@inline bend(::TopLeftIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(TopIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(LeftIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)
#TopRight -1 1
@inline bend(::TopRightIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(TopIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(RightIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)

#BottomLeft 1 1
@inline bend(::BottomLeftIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(BottomIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(LeftIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)

#BottomRight 1 -1
@inline bend(::BottomRightIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(BottomIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(RightIntersection(),x,0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)


# return both the ray
@inline function _bend_common(origin::Vec2{T},_direction::Vec2{T},
  N::Vec2{T},n₀₁::T,n₀₁²::T) where T
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
###################################################################################################
# VISUAL REPRESENTATION OF TANGENT AND NORMAL VECTORS
###################################################################################################
#          n=(-nx,ny), t=(ny,nx)
#           ^  n                                                                         ^ t(0,1)
#           | (0,1)                                n \    / t                    n(-1,0) |
#           |                              (-0.5,1)   \  /   (1,0.5)               <-----|
#           -------------> t (1,0)                     \/
###################################################################################################
@inline _normal_vector(x0,y0,h,squared_eccentricity_earth)=Vec2(x0/(1+h)^2,y0/(sqrt(1-squared_eccentricity_earth)+h)^2) |> x->x/hypot(x...)
@inline _tangent_vector(args...)=_normal_vector(args...) |> x-> Vec2(x.y,-x.x)

@inline function _normal_vector_new(v::Vec2{T},b²::T)::Vec2{T} where T
   yn=v[2]/b²
  _magnitude=hypot(v[1],yn)
  return Vec2(v[1]/_magnitude,yn/_magnitude)
end

@inline function _tangent_vector_new(v::Vec2{T},b²::T)::Vec2{T} where T
  yn=v[2]/b²
  _magnitude=hypot(yn,v[1])
  return Vec2(yn/_magnitude,-v[1]/_magnitude)
end

# Create Rays from the orbit

"""
    create_rays(orbit::SatOrbit)::Ray2D{T}

Create a ray from the datum and the orbit. If the datum is not provided, it will use the NormalizedEarth datum
"""
function create_rays(orbit::SatOrbit{T}) where T
  (w,z)=(orbit.w,orbit.z)
  # LIMB ANGLE WITH RESPECT TO CENTER TO THE EARTH
  # θ = atan(z/w)
  (tx,ty)=(z,-w)|> x-> x./hypot(x...) .*INWARD_NORMAL
  #################################

  angle= orbit.ang*INWARD_NORMAL

  @inline _rotation_matrix(θ)= SMatrix{2,2}([cosd(θ) sind(θ);-sind(θ) cosd(θ)] )

  return Ray2D{T}(Vec2(w,z),_rotation_matrix(angle)*Vec2(tx,ty),orbit._normalized)
end

###################################################################################
@inline __using_limb_angle(w′,z′)= (z′,-w′)
@inline __using_nadir_angle(w′,z′)= (w′,z′)

@inline function __from_orbits_to_rays(z::T,w::T,angle::T,islimb::Bool)::Vec2{T} where T<:IEEEFloat
  tangent_magnitude=hypot(z,w)
  angle= angle*INWARD_NORMAL
  cosθ=cosd(angle)
  sinθ=sind(angle)

  z′=  z/tangent_magnitude .*INWARD_NORMAL
  w′=  w/tangent_magnitude .*INWARD_NORMAL

  (tx,ty)=islimb ? __using_limb_angle(w′,z′) : __using_nadir_angle(w′,z′)

  dirx=cosθ*tx+sinθ*ty
  diry=-sinθ*tx+cosθ*ty
  return Vec2(dirx,diry)
end

"""
    rays_from_orbit!(rays::StructArray{Ray2D{T}},orbits) where T

Create a ray from the datum and the orbit. Using a preallocated StructArray.
The data type is given by the type of the StructArray.

!!! note
    The rays and the orbits must have the same size.

See also: [`create_bundle_rays`](@ref)
"""
function rays_from_orbit!(rays::Ar,orbits::SA) where {Ar<:StructArray{Ray2D{T}},SA<:AbstractArray{O}} where {O<:SatOrbit{T}} where {T}
  @assert(length(rays)==length(orbits), "The rays and the orbits must have the same size")
  for i in eachindex(orbits)
    rays.direction[i]=__from_orbits_to_rays(orbits[i].z,orbits[i].w,orbits[i].ang,orbits[i]._islimb)
    rays.origin[i]=Vec2(orbits[i].w,orbits[i].z)
    rays._normalized[i]=orbits[i]._normalized
  end
  return nothing
end

"""
    create_bundle_rays([::Type=Float64],orbits)::StructArray{Ray2D{T}}  where T

Create a bundle of rays from the datum and an array of orbits orbit. If no type is provided, it will use Float64.

See also: `rays_from_orbit!`](@ref)
"""
function create_bundle_rays(::Type{T},orbit::Ao) where {T<:IEEEFloat,Ao<:AbstractArray{O}} where {O<:SatOrbit}
  # Preallocation
  #Rays=StructArray(fill(Ray2D(ones(T,4)...),size(orbit)))
  Rays=StructArray{Ray2D{T}}(undef,size(orbit))

  for i in eachindex(orbit)
    # use accessor to write static arrays
    #θ = atan(orbit.z,orbit.w)
    tangent_magnitude=hypot(orbit.z[i],orbit.w[i])
    angle= orbit.ang[i]*INWARD_NORMAL
    cosθ=cosd(angle)
    sinθ=sind(angle)
    tx=   orbit.z[i]/tangent_magnitude .*INWARD_NORMAL
    ty=  -orbit.w[i]/tangent_magnitude .*INWARD_NORMAL
    dirx=cosθ*tx+sinθ*ty
    diry=-sinθ*tx+cosθ*ty
    Rays.origin[i]=Vec2(orbit.w[i],orbit.z[i])
    Rays.direction[i]=Vec2(dirx,diry)
    Rays._normalized[i]=orbit._normalized[i]
    #Staticvector changes debug
    #Rays[i]=Ray2D{T}(Vec2(orbit.w[i],orbit.z[i]),Vec2(dirx,diry),orbit._normalized[i])
    ########################################

  end
  return Rays
end

create_bundle_rays(allorbits::Ao) where {Ao<:AbstractArray{O}} where {O<:SatOrbit{T}} where {T}=create_bundle_rays(T,allorbits)

@inline function _scale_earth_by_h(datum::Datum,h::T) where {Datum,T<:Real}
Ellipsoid(
  majoraxis(ellipsoid(datum))+h,
  minoraxis(ellipsoid(datum))+h
)
end

@inline function _scale_earth_by_h1(datum::Datum,h::T) where {Datum,T<:Real}
Ellipsoid(
  majoraxis(ellipsoid(datum))+h,
  minoraxis(ellipsoid(datum))+h
).scale
end


@inline _scale_earth_from_h(h::T) where T<:Real=_scale_earth_by_h(NormalizedEarth,h)

@inline function _create_radii_from_θ(datum::Datum, θ::T) where {Datum,T<:Real}
  eccentricity_squared= eccentricity²(ellipsoid(datum))
  origin=convert(ECEF2D{datum},LLA2D{datum}(0.,θ)) |> x-> Vec2([x.w,x.z])
  direction=_normal_vector(origin...,0,eccentricity_squared)
  return Radius(Ray2D(origin,direction))
end
@inline _create_radii_from_θ(θ::T) where T<:Real=_create_radii_from_θ(NormalizedEarth,θ)

create_radii_from_θ(datum::Datum,θ::T) where {Datum,T<:Real}=_create_radii_from_θ(datum,θ)
create_radii_from_θ(θ::T) where T<:Real=_create_radii_from_θ(θ)
scale_earth_by_h(datum::Datum,h::T) where {Datum,T<:Real}=_scale_earth_by_h(datum,h)
scale_earth_by_h(h::T) where T<:Real=_scale_earth_by_h(NormalizedEarth,h)

#########################
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
    line_radii[j]=_create_radii_from_θ(Datum,lla[1,j].θ)
  end

  return h_levels,StructArray(scale_levels).scale,SemiCircularArray(θ_radii),SemiCircularArray(line_radii)
end
#
function getIntersectionObjects(::Datum,h::V1,θ::V2) where {Datum,V1<:AbstractVector{T},V2<:AbstractVector{T}} where T

   scale_levels = Vector{LinearMap}(undef,length(h))
   line_radii   = Vector{Radius{T}}(undef,length(θ))
   map!(x->_scale_earth_by_h1(Datum,x) ,scale_levels,h)
   map!(x->create_radii_from_θ(Datum,x),line_radii,θ)

  return scale_levels,SemiCircularArray(line_radii)
end


const NEARBY_LEFT_INDEX=+1
const NEARBY_RIGHT_INDEX=-1
const NEARBY_TOP_INDEX=-1
const NEARBY_BOTTOM_INDEX=+1
function _intersection_type_radii(lr::LR,position,wedge_index,h_levels,refractive_map) where LR<:RadiusIntersection

  # left is to the left +1 and right is to the right -1
  position_interface_index2= lr==LeftIntersection() ? wedge_index[2]+NEARBY_LEFT_INDEX : wedge_index[2]+NEARBY_RIGHT_INDEX

  index=(wedge_index[1],position_interface_index2)

  h=convert(LLA2D{NormalizedEarth},ECEF2D{NormalizedEarth}(position...)) |> x-> x.h
  n₁=refractive_map[index...]
  # if the last radius was ABOVE the maximum level of stratification,
  # then the ray is in free space and you stop the computation
  # TO DO maybe find a better breaking condition
  # like checking the altitude of the ray outside instead of basing it on
  # the index

  index=ifelse(h<=h_levels[1],index,(0,index[2]))
  (n₁,index,h)
end




function _intersection_type_levels(tb::TB,wedge_index,h_levels,refractive_map) where TB<:LevelIntersection
  ##############################################
  # The -1 is because the number of wedges is one less the number of levels
  #
  #     ------------- h[1]
  #         n[1]
  #     ------------- h[2]
  #         n[2]
  #     ------------- h[3]
  ##############################################
  nlevels=length(h_levels)-1
  ##############################################
  # bottom is to the top +1 and top is to the top -1
  # it is confusing but the reason is that the altitude is ordered from the top to the bottom

  position_interface_index1= tb==BottomIntersection() ? wedge_index[1]+NEARBY_BOTTOM_INDEX : wedge_index[1]+NEARBY_TOP_INDEX
  # this is important, since the altitude uses the h_levels
  # and the h levels are one less than the number of levels
  # the first index coincide with the highest wedge altitude TOP
  # and i+1 index coincide with the lowest wedge altitude BOTTOM
  # since the index uses was i+1 for bottom and i-1 for TOP
  # i need to ADD +1 to get the correct altitude for the TOP
  # and just use the index for the top

  index=(position_interface_index1,wedge_index[2])

  h= tb==BottomIntersection() ? h_levels[position_interface_index1] : h_levels[position_interface_index1+1]

  # the altitude is not circular and it is used to stop the computation of ray tracing
  # thus, the extremes are not considered
  if position_interface_index1<1  # free space
    n₁=1.0
  elseif position_interface_index1>nlevels  # ground full reflection
    n₁=0.01
  else            # atmosphere
    n₁=refractive_map[index...]
  end

   (n₁,index,h)
end

const DEBUG_INTERSECTION=Ref{Int}(0)
setDebugIntersection(x::Int)=DEBUG_INTERSECTION[]=x
setDebugIntersection()=DEBUG_INTERSECTION[]=0
getDebugIntersection()=DEBUG_INTERSECTION[]

## USEFUL CONSTANTS FOR READABILITY AND AVOIDING MISTAKES
const LEFT_RADIUS_INDEX=1
const RIGHT_RADIUS_INDEX=0
const BOTTOM_LEVEL_INDEX=1
const TOP_LEVEL_INDEX=0
######################################################
@inline function new_intersection(ray::Ray2D, # ray
  input_index; # wedge index
  refractive_map, # map of refractive index
  h_levels, # levels for i
  θ_radii, # radii angles for j
  scale_levels, # scaling mapping to unit circle
  line_radii, # line radii
  tangent_quote, # tangent quote
  register,
  #for debug
  previous_intersection=nothing) # register

  ###########################
  n₀= refractive_map[input_index...]
  ###########################

  t_radius_l=advance(LeftIntersection(),ray,line_radii[input_index[2]+LEFT_RADIUS_INDEX])
  t_radius_r=advance(RightIntersection(),ray,line_radii[input_index[2]+RIGHT_RADIUS_INDEX])

  t_level_b=advance(BottomIntersection(),ray,scale_levels[input_index[1]+BOTTOM_LEVEL_INDEX])
  t_level_t=advance(TopIntersection(),ray,scale_levels[input_index[1]+TOP_LEVEL_INDEX])

  t_radius= min(t_radius_l,t_radius_r)
  leftright= t_radius_l<t_radius_r ? LeftIntersection() : RightIntersection() # left or right
  t_level= min(t_level_b,t_level_t)

  bottomtop= t_level_b<t_level_t ? BottomIntersection() : TopIntersection() # bottom or top


  radius_or_level=t_radius<t_level ? leftright : bottomtop
  t_real= min(t_radius,t_level)

  n₀=refractive_map[input_index...]

  if isa(radius_or_level,LevelIntersection)
    n₁,index,h=_intersection_type_levels(bottomtop,input_index,h_levels,refractive_map)
    normal_direction=(radius_or_level==TopIntersection() ? INWARD_NORMAL : OUTWARD_NORMAL)
    ray_out,isReflected = _bend_ellipse(ray,t_real,n₀,n₁,h,normal_direction)
  else
    n₁,index,h=_intersection_type_radii(leftright,ray(t_real),input_index,h_levels,refractive_map)
    normal_direction=(radius_or_level==RightIntersection() ? INWARD_NORMAL : OUTWARD_NORMAL)
    ray_out,isReflected = _bend_radii(ray,t_real,n₀,n₁,h,normal_direction)
  end
  index= isReflected ? input_index : index
  # condition to stop ray tracing
  #target = 1<index[1]<(length(h_levels)-1)
  target = false
  return (ray_out,target,h,index,radius_or_level)
end



function minimizing_function_distance_h²_from_surface!(
  t::T,θ::T,origin::V2,
  direction::V3,b::T,h²::T)::Tuple{Vec2{T},T,T}  where {
    V2<:AbstractVector{T},
    V3<:AbstractVector{T}
  } where T<:Real

  cosθ=cos(θ)
  sinθ=sin(θ)
  bsinθ=b*sinθ
  bcosθ=b*cosθ
  #b²=b*b
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
  #Dᵣ² = direction[1]*direction[1]+direction[2]*direction[2]
  Tₑ² = bcosθ*bcosθ+sinθ*sinθ
  Dᵣ_Tₑ = direction[1]*sinθ-direction[2]*bcosθ
  AB_A=X*cosθ+Y*bsinθ
  eight_AB_Dᵣ_AB_Tₑ=8*AB_Dᵣ*AB_Tₑ

  g1=2*AB_Dᵣ*(twice_AB²_minus_h²)
  g2=2*AB_Tₑ*(1+twice_AB²_minus_h²)

  #∇²f_11= 2*(Dᵣ²)*twice_AB²_minus_h² + 8*(AB_Dᵣ*AB_Dᵣ)
  ∇²f_11= 2*twice_AB²_minus_h² + 8*(AB_Dᵣ*AB_Dᵣ)
  ∇²f_22= 2*(AB_A+Tₑ²)*twice_AB²_minus_h² + 8*(AB_Tₑ*AB_Tₑ)+2*(AB_A+Tₑ²)
  ∇²f_12= Dᵣ_Tₑ*twice_AB²_minus_h² + eight_AB_Dᵣ_AB_Tₑ
  Denom=∇²f_11*∇²f_22-∇²f_12*∇²f_12
  Denom+= (Denom<atol(T))*atol(T)




  return (Vec2(( ∇²f_22*g1-∇²f_12*g2)./Denom,(-∇²f_12*g1+∇²f_11*g2)./Denom),
    distance²+AB²  ,h²+(distance²>atol(T))*distance²)

end



@inline __h_t(x)=max(-x,0)
@inline __h_theta_min(x,θmin)= max(θmin-x,0)
@inline __h_theta_max(x,θmax)= max(x-θmax,0)


function optimize_distance_h²_from_surface_lagrangian2!(x::V,origin::V2,direction::V2,b::T,h²::T,θmin::T,θmax::T;rho=T(1.0),gamma=T(1.5))::T where {V<:AbstractArray{T},V2<:Vec2{T}} where T
    _old=0.0
    _h2=999.0
    kmax=10
    k=0
    @debug "θ_0: $(x[2]) θmin: $θmin θmax: $θmax"
    while k<kmax
      k+=1
      while true
        _fun,_h2=minimizing_function_distance_h²_from_surface_constrained_lagrangian!(x,origin,direction,b,h²,θmin,θmax,rho)
        err=abs(_fun-_old)
        _old=_fun
        if err<10^-15
          #@debug "err: $err h: $(sqrt(_h2))"
          #@debug "b: $(minoraxis), h2: $h2 x: $x"
          break
        end
      end
      rho*=gamma

      if (x[1]>=0 && x[2]>=θmin && x[2]<=θmax)
        break
      end
      if k>kmax
        break
      end
    end

    return _h2
end

function minimizing_function_distance_h²_from_surface_lagrangian!(
  x::A1,origin::A2,
  direction::A3,b::T,h²::T,θmin::T,θmax::T,rho::T)::Tuple{T,T}  where {
    A1<:AbstractArray,
    A2<:Vec2{T},
    A3<:Vec2{T}
  } where T

  cosθ=cos(x[2])
  sinθ=sin(x[2])
  bsinθ=b*sinθ
  bcosθ=b*cosθ
  #b²=b*b

  X=origin[1] + x[1]*direction[1]-cosθ
  Y=origin[2] + x[1]*direction[2]-bsinθ

  AB² = X*X + Y*Y
  distance=AB²-h²
  distance² =distance*distance

  twice_AB²_minus_h² = 2*(AB²-h²)
  AB_Dᵣ = X*direction[1]+Y*direction[2]
  AB_Tₑ = X*sinθ-Y*bcosθ
  #Dᵣ² = direction[1]*direction[1]+direction[2]*direction[2]
  Tₑ² = bcosθ*bcosθ+sinθ*sinθ
  Dᵣ_Tₑ = direction[1]*sinθ-direction[2]*bcosθ
  AB_A=X*cosθ+Y*bsinθ
  eight_AB_Dᵣ_AB_Tₑ=8*AB_Dᵣ*AB_Tₑ


  _ht=__h_t(x[1])
  _hθ_min=__h_theta_min(x[2],θmin)
  _hθ_max=__h_theta_max(x[2],θmax)

  _dht=x[1]<0
  _dhθ_min=(x[2]<θmin)
  _dhθ_max=(x[2]>θmax)

  penalty=0.5*rho*(_ht*_ht+_hθ_min*_hθ_min+_hθ_max*_hθ_max)#-lambda[1]*_ht #-lambda[2]*_hθ
  ∂penalty_t= rho*_ht*_dht
  ∂penalty_θ= rho*(_hθ_min*_dhθ_min+_hθ_max*_dhθ_max)

  fvalue=distance²+AB²+penalty

  g1=2*AB_Dᵣ*(twice_AB²_minus_h²)+∂penalty_t #-lambda[1]*_dht
  g2=2*AB_Tₑ*(1+twice_AB²_minus_h²)+∂penalty_θ #-lambda[2]*_dhθ

  #∇²f_11= 2*(Dᵣ²)*twice_AB²_minus_h² + 8*(AB_Dᵣ*AB_Dᵣ)
  ∇²f_11= 2*twice_AB²_minus_h² + 8*(AB_Dᵣ*AB_Dᵣ)+rho*_dht
  ∇²f_22= 2*(AB_A+Tₑ²)*twice_AB²_minus_h² + 8*(AB_Tₑ*AB_Tₑ)+2*(AB_A+Tₑ²)+rho*(_dhθ_min+_dhθ_max)
  ∇²f_12= Dᵣ_Tₑ*twice_AB²_minus_h² + eight_AB_Dᵣ_AB_Tₑ
  Denom=∇²f_11*∇²f_22-∇²f_12*∇²f_12
  Denom+= (Denom<atol(T))*atol(T)


  x[1]-=( ∇²f_22*g1-∇²f_12*g2)/Denom
  x[2]-=(-∇²f_12*g1+∇²f_11*g2)/Denom

  return (fvalue  ,h²+(distance²>atol(T))*distance²)

end


# Optimization internal function for t>0 and θmin<θ<θmax
@inline _Δc(higher::T,lower::T) where T=higher-lower


function optimization_augmented_lagrangian!(g::AV,x::AV,r_origin::V2,r_direction::V2,b::T,h²::T,θmin::T,θmax::T,slack::AS;
  rho=T(1.0),gamma=T(1.5),kmax_out::Int=10,kmax_in::Int=10)::T where {AV<:AbstractArray{T},AS<:AbstractArray{T},V2<:Vec2{T}} where T

  func=T(0.0)
  _h²=T(999.0)
  kout=0
  # Initialized Lagrange multipliers
  λ1=T(0.0)
  λ2=T(0.0)
  λ3=T(0.0)
  # initialized slack and x variables for t since it is the only one that is reset

  slack[1]=T(0.0)
  x[1]=T(0.0)
#=
  while kout<kmax_out
=#
    func,_h²=inner_optimization!(g,x,r_origin,r_direction,b,h²,θmin,θmax,rho,λ1,λ2,λ3,slack;kmax=kmax_in)
    # not necessary but for commodity
    #=
    (s1,s2,s3)=Vec3(slack[1],slack[2],slack[3])
    (t,θ)=Vec2(x[1],x[2])
    h1=_Δc(t,T(0))-s1
    h2=_Δc(θ,θmin)-s2
    h3=_Δc(θmax,θ)-s3
    # update Lagrange multipliers
    (λ1,λ2,λ3)=Vec3(λ1-h1*rho,λ2-h2*rho,λ3-h3*rho)
    rho*=gamma
    if norm((h1,h2,h3))<atol(T)
      break
    end
    kout+=1
  end
  =#
  return _h²
end


function inner_optimization!(g::AV,x::AV,r_origin::V2,r_direction::V2,b::T,h²::T,θmin::T,θmax::T,rho::T,λ1::T,λ2::T,λ3::T,slack::AS;kmax::Int=30)::Tuple{T,T} where {AV<:AbstractVector{T},AS<:AbstractVector{T},V2<:Vec2{T}} where T
  k=0
  _h²=999.0
  func=0.0
  k=0
  while k<=kmax
    func,_h²=ray_constrained_lagrangian!(g,x,r_origin,r_direction,b,h²,θmin,θmax,rho,
    λ1,λ2,λ3,
    slack)

    if norm(g)<atol(T)
      @info k
     s @info norm(g)
      break
    end
    k+=1
  end
  return func,_h²
end


#=
function ray_constrained_lagrangian!(g::AV,
  x::AV,r_origin::V2,
  r_direction::V2,b::T,h²::T,θmin::T,θmax::T,p_barrier::T,
  λ1::T,λ2::T,λ3::T,
  slack::AS)::Tuple{T,T}  where {
    AV<:AbstractArray{T},
    AS<:AbstractArray{T},
    V2<:Vec2{T}
  } where T

  inv_p_barrier=1/p_barrier
  cosθ=cos(x[2])
  sinθ=sin(x[2])

  bsinθ=b*sinθ
  bcosθ=b*cosθ
  b²=b*b

  X=r_origin[1] + x[1]*r_direction[1]-cosθ
  Y=r_origin[2] + x[1]*r_direction[2]-bsinθ

  AB² = X*X + Y*Y
  distance=AB²-h²
  distance² =distance*distance


  AB_Dᵣ = X*r_direction[1]+Y*r_direction[2]
  AB_Tₑ = X*sinθ-Y*bcosθ
  AB_Tₑ² = AB_Tₑ*AB_Tₑ
  AB_P  = X*cosθ+Y*bsinθ
  Tₑ² = bcosθ*bcosθ+sinθ*sinθ
  Dᵣ_Tₑ = r_direction[1]*sinθ-r_direction[2]*bcosθ
  four_dist=4*distance


  ∂f∂t=four_dist*AB_Dᵣ
  ∂f∂θ=four_dist*AB_Tₑ

  ∇²f_tt= four_dist+8*(AB_Dᵣ*AB_Dᵣ)

  ∇²f_θθ= (AB_P+Tₑ²)*four_dist+ 8*(AB_Tₑ²)

  ∇²f_tθ= (Dᵣ_Tₑ)*four_dist+8*AB_Dᵣ*AB_Tₑ


  Denom = ∇²f_tt*∇²f_θθ - (∇²f_tθ*∇²f_tθ)
  Denom += (Denom<atol(T))*atol(T)
  N1= (∂f∂t*∇²f_θθ-∂f∂θ*∇²f_tθ)/Denom
  N2= (∂f∂θ*∇²f_tt-∂f∂t*∇²f_tθ)/Denom

  x[1]-=N1
  x[2]-=N2

  s1=slack[1]
  s2=slack[2]
  s3=slack[3]


  c1=_Δc(x[1],T(0))
  c2=_Δc(x[2],θmin)
  c3=_Δc(θmax,x[2])

  h1 = c1-s1
  h2 = c2-s2
  h3 = c3-s3

  g[3]=+λ1-p_barrier*h1
  g[4]=+λ2-p_barrier*h2
  g[5]=+λ3-p_barrier*h3

  g[1]=∂f∂t-g[3]
  g[2]=∂f∂θ-g[4]
  g[3]=∂f∂θ+g[5]

  ds1=+N1+g[3]*inv_p_barrier
  ds2=+N2+g[4]*inv_p_barrier
  ds3=-N2+g[5]*inv_p_barrier

  # Projection into the positive orthant
  @info "N1: $N1 N2: $N2"
  @info "λ1: $λ1 λ2: $λ2 λ3: $λ3"
  @info "ds1: $ds1 ds2: $ds2 ds3: $ds3"
  @info "g3: $(g[3]) g4: $(g[4]) g5: $(g[5])"
  @info "s1: $s1 s2: $s2 s3: $s3"
  @info "c1: $c1 c2: $c2 c3: $c3"
  @info "h1: $h1 h2: $h2 h3: $h3"
  #= Try explicitly impose it
  slack[1]=max(s1-ds1,0)
  slack[2]=max(s2-ds2,0)
  slack[3]=max(s3-ds3,0)

  =#
  # Recalculate the Lagrange multipliers
  h1=_Δc(x[1],T(0))-slack[1]
  h2=_Δc(x[2],θmin)-slack[2]
  h3=_Δc(θmax,x[2])-slack[3]


  penalty=p_barrier/2*(h1*h1+h2*h2+h3*h3)-λ1*h1-λ2*h2-λ3*h3

  @info "penality: $penalty θmin: $θmin θmax: $θmax"
  # Compute new function
  fvalue=distance²+ penalty


  @info "g: $g"
  @info "x: $x s: $slack"
  return (fvalue  ,h²+(distance²>atol(T))*distance²)

end
=#




function ray_constrained_lagrangian!(g::AV,
  x::AV,r_origin::V2,
  r_direction::V2,b::T,h²::T,θmin::T,θmax::T,p_barrier::T,
  λ1::T,λ2::T,λ3::T,
  slack::AS)::Tuple{T,T}  where {
    AV<:AbstractArray{T},
    AS<:AbstractArray{T},
    V2<:Vec2{T}
  } where T

  inv_p_barrier=1/p_barrier
  cosθ=cos(x[2])
  sinθ=sin(x[2])

  bsinθ=b*sinθ
  bcosθ=b*cosθ
  b²=b*b

  X=r_origin[1] + x[1]*r_direction[1]-cosθ
  Y=r_origin[2] + x[1]*r_direction[2]-bsinθ

  AB² = X*X + Y*Y
  distance=AB²-h²
  distance² =distance*distance


  AB_Dᵣ = X*r_direction[1]+Y*r_direction[2]
  AB_Tₑ = X*sinθ-Y*bcosθ
  AB_Tₑ² = AB_Tₑ*AB_Tₑ
  AB_P  = X*cosθ+Y*bsinθ
  Tₑ² = bcosθ*bcosθ+sinθ*sinθ
  Dᵣ_Tₑ = r_direction[1]*sinθ-r_direction[2]*bcosθ
  four_dist=4*distance


  ∂f∂t=four_dist*AB_Dᵣ
  ∂f∂θ=four_dist*AB_Tₑ

  ∇²f_tt= four_dist+8*(AB_Dᵣ*AB_Dᵣ)

  ∇²f_θθ= (AB_P+Tₑ²)*four_dist+ 8*(AB_Tₑ²)

  ∇²f_tθ= (Dᵣ_Tₑ)*four_dist+8*AB_Dᵣ*AB_Tₑ


  Denom = ∇²f_tt*∇²f_θθ - (∇²f_tθ*∇²f_tθ)

  #Denom = (abs(Denom)<atol(T))*atol(T)+(abs(Denom)<atol(T))*Denom
  N1= (∂f∂t*∇²f_θθ-∂f∂θ*∇²f_tθ)/Denom
  N2= (∂f∂θ*∇²f_tt-∂f∂t*∇²f_tθ)/Denom

  @info "Denom: $Denom N1: $N1 N2: $N2"
  x[2]-=N2
  x[1]-=N1

  g[1]=∂f∂t
  g[2]=∂f∂θ

#=
  if θmin<θnew<θmax
    x[1]-=N1
    x[2]=θnew
  else
    x[1]-=∂f∂t/∇²f_tt
    x[2]=min(max(θmin,θnew),θmax)
    g[2]=0.0
  end
=#
  return (distance²  ,h²+(distance²>atol(T))*distance²)

end
