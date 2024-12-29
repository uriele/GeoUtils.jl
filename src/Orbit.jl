
const ATOL64 = ScopedValue(1.0e-10)
const ATOL32 = ScopedValue(1.0f-5)

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

struct Ray2D{T<:IEEEFloat}
  origin::Vec2{T}
  direction::Vec2{T}

  @inline function Ray2D(origin::Vec2{T},direction::Vec2{T}) where T
    _hypothenuse(direction)==0 && throw(ArgumentError("Direction cannot be zero"))
    return new{T}(origin,_normalize(direction))
  end
end

function (lm::LinearMap)(r::Ray2D{T}) where T
  return Ray2D(Vec2(lm(r.origin)...),_normalize(Vec2(lm(r.direction)...)))
end

(r::Ray2D{T})(t) where T = r.origin+T(t)*_normalize(r.direction)

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
@inline _normalize(v::Vec2{T}) where T = v/_hypothenuse(v)


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

distance_from_unit_circle(ray::Ray2D{T}) where T = distance_from_unit_circle(ray.origin,_normalize(ray.direction))


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
  #@debug "A=$A B=$B C=$C D=$D"
  M=[(B-A) (C-D)];
  y=(C-A);
  #@debug "M=$M y=$y"
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
  @debug "advance level"
  @debug "scale: $scale"
  @debug "r: $r"
  _orig=scale(r.origin)
  _direc=scale(r.direction)
  _hyp_direct= hypot(_direc...)
  t=distance_from_unit_circle(_orig,_direc./_hyp_direct) |> x-> x/_hyp_direct
  #@debug "t_level: $t"
  #t<=1e-10 ? T(Inf) : t
  return t
end


@inline function advance(::R, ray::Ray2D{T},radius)::T where {R<:RadiusIntersection,T}
  t=distance_from_radii(ray,radius)
  @debug "t_radius: $t, radius: $radius,ray: $ray"
  #t<=1e-10 ? T(Inf) : t
  return t
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
  squared_eccentricity_earth=eccentricity²(ellipsoid(NormalizeEarth))
  ## shift the origin to avoidnumerical issues

  (_neworigin,_direction,n₀₁,n₀₁²)=_bend_initialize(ray,t,n₀,n₁)
  N=_normal_vector(_neworigin...,h,squared_eccentricity_earth).*outward |> x-> x/hypot(x...)
  ray_out,isReflected=_bend_common(_neworigin,_direction,N,n₀₁,n₀₁²)
  @debug "h: $h, n₀: $n₀, n₁: $n₁"
  @debug "N: $(hypot(N...)) $(atand(-N[2]/N[1]))"
  @debug "direction_new: $(hypot(_direction...)) θ: $(atand(_direction[2]/_direction[1])),"
  @debug "direction_old: $(hypot(ray.direction...)) θ: $(atand(ray.direction[2]/ray.direction[1])),"

  sinθ₀= sqrt(1-(-N⋅_direction)^2)
  sinθ₁= sqrt(1-(-N⋅get_direction(ray_out))^2)
  @debug  "sinθ₀= $(sinθ₀), sinθ₁= $(sinθ₁), sinθ₀-sinθ₁=$(sinθ₀-sinθ₁)"
  @debug  "θ₀=$(asind(sinθ₀)), θ₁= $(asind(sinθ₁))"
  @debug "n₀= $(n₀), n₁= $(n₁)"
  isRising= acosd(-N⋅_direction)>90
  @debug "isRising: $isRising, isReflected: $isReflected"
  return (ray_out,isReflected,isRising)
end


@inline function _bend_radii(ray::Ray2D,t::T, n₀::T,n₁::T, h::T=T(0), outward=1.0) where T
  squared_eccentricity_earth=eccentricity²(ellipsoid(NormalizeEarth))
  (_neworigin,_direction,n₀₁,n₀₁²)=_bend_initialize(ray,t,n₀,n₁)
  N=_tangent_vector(_neworigin...,h,squared_eccentricity_earth).*outward |> x-> x/hypot(x...)
  ray_out,isReflected=_bend_common(_neworigin,_direction,N,n₀₁,n₀₁²)
  isRising= acosd(-N⋅_direction) #|> x-> ifelse(abs(x)==90,x,rem(x,180,RoundNearest) )<0
  @debug "N: $(hypot(N...)) $(atand(-N[2]/N[1])), direction: $(hypot(_direction...)) $(atand(-_direction[2]/_direction[1])),"
  @debug "ray_out: $ray_out"
  @debug "direction_new: $(hypot(_direction...)) θ: $(atand(_direction[2]/_direction[1])),"
  @debug "direction_old: $(hypot(ray.direction...)) θ: $(atand(ray.direction[2]/ray.direction[1])),"
  sinθ₀= sqrt(1-(-N⋅_direction)^2)
  sinθ₁= sqrt(1-(-N⋅get_direction(ray_out))^2)
  @debug  "sinθ₀= $(sinθ₀), sinθ₁= $(sinθ₁), sinθ₀-sinθ₁=$(sinθ₀-sinθ₁)"
  @debug  "θ₀=$(asind(sinθ₀)), θ₁= $(asind(sinθ₁))"
  @debug "isRising: $isRising, isReflected: $isReflected"
  return (ray_out,isReflected,isRising)
end

const OUTWARD_NORMAL=1.0
const INWARD_NORMAL=-1.0


#TopLevelIntersection -1
#--------------------
#         |  n
#         V
####################
@inline bend(::TopIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)    = _bend_ellipse(ray,t,n₀,n₁,h,INWARD_NORMAL) |> x-> (x[1],x[2])
#RightRadiusIntersection -1
#
#      n    |
# <---------|
#           |
####################
@inline bend(::RightIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)  = _bend_radii(ray,t,n₀,n₁,0.0,INWARD_NORMAL)  |> x-> (x[1],x[2])


#BottomLevelIntersection 1
#         ^
#         |  n
#--------------------
####################

@inline bend(::BottomIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) = _bend_ellipse(ray,t,n₀,n₁,h,OUTWARD_NORMAL)  |> x-> (x[1],x[2])
#LeftRadiusIntersection 1
# |
# |----------> n
# |
#####################
@inline bend(::LeftIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)   = _bend_radii(ray,t,n₀,n₁,0.0,OUTWARD_NORMAL)  |> x-> (x[1],x[2])

#LeftTop 1 -1
@inline bend(::LeftTopIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)=
  bend(LeftIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(TopIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)  |> x-> (x[1],x[2])
#LeftBottom 1 1
@inline bend(::LeftBottomIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(LeftIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(BottomIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h) |> x-> (x[1],x[2])

#RightTop -1 -1
@inline bend(::RightTopIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(RightIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(TopIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h) |> x-> (x[1],x[2])
#RightBottom -1 1
@inline bend(::RightBottomIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)=
  bend(RightIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(BottomIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h) |> x-> (x[1],x[2])

#TopLeft -1 1
@inline bend(::TopLeftIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(TopIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(LeftIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h)|> x-> (x[1],x[2])
#TopRight -1 1
@inline bend(::TopRightIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(TopIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(RightIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h) |> x-> (x[1],x[2])

#BottomLeft 1 1
@inline bend(::BottomLeftIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(BottomIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(LeftIntersection(),x[1],0;n₀=n₁,n₁=n₂,n₂=n₂,h=h) |> (x[1],x[2])

#BottomRight 1 -1
@inline bend(::BottomRightIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) =
  bend(BottomIntersection(),ray,t;n₀=n₀,n₁=n₁,n₂=n₂,h=h) |>
  x-> bend(RightIntersection(),x,0;n₀=n₁,n₁=n₂,n₂=n₂,h=h) |> (x[1],x[2])


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
# Create Rays from the orbit

"""
    create_rays(datum::Datum,orb::Orbit)::Ray2D{T}
    create_rays(orb::Orbit)::Ray2D{T}

Create a ray from the datum and the orbit. If the datum is not provided, it will use the NormalizeEarth datum
"""
function create_rays(datum::Datum,orb::Orbit) where Datum
  squared_eccentricity_earth=eccentricity²(ellipsoid(datum))
  (w,z)=(orb.w,orb.z)
  #################MARCO CHANGE BACK
  #= LIMB ANGLE WITH RESPECT TO NORMAL TO EARTH
  θ,h=geocentric_to_geodesic_θ(w,z)
  (w_real,z_real)=convert(ECEF2D{datum},LLA2D{datum}(0.,θ)) |> x-> (x.w,x.z)
  #@info "θ: $θ, h: $h"
  (tx,ty)=_tangent_vector(w_real,z_real,0,squared_eccentricity_earth)*INWARD_NORMAL
  =#########################
  # LIMB ANGLE WITH RESPECT TO CENTER TO THE EARTH
   θ = atan(z/w)
  (tx,ty)=(z,-w)|> x-> x./hypot(x...) .*INWARD_NORMAL
  #################################

  angle= orb.ang*INWARD_NORMAL

  @inline _rotation_matrix(θ)= SMatrix{2,2}([cosd(θ) sind(θ);-sind(θ) cosd(θ)] )
  return Ray2D(Vec2(w,z),_rotation_matrix(angle)*Vec2(tx,ty))
end

create_rays(orb::Orbit)=create_rays(NormalizeEarth,orb)


@inline function _scale_earth_by_h(datum::Datum,h::T) where {Datum,T<:Real}
Ellipsoid(
  majoraxis(ellipsoid(datum))+h,
  minoraxis(ellipsoid(datum))+h
)
end

@inline _scale_earth_from_h(h::T) where T<:Real=_scale_earth_by_h(NormalizeEarth,h)

@inline function _create_radii_from_θ(datum::Datum, θ::T) where {Datum,T<:Real}
  eccentricity_squared= eccentricity²(ellipsoid(datum))
  origin=convert(ECEF2D{datum},LLA2D{datum}(0.,θ)) |> x-> Vec2([x.w,x.z])
  direction=_normal_vector(origin...,0,eccentricity_squared)
  return Radius(Ray2D(origin,direction))
end
@inline _create_radii_from_θ(θ::T) where T<:Real=_create_radii_from_θ(NormalizeEarth,θ)

create_radii_from_θ(datum::Datum,θ::T) where {Datum,T<:Real}=_create_radii_from_θ(datum,θ)
create_radii_from_θ(θ::T) where T<:Real=_create_radii_from_θ(θ)
scale_earth_by_h(datum::Datum,h::T) where {Datum,T<:Real}=_scale_earth_by_h(datum,h)
scale_earth_by_h(h::T) where T<:Real=_scale_earth_by_h(NormalizeEarth,h)

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

const NEARBY_LEFT_INDEX=+1
const NEARBY_RIGHT_INDEX=-1
const NEARBY_TOP_INDEX=-1
const NEARBY_BOTTOM_INDEX=+1
function _intersection_type(lr::LR,position,wedge_index,h_levels,refractive_map) where LR<:RadiusIntersection

  # left is to the left +1 and right is to the right -1
  position_interface_index2= lr==LeftIntersection() ? wedge_index[2]+NEARBY_LEFT_INDEX : wedge_index[2]+NEARBY_RIGHT_INDEX

  index=(wedge_index[1],position_interface_index2)
  @debug "position: $position, index: $index"
  h=convert(LLA2D{NormalizeEarth},ECEF2D{NormalizeEarth}(position...)) |> x-> x.h
  n₁=refractive_map[index...]
  # if the last radius was ABOVE the maximum level of stratification,
  # then the ray is in free space and you stop the computation
  # TO DO maybe find a better breaking condition
  # like checking the altitude of the ray outside instead of basing it on
  # the index
  @debug "h: $h, h_level_max: $(h_levels[1])"
  index=ifelse(h<=h_levels[1],index,(0,index[2]))
  (n₁,index,h)
end




function _intersection_type(tb::TB,position,wedge_index,h_levels,refractive_map) where TB<:LevelIntersection
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
  previous_index=input_index
  ###########################

  t_radius_l=advance(LeftIntersection(),ray,line_radii[input_index[2]+LEFT_RADIUS_INDEX])
  t_radius_r=advance(RightIntersection(),ray,line_radii[input_index[2]+RIGHT_RADIUS_INDEX])

  @debug "BottomIntersection"
  t_level_b=advance(BottomIntersection(),ray,scale_levels[input_index[1]+BOTTOM_LEVEL_INDEX])
  t_level_t=advance(TopIntersection(),ray,scale_levels[input_index[1]+TOP_LEVEL_INDEX])

  t_radius= min(t_radius_l,t_radius_r)
  leftright= t_radius_l<t_radius_r ? LeftIntersection() : RightIntersection() # left or right
  t_level= min(t_level_b,t_level_t)

  @debug "t_radius: $t_radius,t_left:$t_radius_l,t_right: $t_radius_r"
  @debug "t_level: $t_level,t_bottom:$t_level_b,t_top: $t_level_t"
  bottomtop= t_level_b<t_level_t ? BottomIntersection() : TopIntersection() # bottom or top


  debug_intesection=getDebugIntersection()
  @debug "DEBUG_INTERSECTION: $debug_intesection"
  if debug_intesection==0
    #######################
    #REAL
    #######################
    radius_or_level=t_radius<t_level ? leftright : bottomtop
    t_real= min(t_radius,t_level)
  elseif debug_intesection==1
    #######################
    # DEBUG ONLY RADII
    #######################
    @debug "DEBUG ONLY RADII"
    radius_or_level= leftright #t_radius<t_level ? leftright : bottomtop
    t_real= t_radius #min(t_radius,t_level)
  elseif debug_intesection==2
    #######################
    # DEBUG ONLY LEVELS
    #######################
    @debug "DEBUG ONLY LEVELS"
    radius_or_level= bottomtop
    t_real= t_level
    #######################
  else
    throw(ArgumentError("Invalid debug intersection, use setDebugIntersection(x) x=0,1,2 \n 0: Real, 1: Debug only radii, 2: Debug only levels"))
  end
  @debug "radius_or_level: $radius_or_level"
  @debug "t_real: $t_real"
  @debug "t_level: $t_level"
  @debug "t_radius: $t_radius"


  (n₁,index,h)=_intersection_type(radius_or_level,ray(t_real),input_index,h_levels,refractive_map)


  n₀=refractive_map[input_index...]
  @debug "t: $t_real, n₀: $n₀, n₁: $n₁, h: $h, whatMatch: $radius_or_level, index: $index"

  ray_out,isReflected = bend(radius_or_level,ray,t_real;n₀=n₀,n₁=n₁,h=h)
  index= isReflected ? input_index : index
  # condition to stop ray tracing
  target = 1<index[1]<(length(h_levels)-1)
  return (ray_out,target,h,index,radius_or_level)
end
