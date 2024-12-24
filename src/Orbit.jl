
const ATOL64 = ScopedValue(1.0e-10)
const ATOL32 = ScopedValue(1.0f-5)

"""
    atol(T)
    atol(x::T)

Absolute tolerance used in algorithms for approximate
comparison with numbers of type `T`. It is used in the
source code in calls to the [`isapprox`](@ref) function:

```julia
isapprox(a::T, b::T, atol=atol(T))
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

  @inline function Ellipsoid(a::T,b::T) where T
    radius=Vec2{T}(a,b)
    unscale=LinearMap(SDiagonal(radius))
    scale=inv(unscale)
    return new{T}(radius,scale,unscale)
  end
end
Ellipsoid(a::T) where T = Ellipsoid(a,a)

struct Radius{T<:IEEEFloat}
  a::T
  b::T
  c::T
end
Radius(r::Ray2D{T}) where T = Radius{T}(_coefficients(r)...)


@inline _coefficients(r::Ray2D{T}) where T = _coefficients(r.origin...,r.direction...)
@inline function _coefficients(x0::T,y0::T,dx::T,dy::T) where T
  a = dy
  b = dx
  c = -dx * y0 - dy * x0
  return (a,b,c)

end


#r=Ray2D(Vec2(-2.,343.),Vec2(-4.,2.43))

direction(r::Ray2D)=r.direction
origin(r::Ray2D)=r.origin



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
  # -1+x^2/a^2+y^2/b^2=0
  #
  # C= x0^2/a^2+y0^2/b^2-1
  # B= x0*dx/a^2+y0*dy/b^2
  # A= dx^2/a^2+dy^2/b^2
  #t^2+(mu-B/2)(-mu-B/2)*t+(mu^2-C)=0
  #t^2-(B/2)^2*t+(mu^2-C)=0
  C=origin.x*origin.x+origin.y*origin.y-1
  halfB=-(origin.x*direction.x+origin.y*direction.y)
  A= max(hypot(direction.x,direction.y),1e-10)
  C/=A
  halfB/=A

  b=halfB*halfB-C

  if b<=0
    return T(Inf)
  end

  b_sqrt =sqrt(b)
  return filter(x-> x>=0, [b_sqrt,-b_sqrt].+halfB) |> x-> reduce(min,x;init=Inf)

end

"""
    distance_from_unit_ellipse(origin::Vec2{T},direction::Vec2{T})::T where T
    distance_from_unit_elllipse(ray::Ray2D{T})::T where T

Calculate the distance from the unit circle to a ray or a point
"""
function distance_from_ellipse(origin::Vec2{T},direction::Vec2{T};b::T=1.0,h::T=0.0)::T where T
  # it is supposed to be divided by the square of the direction, but assuming normalized
  #
  # -1+x^2/a^2+y^2/b^2=0
  #
  # C= x0^2/a^2+y0^2/b^2-1
  # B= x0*dx/a^2+y0*dy/b^2
  # A= dx^2/a^2+dy^2/b^2
  #t^2+(mu-B/2)(-mu-B/2)*t+(mu^2-C)=0
  #t^2-(B/2)^2*t+(mu^2-C)=0
  invb²=1/(b+h) |> x-> x*x
  inva²=1/(1+h) |> x-> x*x

  scale=[inva²; invb²]

  origin=scale.*origin
  direction=scale.*direction

  C=origin.x*origin.x+origin.y*origin.y-1
  halfB=-(origin.x*direction.x+origin.y*direction.y)
  A= max(hypot(direction.x,direction.y),1e-10)
  C/=A
  halfB/=A

  b=halfB*halfB-C

  if b<=0
    return Inf
  end

  b_sqrt =sqrt(b)
  return filter(x-> x>=0, [b_sqrt,-b_sqrt].+halfB) |> x-> reduce(min,x;init=Inf)

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
  distance_from_radii(origin::Vec2{T},direction::Vec2{T},a::T,b::T,c::T)::T where T
  distance_from_radii(ray::Ray2D{T},r::Radius{T})::T where T

Calculate the distance from a ray to one of Earth Radii or from a line with coefficients a,b,c
"""
function distance_from_radii(origin::Vec2{T},direction::Vec2{T},a::T,b::T,c::T) where T
  # a*x+b*y+c=0
  # (a*x0+b*y0+c)+t*(a*dx+b*dy)=0
  # t=-(a*x0+b*y0+c)/(a*dx+b*dy)
  t= -((a*origin.x+b*origin.y+c)/(a*direction.x+b*direction.y))
  return t>=0 ? t : T(Inf)
end
distance_from_radii(ray::Ray2D{T},r::Radius{T}) where T = distance_from_radii(ray.origin,_normalize(ray.direction),r.a,r.b,r.c)


"""
    distance_from_segment(A::Vec2{T},B::Vec2{T},C::Vec2{T},D::Vec2{T})::T where T
    distance_from_segment(ray::Ray2D{T},C::Vec2{T},D::Vec2{T})::T where T

Calculate the distance from a segment AB to a segment CD
"""
function distance_from_segment(A::Vec2{T},B::Vec2{T},C::Vec2{T},D::Vec2{T})::T where T
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
  ((r>=2) && (M\y))[1] |> x-> (0<x<=1) ? x : T(Inf)
end

distance_from_segment(ray::Ray2D{T},C::Vec2{T},D::Vec2{T}) where T = distance_from_segment(ray(0),ray(1),C,D)



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
  orig=scale(r.origin)
  direc=scale(r.direction)
  hyp_direct= hypot(direc...)
  t=distance_from_unit_circle(orig,direc./hyp_direct) |> x-> x/hyp_direct
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

#=
@inline function bend(r::Ray2D,t,ii)
  neworigin=r(t)
  h=convert(LLA2D{NormalizeEarth},ECEF2D{NormalizeEarth}(neworigin[1],neworigin[2])) |> x-> x.h
  N=_normal_vector(neworigin...,h,e²) |> x-> x/hypot(x...)
  a=-N⋅r.direction
  b=ii.n²*(1-a^2)
  newdir= if b≤1 # refract or reflect
    ii.n*r.direction+(ii.n*a-sqrt(1-b))*N
  else
    r.direction-2*a*N
  end
  return Ray2D(neworigin,newdir)
end
=#


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


const SHIFTORIGIN=1e-8

@inline function _bend_initialize(ray::Ray2D,t::T, n₀::T,n₁::T) where T
  neworigin=ray(t)
  direction=ray.direction
  n₀₁=n₀/n₁
  return (neworigin,direction,n₀₁,n₀₁*n₀/n₁)
end

@inline function _bend_ellipse(ray::Ray2D,t::T, n₀::T,n₁::T, h::T=T(0), outward=1.0) where T
  e²=eccentricity²(ellipsoid(NormalizeEarth))
  ## shift the origin to avoidnumerical issues

  (neworigin,direction,n₀₁,n₀₁²)=_bend_initialize(ray,t,n₀,n₁)
  N=_normal_vector(neworigin...,h,e²) |> x-> x/hypot(x...).*outward
  ray_out,isReflected=_bend_common(neworigin,direction,N,n₀₁,n₀₁²)


  isRising= acosd(-N⋅direction)>90
  @debug "isRising: $isRising, isReflected: $isReflected"
  return (ray_out,isReflected,isRising)
end


@inline function _bend_radii(ray::Ray2D,t::T, n₀::T,n₁::T, h::T=T(0), outward=1.0) where T
  e²=eccentricity²(ellipsoid(NormalizeEarth))
  (neworigin,direction,n₀₁,n₀₁²)=_bend_initialize(ray,t,n₀,n₁)
  N=_tangent_vector(neworigin...,h,e²) |> x-> x/hypot(x...).*outward
  ray_out,isReflected=_bend_common(neworigin,direction,N,n₀₁,n₀₁²)
  isRising= acosd(-N⋅direction) #|> x-> ifelse(abs(x)==90,x,rem(x,180,RoundNearest) )<0
  @debug "θ: $isRising"
  isRising= isRising>90
  @debug "isRising: $isRising, isReflected: $isReflected"
  return (ray_out,isReflected,isRising)
end


#TopLevelIntersection -1
@inline bend(::TopIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)    = _bend_ellipse(ray,t,n₀,n₁,h,-1.0) |> x-> (x[1],x[2])
#RightRadiusIntersection -1
@inline bend(::RightIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)  = _bend_radii(ray,t,n₀,n₁,0.0,1.0)  |> x-> (x[1],x[2])


#BottomLevelIntersection 1
@inline bend(::BottomIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0) = _bend_ellipse(ray,t,n₀,n₁,h,1.0)  |> x-> (x[1],x[2])
#LeftRadiusIntersection 1
@inline bend(::LeftIntersection,ray::Ray2D,t;n₀=1.0,n₁=1.0,n₂=1.0,h=0.0)   = _bend_radii(ray,t,n₀,n₁,0.0,-1.0)  |> x-> (x[1],x[2])

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
@inline function _bend_common(origin::Vec2{T},direction::Vec2{T},
  N::Vec2{T},n₀₁::T,n₀₁²::T) where T

  # Base condition
  isReflected=false

  a=-N⋅direction

  b=n₀₁²*(1-a^2)
  @debug "b: $b, θ: $(acosd(a))"
  @debug "θ: $(acosd(a)),sinθ: $(a), n₀₁²: $(n₀₁²)"
  newdir= if b≤1 # refract or reflect
    n₀₁*direction+(n₀₁*a-sqrt(1-b))*N
  else
    direction-2*a*N
    isReflected=true
    @debug "is reflected"
    #readline()
  end
  # SHIFT ORIGIN TO AVOID NUMERICAL ISSUES
  return Ray2D(origin+newdir*SHIFTORIGIN,newdir),isReflected
end

@inline _normal_vector(x0,y0,h,e²)=Vec2(x0/(1+h)^2,y0/(sqrt(1-e²)+h)^2) |> x->x/hypot(x...)
@inline _tangent_vector(args...)=_normal_vector(args...) |> x-> Vec2(-x.y,x.x)
# Create Rays from the orbit


function create_rays(datum::Datum,orb::Orbit) where Datum
  e²=eccentricity²(ellipsoid(datum))
  (w,z)=(orb.w,orb.z)
  (nx,ny)=_normal_vector(w,z,0,e²)
  #angle computed with respect to the tangent plane
  (tx,ty)=(-ny,-nx)
  angle= orb.ang
  @inline _rotation_matrix(θ)= SMatrix{2,2}([cosd(θ) sind(θ);-sind(θ) cosd(θ)] )
  return Ray2D(Vec2(w,z),_rotation_matrix(-angle)*Vec2(tx,ty))
end

create_rays(orb::Orbit)=create_rays(NormalizeEarth,orb)


function getIntersectionObjects(lla,ecef,b::T) where T
  (nlevels,nradii)=size(lla)
  h_levels = Vector{Float64}(undef,nlevels)
  scale_levels = Vector{Ellipsoid{Float64}}(undef,nlevels)
  θ_radii  = Vector{Float64}(undef,nradii)
  line_radii = Vector{Radius{Float64}}(undef,nradii)
  invb²=1/(b*b)
  [
    let
      h=lla[i,1][1]
      h_levels[i,1]=h
      scale_levels[i]=Ellipsoid(1.0+h,b+h)
    end
    for i in axes(lla,1)
  ]
  [
    let
      origin=Vec2(ecef[end,j][1],ecef[end,j][2])
      direction=Vec2(origin.x,origin.y*invb²)
      θ_radii[j]=lla[end,j][2]
      line_radii[j]=Radius(Ray2D(origin,direction))
    end
    for j in axes(ecef,2)
  ]
  return h_levels,StructArray(scale_levels).scale,SemiCircularArray(θ_radii),SemiCircularArray(line_radii)
end

@inline function new_intersection(ray::Ray2D, # ray
  wedge_index=wedge_index; # wedge index
  refractive_index_map=refractive, # map of refractive index
  h_levels=h_levels, # levels for i
  θ_radii=θ_radii, # radii angles for j
  scale_levels=scale_levels, # scaling mapping to unit circle
  line_radii=line_radii, # line radii
  tangent_quote=tangent_quote, # tangent quote
  register=register,
  #for debug
  previous_intersection=nothing) # register

  nlevels=size(refractive_index_map,1) #number of levels to check for the bottom intersection

  n₀= refractive_index_map[wedge_index...]
  @debug "n₀: $n₀"
  @debug "wedge_index: $(wedge_index)"
  @debug "line_radii: $(line_radii[wedge_index[2]])"
  @debug "scale_levels: $(scale_levels[wedge_index[1]+1])"
  # find the next intersection
  # for simplicity, I only need to check only Left and Bottom
  @debug "LeftIntersection"
#############################
previous_index=wedge_index
############################

  t_radius_l=advance(LeftIntersection(),ray,line_radii[wedge_index[2]])
  t_radius_r=advance(RightIntersection(),ray,line_radii[wedge_index[2]+1])

  t_radius= min(t_radius_l,t_radius_r)
  leftright= t_radius_r>t_radius_l ? 1 : -1 # left or right
  @debug "BottomIntersection"
  t_level_b=advance(BottomIntersection(),ray,scale_levels[wedge_index[1]+1])
  t_level_t=advance(TopIntersection(),ray,scale_levels[wedge_index[1]])
  t_level= min(t_level_b,t_level_t)
  bottomtop= t_level_b<t_level_t ? 1 : -1 # bottom or top

  if isinf(t_level)
    @debug "level failed"
    @debug  ray
    @debug wedge_index[1]+1,scale_levels[wedge_index[1]+1]
    @debug wedge_index[1],scale_levels[wedge_index[1]]
    readline()
  end

  @debug "t_level: $t_level,$t_level_b,$t_level_t, t_radius: $t_radius,$t_radius_l,$t_radius_r"
  @debug "Δt: $(t_radius-t_level)"

  Δt=(t_radius-t_level)

  @inline function _radius_intersection(ray::Ray2D,t,wedge_index,nlevels,leftright::Int=1)
    whatMatch= if leftright==1
      LeftIntersection()
    else
      RightIntersection()
    end
    index=(wedge_index[1],wedge_index[2]+leftright)
    h=convert(LLA2D{NormalizeEarth},ECEF2D{NormalizeEarth}(ray(t)...)) |> x-> x.h
    @debug "$whatMatch input: $wedge_index, output: $index, h: $h"
    @debug whatMatch
    (whatMatch,
    refractive_index_map[index...],
    index,
    h
    )
  end

  @inline function _ellipse_intersection(ray::Ray2D,t,wedge_index,nlevels,bottomtop::Int=1)
    whatMatch= if bottomtop==1
      BottomIntersection()
    else
      TopIntersection()
    end

    hindex=wedge_index[1]+bottomtop
    index=(hindex,wedge_index[2])

    @debug "bottomtop: $bottomtop"
    @debug "whatMatch: $whatMatch"
    @debug "index: $index"
    @debug "hindex: $hindex"
    @debug "nlevels: $nlevels"
    if hindex<1  # free space
      n₁=1.0
    elseif hindex>=nlevels  # ground full reflection
      n₁=Inf
    else            # atmosphere
      n₁=refractive_index_map[index...]
    end
    h= if whatMatch isa BottomIntersection
      h_levels[wedge_index[1]+1]             # level of the bottom intersection
    else # whatMatch isa TopIntersection
      h_levels[wedge_index[1]]               # level of the top intersection
    end
    @debug "$whatMatch input: $wedge_index, output: $index, h: $h"
    (whatMatch,n₁,index,h)
  end


  (t,(whatMatch,n₁,index,h))= if Δt<0
    (t_radius,_radius_intersection(ray,t_radius,wedge_index,nlevels,leftright))
  else
    (t_level,_ellipse_intersection(ray,t_level,wedge_index,nlevels,bottomtop))
  end

  n₀=refractive_index_map[wedge_index...]
  @debug "t: $t, n₀: $n₀, n₁: $n₁, h: $h, whatMatch: $whatMatch, index: $index"

  b,isReflected = bend(whatMatch,ray,t;n₀=n₀,n₁=n₁,h=h)
  index= isReflected ? wedge_index : index
  # condition to stop ray tracing
  target = 1<index[1]<nlevels
  return (b,target,h,index,whatMatch)
end
