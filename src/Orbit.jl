
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
  # direction it is not necessary and save a division
  C=_dot(origin)-T(1)  # distance of the ray source from the origin
  halfB=-_dot(origin,direction) # projection of the origin on the direction
  μ²=halfB*halfB-C
  sqrt_llvm(μ²) |>             # use llvm intrinsic for sqrt for fast math,if negative it will return NaN
  x -> (halfB*ones(2)+[x,-x]) |>
  x -> reduce(min,filter(x->x>=0,x);init=T(Inf))  # filter takes advantage of the fact that NaN return false for >0 <0 or ==0
end
distance_from_unit_circle(ray::Ray2D{T}) where T = distance_from_unit_circle(ray.origin,_normalize(ray.direction))

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
  #@info "A=$A B=$B C=$C D=$D"
  M=[(B-A) (C-D)];
  y=(C-A);
  #@info "M=$M y=$y"
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

 function distance_from_radius(origin::Vec2{T},direction::Vec2{T},angle)::T where T
  sinθ= sin(angle)
  cosθ= cos(angle)
  v=[sinθ;-cosθ]
  N=dot(origin,v)
  D=-dot(direction,v)
  τ=atol(T)
  D=-dot(direction,v) |> x-> !(-τ<x<τ)*x
  reduce(min,filter(x->x>0,[N/D]);init=T(Inf))
end

function distance_from_radius_new(origin::Vec2{T},direction::Vec2{T},angle)::T where T
  sinθ= sin(angle)
  cosθ= cos(angle)
  v=[sinθ;-cosθ]
  N=dot(origin,v)
  τ=atol(T)
  D=-dot(direction,v) |> x-> !(-τ<x<τ)*x
  _distance=N/D
  #reduce(min,filter(x->x>=0,[N/D]);init=T(Inf))
  _distance > 0 ? _distance : T(Inf)
end
distance_from_radius(ray::Ray2D{T},angle) where T = distance_from_radius(ray.origin,ray.direction,angle)
distance_from_radius_new(ray::Ray2D{T},angle) where T = distance_from_radius(ray.origin,ray.direction,angle)
