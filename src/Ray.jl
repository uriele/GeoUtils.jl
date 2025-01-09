
const SVec3{T} = SVector{3,T}

const SVec2{T} = SVector{2,T}

const SVec2(x::T,y::T=T(0)) where T<:Number=SVec2{T}(x,y)
const SVec3(x::T,y::T=T(0),z::T=T(0)) where T<:Number=SVec3{T}(x,y,z)

const Vec3F64 = SVec3{Float64}

const Vec2F64 = SVec2{Float64}

const Vec3F32 = SVec3{Float32}

const Vec2F32 = SVec2{Float32}

SDiagonal(v::SVec2{T}) where T = Diagonal(v)
SDiagonal(v::SVec3{T}) where T = Diagonal(v)

const MVec2{T} = MVector{2,T}
const MVec3{T} = MVector{3,T}

const MVec2(x::T,y::T=T(0)) where T<:Number=MVec2{T}(x,y)
const MVec3(x::T,y::T=T(0),z::T=T(0)) where T<:Number=MVec3{T}(x,y,z)

const MVec3F64 = MVec3{Float64}
const MVec2F64 = MVec2{Float64}
const MVec3F32 = MVec3{Float32}
const MVec2F32 = MVec2{Float32}

"""
  AbstractRay{T<:IEEEFloat,N}

  An abstract ray with origin and direction. The origin and direction can either be mutable or static vectors.

  # Interface

  - `get_origin(r)`: return the origin
  - `get_direction(r)`: return the direction

  # Methods

  - `r(t::T)::SVN`: Returns the point at distance `t` from the origin in the direction of the ray.
"""
abstract type AbstractRay{T<:IEEEFloat,N} end
(r::AbstractRay{T})(t) where T= r.origin+T(t)*_normalize(r.direction)

"""
  AbstractRay2D{T}

  Alias of AbstractRay{T,2}.
"""
const AbstractRay2D{T}=AbstractRay{T,2}
"""
  AbstractRay3D{T}

  Alias of AbstractRay{T,3}
"""
const AbstractRay3D{T}=AbstractRay{T,3}



"""
  Ray2D{T,SV2<:StaticVector{2,T}}(origin::SV2,direction::SV2)

  A 2D ray with origin and direction. The origin and direction can either be mutable or static vectors.
"""
struct Ray2D{T,SV2<:StaticVector{2,T}}<:AbstractRay2D{T}
  origin::SV2
  direction::SV2

  @inline function Ray2D(origin::SV,direction::SV) where {SV<:StaticVector{2,T}} where T
    _hypothenuse(direction)==0 && throw(ArgumentError("Direction cannot be zero"))
    return new{T,SV}(origin,_normalize(direction))
  end
end

"""
  Ray3D{T,SV3<:StaticVector{3,T}}(origin::SV3,direction::SV3)

  A 3D ray with origin and direction. The origin and direction can either be mutable or static vectors.
"""
struct Ray3D{T,SV3<:StaticVector{3,T}}<:AbstractRay3D{T}
  origin::SV3
  direction::SV3

  @inline function Ray3D(origin::SV,direction::SV) where {SV<:StaticVector{3,T}} where T
    _hypothenuse(direction)==0 && throw(ArgumentError("Direction cannot be zero"))
    return new{T,SV}(origin,_normalize(direction))
  end
end


const SRay2D{T} = Ray2D{T,SVec2{T}}
const SRay3D{T} = Ray3D{T,SVec3{T}}
const MRay2D{T} = Ray2D{T,MVec2{T}}
const MRay3D{T} = Ray3D{T,MVec3{T}}

"""
  SRay2D(x::T,y::T,dx::T,dy::T)

  Create a 2D ray with origin (x,y) and direction (dx,dy) with static vectors.
"""
SRay2D(x::T,y::T,dx::T,dy::T) where T = Ray2D(SVec2(x,y),SVec2(dx,dy))
"""
  MRay2D(x::T,y::T,dx::T,dy::T)

  Create a 2D ray with origin (x,y) and direction (dx,dy) with mutable vectors.
"""
MRay2D(x::T,y::T,dx::T,dy::T) where T = Ray2D(MVec2(x,y),MVec2(dx,dy))
"""
  SRay3D(x::T,y::T,z::T,dx::T,dy::T,dz::T)

  Create a 3D ray with origin (x,y,z) and direction (dx,dy,dz) with static vectors.
"""
SRay3D(x::T,y::T,z::T,dx::T,dy::T,dz::T) where T = Ray3D(SVec3(x,y,z),SVec3(dx,dy,dz))
"""
  MRay3D(x::T,y::T,z::T,dx::T,dy::T,dz::T)

  Create a 3D ray with origin (x,y,z) and direction (dx,dy,dz) with mutable vectors.
"""
MRay3D(x::T,y::T,z::T,dx::T,dy::T,dz::T) where T = Ray3D(MVec3(x,y,z),MVec3(dx,dy,dz))

# Interfaces
get_direction(r::R) where {R<:AbstractRay}=r.direction
get_origin(r::R) where {R<:AbstractRay}=r.origin

@inline _dot(v::SV2) where {SV2<:StaticVector{2,T}} where T = dot(v,v)::T
@inline _dot(v1::SV21,v2::SV22) where {SV21<:StaticVector{2,T},SV22<:StaticVector{2,T}} where T = dot(v1,v2)::T
@inline _dot(v::SV3) where {SV3<:StaticVector{3,T}} where T = dot(v,v)::T
@inline _dot(v1::SV31,v2::SV32) where {SV31<:StaticVector{3,T},SV32<:StaticVector{3,T}} where T = dot(v1,v2)::T
@inline _hypothenuse(v::SV2) where {SV2<:StaticVector{2,T}} where T = hypot(v.x,v.y)
@inline _normalize(v::SV2)  where {SV2<:StaticVector{2,T}} where T = v/_hypothenuse(v)
@inline _hypothenuse(v::SV3) where {SV3<:StaticVector{3,T}} where T = norm(v)
@inline _normalize(v::SV3)  where {SV3<:StaticVector{3,T}} where T = v/_hypothenuse(v)
