######################################################
using DataAPI
using DataAPI: colmetadatasupport, colmetadata, colmetadatakeys,colmetadata!
using StaticArrays
using Lazy
using GeoUtils
using StructArrays

abstract type AbstractMetaArray{T,N} <: AbstractArray{T,N} end
const AbstractMetaVector{T}=AbstractMetaArray{T,1}
const AbstractMetaMatrix{T}=AbstractMetaArray{T,2}



@forward AbstractMetaArray.data Base.getindex, Base.setindex!, Base.size, Base.eltype, Base.parent
@forward AbstractMetaArray.data  Base.size, Base.getindex, Base.setindex!, Base.similar, Base.axes
@forward AbstractMetaArray.data Base.iterate,Base.eltype, Base.parent
DataAPI.colmetadatasupport(::M) where {M<:AbstractMetaArray{T,N}} where {T,N} = (read=true, write=true)

function metadata(x::M,key::AbstractString; style::Bool=false) where M<:AbstractMetaArray
  if haskey(x.meta,key)
    return x.meta[key]
  end
end
metadata(x::M; style::Bool=false) where M<:AbstractMetaArray = x.meta



# Define a custom `similar` method for MetaArray
function Base.similar(x::A, ::Type{S}, dims::Dims{N}) where A<:AbstractMetaArray{T,N} where {T,S, N}
  # Create a new MetaArray with the same metadata and the specified element type and dimensions
  A(Array{T}(undef, dims), deepcopy(x.meta))
end

# Overload `similar` for cases without dimensions
function Base.similar(x::A, ::Type{S}) where A<:AbstractMetaArray{T,N} where {T, S, N}
  similar(x, S, size(x))
end


metacomponent(m::AbstractMetaArray,key) =  getfield(m,key)

metacomponent(x::StaticArray, i::Int)= getindex(x,i)
function metacomponent(m::AbstractMetaArray{<:Union{SVector,MVector}},key::Symbol)
  i = key == :x ? 1 :
      key == :y ? 2 :
      key == :z ? 3 :
      key == :w ? 4 :
      throw(ArgumentError("Invalid key for metacomponent: $key"))
      metacomponent(m,i)
end


metacomponents(m::AbstractMetaArray)=parent(m)
function Base.getproperty(m::AbstractMetaArray, key::Symbol)
  if key === :data || key === :meta
      # Directly access the fields of MetaArray
      return getfield(m, key)
  elseif isa(getfield(m, :data), StructArray)
      # Delegate to StructArray's component access
      return StructArrays.component(getfield(m, :data), key)
  else
      # Handle other cases (e.g., metadata or custom components)
      return metacomponent(m, key)
  end
end
Base.getproperty(m::AbstractMetaArray,key::Int) = metacomponent(m,key)
Base.propertynames(m::AbstractMetaArray,bool::Bool=true)= propertynames(metacomponents(m))
metadatakeys(x::M) where M<:AbstractMetaArray = keys(x.meta)

Base.convert(::Type{A},v::AbstractArray) where A<:AbstractMetaArray = A(v)
Base.convert(::Type{A},v::A) where A<:AbstractMetaArray = v

@forward AbstractMetaArray.data Base.resize!

Base.reshape(m::A, d::Vararg{Int64}) where A<:AbstractMetaArray{T} where T = A(reshape(parent(m),d), m.meta)

struct MetaArray{T,N}<:AbstractMetaArray{T,N}
  data
  meta::Dict

  function MetaArray{T,N}(data::A, meta::Dict=Dict()) where {T,N,A<:AbstractArray{T,N}}
    new{T,N}(data,meta)
  end

  function MetaArray{T}(data::A,meta::Dict=Dict()) where {T,A<:AbstractArray{T,N}} where N
    new{T,N}(data,meta)
  end
end
MetaArray(data::AbstractArray{T,N}, meta::Dict=Dict()) where {T,N} = MetaArray{T}(data,meta)

a=MetaArray(rand(3,3),Dict(:x=>1,:y=>2,:z=>3))
b=MetaArray(StructArray{Orbit{Float64}}(undef,3,3),Dict(:x=>1,:y=>2,:z=>3))


struct GeoMetaArray{T,N} <: AbstractMetaArray{T,N}
  data
  meta::Dict
  function GeoMetaArray{T,N}(data::A, meta::Dict=Dict()) where {T,N,A<:AbstractArray{T,N}}
    _meta= fieldnames(eltype(data))
    for key in _meta
      if !haskey(meta,key)
        meta[key]=unit()
      end
    end

    new{T,N}(data,meta)
  end

  function GeoMetaArray{T}(data::A,meta::Dict=Dict()) where {T,A<:AbstractArray{T,N}} where N
    GeoMetaArray{T,N}(data,meta)
  end
end
