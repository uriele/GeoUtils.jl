
"""
  SemiCircularArray{T,N,A} <: AbstractArray{T,N}

  'N'-dimensional array with circular indexing.
  If N==1, the array is indexed as a circular array.
  If N>1, all the dimension of the array except the first one are indexed as a circular array.
"""
struct SemiCircularArray{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
  data::A
  SemiCircularArray{T,N}(data::A) where A<: AbstractArray{T,N} where {T,N} = new{T,N,A}(data)
  SemiCircularArray{T,N,A}(data::A) where A<: AbstractArray{T,N} where {T,N} = new{T,N,A}(data)
end

#const SemiCircularArray{T,1}(data::A) where A<:AbstractVector{T} where T = CircularVector{T}(data)


"""
  SemiCircularVector{T} <: AbstractVector{T}

  1-dimensional array with circular indexing.
"""
const SemiCircularVector{T} = SemiCircularArray{T,1}

"""
  SemiCircularArray{T} <: AbstractMatrix{T}

  2-dimensional array with circular indexing in the second dimension.
"""
const SemiCircularMatrix{T} = SemiCircularArray{T,2}

"""
  SemiCircularArray(data::A) where A<: AbstractArray{T,N}

Create a SemiCircularArray from an AbstractArray.
"""
SemiCircularArray{T}(data::AbstractArray{T,N}) where {T,N} = SemiCircularArray{T,N}(data)
SemiCircularArray(data::AbstractArray{T,N}) where {T,N} = SemiCircularArray{T}(data)

"""
  SemiCircularArray(def::T,size) where T

Create a SemiCircularArray with a default value and a size.
"""
SemiCircularArray(def::T,size) where T = SemiCircularArray{T}(fill(def,size))

"""
  SemiCircularVector(data::A) where A<: AbstractVector{T}

Create a SemiCircularVector from an AbstractVector.
"""
SemiCircularVector(data::AbstractVector{T}) where T = SemiCircularVector{T}(data)

"""
  SemiCircularVector(def::T,size) where T

Create a SemiCircularVector with a default value and a size.
"""
SemiCircularVector(def::T,size) where T = SemiCircularVector{T}(fill(def,size))

"""
  SemiCircularMatrix(data::A) where A<: AbstractMatrix{T}

Create a SemiCircularArray from an AbstractMatrix.
"""
SemiCircularMatrix(data::AbstractMatrix{T}) where T = SemiCircularMatrix{T}(data)

"""
  SemiCircularMatrix(def::T,size) where T

Create a SemiCircularMatrix with a default value and a size.
"""
SemiCircularMatrix(def::T,size) where T = SemiCircularArray{T}(fill(def,size))


Base.IndexStyle(::Type{SemiCircularArray{T,N,A}}) where {T,N,A} = IndexCartesian()
Base.IndexStyle(::Type{SemiCircularVector}) = IndexLinear()


@inline _dimension_text(arr::AbstractArray)= join(string.(size(arr)),"×")


@inline Base.getindex(arr::SemiCircularArray,i::Int) = @inbounds getindex(arr.data,mod(i,eachindex(IndexLinear(),arr.data)))

# handle the ambiguities of getindex when N==1
@inline Base.getindex(arr::SemiCircularVector,i::Int,j::Vararg{Int,N}) where N = @inbounds getindex(arr.data,mod(i,eachindex(IndexLinear(),arr.data)))

@inline function Base.getindex(arr::SemiCircularArray{T,N},i::Int,j::Vararg{Int,M}) where {T,N,M}
  try
    (M==N-1) || error()
    return getindex(arr.data,
        i,mod.(j,axes(arr.data)[2:end])...)
  catch
    _j=[j...]
    rethrow(BoundsError(arr,pushfirst!(_j,i)))
  end
end

@inline Base.setindex!(arr::SemiCircularArray,v,i::Int) = @inbounds setindex!(arr.data,v,mod(i,eachindex(IndexLinear(),arr.data)))

# handle the ambiguities of set index when N==1
@inline Base.setindex!(arr::SemiCircularVector,v,i::Int,j::Vararg{Int,N}) where N = @inbounds setindex!(arr.data,v,mod(i,eachindex(IndexLinear(),arr.data)))

@inline function Base.setindex!(arr::SemiCircularArray{T,N},v,i::Int,j::Vararg{Int,M}) where {T,N,M}
  try
    (M==N-1) || error()
    return setindex!(arr.data,v,
          i,mod.(j,axes(arr.data)[2:end])...)
  catch
    _j=[j...]
    rethrow(BoundsError(arr, pushfirst!(_j,i)))
  end
end


@inline Base.size(arr::SemiCircularArray) = size(arr.data)
@inline Base.axes(arr::SemiCircularArray) = axes(arr.data)
@inline Base.parent(arr::SemiCircularArray) = arr.data

@inline Base.iterate(arr::SemiCircularArray, i...) = iterate(parent(arr), i...)

@inline Base.in(x, arr::SemiCircularArray) = in(x, parent(arr))
@inline Base.copy(arr::SemiCircularArray) = SemiCircularArray(copy(parent(arr)))

@inline function Base.checkbounds(arr::SemiCircularArray, I...)
  J=Base.to_indices(arr,I)
  length(J)==1 || length(J)>= ndims(arr) || throw(BoundsError(arr,I))
  return nothing
end

@inline _similar(arr::SemiCircularArray, ::Type{T}, dims) where T = SemiCircularArray(similar(parent(arr),T,dims))
@inline Base.similar(arr::SemiCircularArray, ::Type{T},  dims::Tuple{Base.DimOrInd, Vararg{Base.DimOrInd}}) where T = _similar(arr,T,dims)
# ambiguities resolved with base
@inline Base.similar(arr::SemiCircularArray, ::Type{T}, dims::Dims) where T = _similar(arr, T, dims)
@inline Base.similar(arr::SemiCircularArray, ::Type{T}, dims::Tuple{Integer, Vararg{Integer}}) where T = _similar(arr, T, dims)
@inline Base.similar(arr::SemiCircularArray, ::Type{T}, dims::Tuple{Union{Integer, Base.OneTo}, Vararg{Union{Integer, Base.OneTo}}}) where T = _similar(arr, T, dims)


@inline _similar(::Type{SemiCircularArray{T,N,A}}, dims) where {T,N,A} = SemiCircularArray{T,N}(similar(A,dims))
@inline Base.similar(SCA::Type{SemiCircularArray{T,N,A}}, dims) where {T,N,A} = _similar(SCA,dims)
# Ambiguity resolution with Base
@inline Base.similar(SCA::Type{SemiCircularArray{T,N,A}}, dims::Dims) where {T,N,A} = _similar(SCA, dims)
@inline Base.similar(SCA::Type{SemiCircularArray{T,N,A}}, dims::Tuple{Union{Integer, Base.OneTo}, Vararg{Union{Integer, Base.OneTo}}}) where {T,N,A} = _similar(SCA, dims)


@inline Broadcast.BroadcastStyle(::Type{SemiCircularArray{T,N,A}}) where {T,N,A} = Broadcast.ArrayStyle{SemiCircularArray{T,A}}()
@inline Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SemiCircularArray{T,N,A}}}, ::Type{ElType}) where {T,N,A,ElType} = SemiCircularArray(similar(convert(Broadcast.Broadcasted{typeof(Broadcast.BroadcastStyle(A))},bc),ElType))

@inline Base.dataids(arr::SemiCircularArray) = Base.dataids(arr.data)

function Base.showarg(io::IO, arr::SemiCircularArray, toplevel)
  nd=ndims(arr)
  stype="SemiCircularArray"
  if nd==1
    stype="SemiCircularVector"
  end
  if nd==2
    stype="SemiCircularMatrix"
  end
  print(io, stype,"(")
  Base.showarg(io, parent(arr), false)
  print(io, ')')
  # toplevel && print(io, " with eltype ", eltype(arr))
end


Base.empty(::SemiCircularVector{T}, ::Type{U}=T) where {T,U} = SemiCircularVector{U}(U[])
Base.empty!(arr::SemiCircularVector) = (empty!(parent(arr)); arr)
Base.push!(arr::SemiCircularVector, item...) = (push!(parent(arr), item...); arr)
Base.append!(arr::SemiCircularVector, items) = (append!(parent(arr), items); arr)
Base.resize!(arr::SemiCircularVector, n1::Integer) = (resize!(parent(arr), n1); arr)
Base.pop!(arr::SemiCircularVector) = pop!(parent(arr))
Base.sizehint!(arr::SemiCircularVector, sz::Integer) = (sizehint!(parent(arr), sz);arr)

function Base.deleteat!(arr::SemiCircularVector, i::Integer)
  deleteat!(arr.data, mod(i,eachindex(IndexLinear(),arr.data)))
  return arr
end

function Base.insert!(arr::SemiCircularVector, i::Integer, item)
  insert!(arr.data, mod(i,eachindex(IndexLinear(),arr.data)), item)
  return arr
end

#=
N=20
M=36

data=Matrix{Float64}(undef,N,M)

[ data[i,j]= for i in 1:N, j in 1:M]



OpticalInfo[1,1]=1.0
=#
