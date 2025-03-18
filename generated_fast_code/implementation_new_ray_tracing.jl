# load packages
using Base: IEEEFloat
using Makie,WGLMakie
using CoordRefSystems # for information on earth

using MacroTools
using Unitful
#using Unrolled
using BenchmarkTools #to benchmark
include("../generated_fast_code/global_constants.jl")
include("../generated_fast_code/quoted_functions.jl")

# generate setter and getters functions for the minoraxis to test different Earth
"""
  set_minoraxis([b=DEFAULT_MINORAXIS])

Set the normalized minor_axis of Earth ,`b`, and the relative constants `b²` and `e²=1-b²`. If the
field is empty it sets to the value of the `WGS84` Earth model. The value has to be b ∈ (0,1].

The function set the value for `Float16`, `Float32`, `Float64`.

See: [`get_minoraxis`](@ref), [`get_minoraxis²`](@ref), [`get_e²`](@ref)

"""  function set_minoraxis end
@eval @generated function set_minoraxis(minor_axis=DEFAULT_MINORAXIS)::Nothing
  if minor_axis<:Real
    types=(:F16,:F32,:F64)
    refs=(:MINORAXIS,:MINORAXIS²,:E²)
    vals=(:minor_axis,:minor_axis_squared,:e_squared)
    exps=Expr[]
    push!(exps,:(minor_axis_squared=minor_axis^2))
    push!(exps,:(e_squared=1-minor_axis_squared))
    for typ in types
      for (ref,val) in zip(refs,vals)
        symbol=Symbol("$(ref)_$typ")
        push!(exps,:( $symbol[]=$val))
      end
    end

    quote
      @assert 0<minor_axis<1 "The minor axis must be ∈ (0,1]"

      $(Expr(:block,exps...))
      return nothing
    end
  else
    return :(error("Unsupported type $minor_axis"))
  end
end


for (infoconstant,constantname) in zip((:minoraxis,:minoraxis²,:e²),(:b,:b²,:e²))
  lower_const_name=String(infoconstant)
  upper_const_name=uppercase(lower_const_name)
  upper_const_symbol_F16=Symbol("$(upper_const_name)_F16")
  upper_const_symbol_F32=Symbol("$(upper_const_name)_F32")
  upper_const_symbol_F64=Symbol("$(upper_const_name)_F64")
  fn=Symbol("get_$lower_const_name")


  docstring="""
    $(fn)(::Type{T}=Float64) where T

  Generated functions to retrieve the $(infoconstant) constants `$(constantname)` based on the type `T`.

  # Arguments
  - `::Type{T}`: The type of the constant to retrieve. Supported types are `Float16`, `Float32`, and `Float64`. Defaults to `Float64`.

  # Returns
  - The corresponding constant value for the specified type `T`.

  See: [`set_minoraxis`](@ref)
  """

  @eval begin

    @doc $docstring function $fn end
    @generated function $fn(::Type{T}=Float64) where T
      if T==Float16
        quote
          return $($upper_const_symbol_F16)[]
        end
      elseif T==Float32
        quote
          return $($upper_const_symbol_F32)[]
        end
      elseif T==Float64
        quote
          return $($upper_const_symbol_F64)[]
        end
      else
        :(error("Unsupported type $T"))
      end
    end
  end
end






@generated function _unrolled_initialize_theta(px1::T,py1::T,θ₀::T;ϵ=1e-2) where T
  initial_quote= quote
    b=get_minoraxis(T)
    e²=get_e²(T)
    θ=mod2pi(θ₀)
    fold=0.0
  end

  loop_quote= quote

    cosθ= cos(θ)
    sinθ= sin(θ)

    fx=(cosθ-px1)
    fy=(b*sinθ-py1)

    f= fx*fx+fy*fy

    px2    = cosθ
    py2    = b*sinθ

    dpx2dθ = -sinθ
    dpy2dθ = b*cosθ

    d²px2dθ² = -px2
    d²py2dθ² = -py2

    dpdθ_squared= 1- e²*cosθ*cosθ

    g = fx*dpx2dθ+fy*dpy2dθ

    h = dpdθ_squared+fx*d²px2dθ²+fy*d²py2dθ²

    p=-g/h
    θ+=p

    if abs(f-fold)<ϵ
      return mod2pi(θ)
    end
    fold=f
  end
  unrolled_loop=Expr[]

  for k in 1:10
    push!(unrolled_loop,loop_quote)
  end
  if T<:IEEEFloat
    return quote
      $(initial_quote)
      $(unrolled_loop...)
      return mod2pi(θ)
    end
  end
end
"""
  unrolled_initialize_theta(px1,px2,θ₀=0.0)

Given the normalized coordinates of the satellite, and an optional initial angle (default=0.0).
The function generates the inital θ to speedup the raytracing

"""
unrolled_initialize_theta(px::T,py::T,θ₀::T) where T<:IEEEFloat =_unrolled_initialize_theta(px,py,θ₀)
unrolled_initialize_theta(px1::I,py1::I,θ₀::I) where I<:Integer = unrolled_initialize_theta(float.(px1,py1,θ₀)...)
unrolled_initialize_theta(px1,py1,θ₀) =unrolled_initialize_theta(promote(px1,py1,θ₀)...)
unrolled_initialize_theta(px,py) =unrolled_initialize_theta(px,py,0)
unrolled_initialize_theta(args::A) where A<:NTuple{3,T} where T  =unrolled_initialize_theta(args[1],args[2],args[3])



v=Vector{Float64}(undef,10^5)
v1=similar(v)
init_points=[Tuple(x) for x in eachcol(rand(3,10^5))]

@benchmark map!(x-> unrolled_initialize_theta(x...),$v,$init_points)

@eval @generated function _unrolled_minimization_distance_from_ellipse(::Val{Niter},t::T,θ::T, s::T, px1::T, py1::T, dx1::T, dy1::T,θmin,θmax; δ=1e-15)::NTuple{3,T} where {T,Niter}

    if Niter<=0
      return :(error("Minimum number of loop is 1, kmax  provided was $Niter "))
    end

    exps=Expr[]

    for i in 1:10
      push!(exps,return_loop)
    end


    Nloops,Rem=divrem(Niter,10)


    if Niter>10
      loop_quote=quote
        for i in 1:div(Niter,10)
          $(exps...)
        end
      end
      push!(exps,loop_quote)
    end

    if rem(Niter,10)≠0
      for i in 1:Rem
        push!(exps,return_loop)
      end
    end

    return quote
      $(initial_code_quote)
      if abs(f)<δ
        return (t,θ,T(0))
      end
      $(Expr(:block,exps...))

      $(clamped_f)
      return (t,θ, s1-s)
    end
end

@eval  @generated function minimization_distance_ellipse!(t_out::A,θ_out::A,s_out::A,apx::A,apy::A,adx::A,ady::A,s_bottom::A,s_top::A,aθmin::A,aθmax::A,ascending::AB;δ=1e-10,kmax::Int=30, kwargs...)::Nothing where {A<:AbstractArray{T},AB<:AbstractArray{Bool}} where T<:IEEEFloat
    return quote
      @assert size(adx)==size(ady) "direction coordinates has to be the same size"
      @assert size(apx)==size(apy) "position coordinates has to be the same size"
      @assert size(apx)==size(adx) "position and direction coordinates has to be the same size"
      @assert size(θ_out) ==size(t_out)  "the optimization parameters t and θ had to be the same size"
      @assert size(t_out) ==size(adx) "optimization parameter and problem size should be the same"
      @assert  size(t_out)==size(s_out)   "optimization parameter and problem size should be the same"
      @assert δ>0 "δ has to be positive"
      @inbounds for i in eachindex(t_out)
        local px1=apx[i]
        local py1=apy[i]
        local dx1=adx[i]
        local dy1=ady[i]
        local θ=θ_out[i]
        # check if the ray is descending or ascending
        local s= ascending[i]==false ? s_bottom[i] : s_top[i]
        # If ascending and Δf>0, then the real position of the tangent quote would be s-f (below the expected level)
        # If descending and Δf>0, then the real position of the tangent quote would be s+f (above the expected level)
        local directional_sign= ascending[i] ? -1.0 : 1.0
        #############################################

        local θmin=aθmin[i]
        local θmax=aθmax[i]
        $(_unrolled_minimization_distance_from_ellipse²_inplace_quote)

        if t<0 && ascending[i]==false
          ascending[i]=true
          θ=θ_out[i]
          s=s_top[i]
          $(_unrolled_minimization_distance_from_ellipse²_inplace_quote)
        end

        t_out[i]=t
        θ_out[i]=θ
        s_out[i]=s
      end
    end
end

unrolled_minimization_distance_from_ellipse(t::T,θ::T,s::T,px::T,py::T,dx::T,dy::T;kmax::Int=18,θmin=-Inf,θmax=Inf) where T<:IEEEFloat=@inline  _unrolled_minimization_distance_from_ellipse(Val{kmax}(),t,θ,s,px,py,dx,dy,θmin,θmax)
unrolled_minimization_distance_from_ellipse(t,θ,s,px,py,dx,dy;kwargs...)=@inline unrolled_minimization_distance_from_ellipse(promote(t,θ,s,px,py,dx,dy)...;kwargs...)
unrolled_minimization_distance_from_ellipse(x::A;kwargs...) where A<:NTuple{7} =@inline unrolled_minimization_distance_from_ellipse(x[1],x[2],x[3],x[4],x[5],x[7],x[7];kwargs...)

function distance_from_ellipse!(t::A,θ::A,px::A,py::A,dx::A,dy::A;kmin::Int=10,θmin=-Inf,θmax=Inf, kwargs...)::Nothing where A<:AbstractArray{T} where T<:IEEEFloat
  #  @assert size(dx)==size(dy) "direction coordinates has to be the same size"
  #  @assert size(px)==size(py) "position coordinates has to be the same size"
  #  @assert size(px)==size(dx) "position and direction coordinates has to be the same size"
  #  @assert size(θ) ==size(t)  "the optimization parameters t and θ had to be the same size"
  #  @assert size(t) ==size(dx) "optimization parameter and problem size should be the same"
  s₀=1.0
  @inbounds for i in eachindex(t)

    local px_i = px[i]
    local py_i = py[i]
    local dx_i = dx[i]
    local dy_i = dy[i]
    local θ_i = θ[i]
    # Prepare a Ref to store the result:


    θ_i=@inline _unrolled_initialize_theta(px_i,py_i,θ_i)
     #@debug "i:$i θ: $(θ[i]) px: $(px[i]) py: $(px[i])  dx: $(dx[i]) dy: $(dx[i]) "

    (t_i,θ_i,f)= _unrolled_minimization_distance_from_ellipse(Val{kmin}(),0.0,θ_i,s₀,
      px_i,py_i,dx_i,dy_i,θmin,θmax)

    px[i]=px_i
    py[i]=py_i
    dx[i]=dx_i
    dy[i]=dy_i
    θ[i]=θ_i
    t[i]=t_i

  end
  return nothing

end


@eval function test_unrolled_initialize_theta(px1::T,py1::T,θ₀::T;ϵ=1e-2) where T
  initial_quote= quote
    b=get_minoraxis(T)
    e²=get_e²(T)
    θ=mod2pi(θ₀)
    fold=0.0
  end

  loop_quote= quote

    cosθ= cos(θ)
    sinθ= sin(θ)

    fx=(cosθ-px1)
    fy=(b*sinθ-py1)

    f= fx*fx+fy*fy

    px2    = cosθ
    py2    = b*sinθ

    dpx2dθ = -sinθ
    dpy2dθ = b*cosθ

    d²px2dθ² = -px2
    d²py2dθ² = -py2

    dpdθ_squared= 1- e²*cosθ*cosθ

    g = fx*dpx2dθ+fy*dpy2dθ

    h = dpdθ_squared+fx*d²px2dθ²+fy*d²py2dθ²

    p=-g/h
    θ+=p

    if abs(f-fold)<ϵ
      return mod2pi(θ)
    end
    fold=f
  end
  unrolled_loop=Expr[]

  for k in 1:10
    push!(unrolled_loop,loop_quote)
  end
  if T<:IEEEFloat
    return quote
      $(initial_quote)
      $(unrolled_loop...)
      return mod2pi(θ)
    end
  end
end

#=
@eval  function test_minimization_distance_ellipse!(t_out::A,θ_out::A,s_out::A,apx::A,apy::A,adx::A,ady::A,s_bottom::A,s_top::A,aθmin::A,aθmax::A,ascending::AB;δ=1e-10,kmax::Int=30, kwargs...) where {A<:AbstractArray{T},AB<:AbstractArray{Bool}} where T<:IEEEFloat
  @assert size(adx)==size(ady) "direction coordinates has to be the same size"
  @assert size(apx)==size(apy) "position coordinates has to be the same size"
  @assert size(apx)==size(adx) "position and direction coordinates has to be the same size"
  @assert size(θ_out) ==size(t_out)  "the optimization parameters t and θ had to be the same size"
  @assert size(t_out) ==size(adx) "optimization parameter and problem size should be the same"
  @assert  size(t_out)==size(s_out)   "optimization parameter and problem size should be the same"
  @assert δ>0 "δ has to be positive"
  @inbounds for i in eachindex(t_out)
    local px1=apx[i]
    local py1=apy[i]
    local dx1=adx[i]
    local dy1=ady[i]
    local θ=θ_out[i]

    # check if the ray is descending or ascending
    local s= ascending[i]==false ? s_bottom[i] : s_top[i]

    # If ascending and Δf>0, then the real position of the tangent quote would be s-f (below the expected level)
    # If descending and Δf>0, then the real position of the tangent quote would be s+f (above the expected level)
    local directional_sign= ascending[i] ? -1.0 : 1.0
    #############################################

      local θmin=aθmin[i]
      local θmax=aθmax[i]
      $(_unrolled_minimization_distance_from_ellipse²_inplace_quote)

      if t<0 && ascending[i]==false
        ascending[i]=true
        θ=θ_out[i]
        s=s_top[i]
        $(_unrolled_minimization_distance_from_ellipse²_inplace_quote)
      end

      t_out[i]=t
      θ_out[i]=θ
      s_out[i]=s
    end
  end
end
=#

#=
test_minimization_distance_ellipse!(t_out,θ_out,s_out,px,py,dx,dy,s_bottom,s_top,θmin,θmax,ascending)

MacroTools.striplines(test_minimization_distance_ellipse!(t_out,θ_out,s_out,px,py,dx,dy,s_bottom,s_top,θmin,θmax,ascending)
)
_unrolled_initialize_theta(1.0,1.0,0.0;ϵ=1e-2)
open("generated_fast_code/generated_initialization_theta.jl","w") do io
  print(io,"function _unrolled_initialize_theta(px1::T,px2::T,θ::T;ϵ=1e-2) where T\n")
  println(io,MacroTools.striplines(test_unrolled_initialize_theta(1.0,1.0,0.0;ϵ=1e-2)))
  println(io,"end")
end
=#

include("../generated_fast_code/generated_fast_implementation.jl")
