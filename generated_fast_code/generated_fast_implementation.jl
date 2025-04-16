using StructArrays
using Base: IEEEFloat
using BenchmarkTools
using Polyester

abstract type AbstractAtmosphere{T<:IEEEFloat} end
@kwdef struct AtmosphereProfile2D{T} <: AbstractAtmosphere{T}
    refraction_index_ave::T
    temperature_ave::T
    pressure_ave::T
    θ_left::T                # I do not need the right becaue cause I can compute them by looking at the neighbors
    s_top::T                 # I do not need the bottom cause I can compute them by looking at the upper neighbors
    AtmosphereProfile2D{T}(a::T,b::T,c::T,d::T,e::T) where T = new{T}(a,b,c,d,e)
    AtmosphereProfile2D(a::T,b::T,c::T,d::T,e::T) where T = new{T}(a,b,c,d,e)
end

function Base.show(io::IO, ::AtmosphereProfile2D{T}) where T
  print(io, "AtmosphereProfile2D{$T}")
end

function Base.show(io::IO, ::MIME"text/plain", p::AtmosphereProfile2D{T}) where T
  print(io, "AtmosphereProfile2D{$T}")
  println(io, "  n: $(p.refraction_index_ave)")
  println(io, "  T: $(p.temperature_ave) K")
  println(io, "  P: $(p.pressure_ave) hPa")
  println(io, "  θ: $(p.θ_left) rad")
  println(io, "  h: $(p.s_top)")
end

# Inlining the function might be better than using semicircular matrix for an effective code in particular if I can really inline it
@inline function get_semicircular_index(i::Int, j::Int, N::Int, M::Int)
      local i_1 = mod1(i, N)    # redial index
      local j_1 = j<=M ? M : -1   # linear index  (so if it is <=0 I know I am out of bound, and if it is negative I am at the end)
      return i_1, j_1
end

@kwdef struct InputRay{T<:IEEEFloat}
  point_x::T
  point_y::T
  direction_x::T
  direction_y::T
  n ::T=1.0
  θmin::T=-Inf
  θmax::T=-Inf
  ascending::Bool=false
  InputRay{T}(point_x::T,point_y::T,direction_x::T,direction_y::T,n::T,θmin::T,θmax::T) where T= new{T}(point_x,point_y,direction_x,direction_y,n,θmin,θmax,false)
  InputRay(point_x::T,point_y::T,direction_x::T,direction_y::T,n::T,θmin::T,θmax::T) where T   = new{T}(point_x,point_y,direction_x,direction_y,n,θmin,θmax,false)
end

function Base.show(io::IO, p::InputRay{T}) where T
   print(io, "$((round(p.point_x;sigdigits=3),round(p.point_y;sigdigits=3))) ")
end

function Base.show(io::IO,  ::MIME"text/plain",p::InputRay{T}) where T
  point_x=round(p.point_x;sigdigits=3)
  point_y=round(p.point_y;sigdigits=3)
  direction_x=round(p.direction_x;sigdigits=3)
  direction_y=round(p.direction_y;sigdigits=3)
  n_minus_1 =round(p.n-1;sigdigits=3)
  θmin= p.θmin
  θmax= p.θmax

  if p.ascending
    print(io, "Ascending ")
  end
  println(io, "Ray{$(T)}")
  println(io, "  point      : $((point_x,point_y))")
  println(io, "  direction  : $((direction_x,direction_y))")
  println(io, "  (n-1)      : $(n_minus_1)")
  if !isnothing(θmin) && isfinite(θmin)
    println(io, "  θmin       : $(θmin)")
  end
  if !isnothing(θmax) && isfinite(θmax)
    println(io, "  θmax       : $(θmax)")
  end
end


@kwdef struct ResultsRayTracing{T<:IEEEFloat}
  i::Int=-1
  j::Int=-1
  θ::T
  t::T
  h::T
  n::T=1.0
  point_x::T
  point_y::T
  direction_x::T
  direction_y::T
  ResultsRayTracing{T}(i::Int,j::Int,θ::T,t::T,h::T,n::T,point_x::T,point_y::T,direction_x::T,direction_y::T) where T= new{T}(i,j,θ,t,h,n,point_x,point_y,direction_x,direction_y)
  ResultsRayTracing(i::Int,j::Int,θ::T,t::T,h::T,n::T,point_x::T,point_y::T,direction_x::T,direction_y::T) where T   = new{T}(i,j,θ,t,h,n,point_x,point_y,direction_x,direction_y)
end

function Base.show(io::IO, p::ResultsRayTracing{T}) where T
  print(io, "$((round(p.point_x;sigdigits=3),round(p.point_y;sigdigits=3))) ")
end

function Base.show(io::IO,  ::MIME"text/plain",p::ResultsRayTracing{T}) where T
  point_x=round(p.point_x;sigdigits=3)
  point_y=round(p.point_y;sigdigits=3)
  direction_x=round(p.direction_x;sigdigits=3)
  direction_y=round(p.direction_y;sigdigits=3)
  h =round(p.h;sigdigits=3)
  t =round(p.t;sigdigits=3)
  θ =round(p.θ;sigdigits=3)
  n_minus_1 =round(p.n-1;sigdigits=3)

  println(io, "RayTracing{$(T)}")
  println(io, "clove     : $((p.i,p.j))")
  println(io, "angle     : $(θ)")
  println(io, "height    : $(h)")
  println(io, "length    : $(t)")
  println(io, "(n-1)     : $(n_minus_1)")
  println(io, "point     : $((point_x,point_y))")
  println(io, "direction : $((direction_x,direction_y))")
end


function initialize_theta(θin::T,point_x::T,point_y::T;ϵ=1e-2,kmax_init::Int=30)::T where T<:IEEEFloat
  @assert ϵ > 0 "δ has to be positive"
  @assert kmax_init > 0 "kmax_init has to be positive"

  local θ   = θin
  local Px  = point_x
  local Py  = point_y
  local fold= 0.0
  local f   = 0.0
  local b_normalized = get_minoraxis(T)
  local e² = get_e²(T)
  local hasconverged = false
  # Initialization function
  begin
    begin
      cosθ = cos(θ)
      sinθ = sin(θ)
      bcosθ = b_normalized * cosθ
      bsinθ = b_normalized * sinθ
      cosθ² = cosθ * cosθ
    end
    point_x2    = cosθ
    point_y2    = bsinθ

    Fx= point_x2
    Fy= point_y2

    begin
      fx=(Fx-Px)
      fy=(Fy-Py)
      f= fx*fx+fy*fy
    end
  end
  fold=f
  # Newton loop
  for k in 1:kmax_init
    begin
      begin
        dpoint_x2dθ = -sinθ
        dpoint_y2dθ = bcosθ

        d²point_x2dθ² = -point_x2
        d²point_y2dθ² = -point_y2

        dpdθ_squared= 1- e²*cosθ*cosθ

        g = fx*dpoint_x2dθ+fy*dpoint_y2dθ
        h = dpdθ_squared+fx*d²point_x2dθ²+fy*d²point_y2dθ²
        p = -g/h
      end
      θ = mod2pi(θ+p)
      begin
        begin
          cosθ = cos(θ)
          sinθ = sin(θ)
          bcosθ = b_normalized * cosθ
          bsinθ = b_normalized * sinθ
          cosθ² = cosθ * cosθ
        end
        point_x2    = cosθ
        point_y2    = bsinθ

        Fx= point_x2
        Fy= point_y2

        begin
          fx=(Fx-Px)
          fy=(Fy-Py)
          f= fx*fx+fy*fy
        end
      end
      if abs(f-fold)<ϵ
        hasconverged = true
        @debug "Converged after $k iterations with Δf: $(abs(f-fold))"
        return θ
      end
      fold=f
    end
  end
  @debug "Did not converged after $kmax_init iterations with Δf: $(abs(f-fold))"
  return θ

end

# Generated at: 2021-09-30T14:00:00.000
# High performance code obtained by interpolating the code from the quoted functions and the generated functions
# Manual inlining of functions and loops
# The code is optimized for the input parameters
 ##########################################################################################
"""
    initialize_theta!(θ_out::A,apoint_x1::A,apoint_y2::A;ϵ=1e-2,kmax_init::Int=30) where A<:AbstractArray{T} where T

  Initialize the initial angle theta as the approximate geocentric angle between the satellite and Earth

  # Arguments
    - `θ_out     <:Array{IEEEFloat}` : the optimization parameter θ (astromonical angle of the Earth)
    - `apoint_x       <:Array{IEEEFloat}` : the x position of the ray
    - `apoint_y       <:Array{IEEEFloat}` : the y position of the ray

  # Optional Arguments
    - `ϵ::IEEEFloat=1e-10` : the tolerance for the minimum distance
    - `kmax::Int=20`       : the maximum number of iterations for the Newton method

  See also: [`fast_ray_bending`](@ref) [`fast_initization_theta`](@ref) [`fast_ray_tracing`](@ref)
"""
function initialize_theta!(θ_out::A,apoint_x1::A,apoint_y2::A;ϵ=1e-2,kmax_init::Int=30) where A<:AbstractArray{T} where T
  @assert size(apoint_x) == size(apoint_y)   "position coordinates has to be the same size"
  @assert size(θ_out) == size(apoint_x) "the optimization parameters θ has to be the same size"

  @assert δ > 0 "δ has to be positive"
  @assert kmax_init > 0 "kmax_init has to be positive"

  b_normalized = get_minoraxis(T)
  e² = get_e²(T)
  kinitloops = div(kmax, 5,RoundUp)
  @inbounds for i = eachindex(t_out)

    local fold=0.0
    local f   = 0.0
    local Px  = apoint_x1[i]
    local Py  = apoint_y1[i]
    local θ   = atan(Py,Px)  # start with the angle of the ray as an initial guess
    local Fx  = 0.0
    local Fy  = 0.0
    local fx  = 0.0
    local fy  = 0.0
    begin
      begin
        cosθ = cos(θ)
        sinθ = sin(θ)
        bcosθ = b_normalized * cosθ
        bsinθ = b_normalized * sinθ
        cosθ² = cosθ * cosθ
      end
      point_x2    = cosθ
      point_y2    = bsinθ

      Fx= point_x2
      Fy= point_y2

      begin
        fx=(Fx-Px)
        fy=(Fy-Py)
        f= fx*fx+fy*fy
      end
    end
    for kin in 1:kinitloops
      # step block 1
      begin
        begin
          dpoint_x2dθ = -sinθ
          dpoint_y2dθ = bcosθ

          d²point_x2dθ² = -point_x2
          d²point_y2dθ² = -point_y2

          dpdθ_squared= 1- e²*cosθ*cosθ

          g = fx*dpoint_x2dθ+fy*dpoint_y2dθ
          h = dpdθ_squared+fx*d²point_x2dθ²+fy*d²point_y2dθ²
          p = -g/h
        end
        θ = mod2pi(θ+p)
        begin
          begin
            cosθ = cos(θ)
            sinθ = sin(θ)
            bcosθ = b_normalized * cosθ
            bsinθ = b_normalized * sinθ
            cosθ² = cosθ * cosθ
          end
          point_x2    = cosθ
          point_y2    = bsinθ

          Fx= point_x2
          Fy= point_y2

          begin
            fx=(Fx-Px)
            fy=(Fy-Py)
            f= fx*fx+fy*fy
          end
        end
        if abs(f-fold)<ϵ
          break
        end
        fold=f
      end
      # step block 2
      begin
        begin
          dpoint_x2dθ = -sinθ
          dpoint_y2dθ = bcosθ

          d²point_x2dθ² = -point_x2
          d²point_y2dθ² = -point_y2

          dpdθ_squared= 1- e²*cosθ*cosθ

          g = fx*dpoint_x2dθ+fy*dpoint_y2dθ
          h = dpdθ_squared+fx*d²point_x2dθ²+fy*d²point_y2dθ²
          p = -g/h
        end
        θ = mod2pi(θ+p)
        begin
          begin
            cosθ = cos(θ)
            sinθ = sin(θ)
            bcosθ = b_normalized * cosθ
            bsinθ = b_normalized * sinθ
            cosθ² = cosθ * cosθ
          end
          point_x2    = cosθ
          point_y2    = bsinθ

          Fx= point_x2
          Fy= point_y2

          begin
            fx=(Fx-Px)
            fy=(Fy-Py)
            f= fx*fx+fy*fy
          end
        end
        if abs(f-fold)<ϵ
          break
        end
        fold=f
      end
      # step block 3
      begin
        begin
          dpoint_x2dθ = -sinθ
          dpoint_y2dθ = bcosθ

          d²point_x2dθ² = -point_x2
          d²point_y2dθ² = -point_y2

          dpdθ_squared= 1- e²*cosθ*cosθ

          g = fx*dpoint_x2dθ+fy*dpoint_y2dθ
          h = dpdθ_squared+fx*d²point_x2dθ²+fy*d²point_y2dθ²
          p = -g/h
        end
        θ = mod2pi(θ+p)
        begin
          begin
            cosθ = cos(θ)
            sinθ = sin(θ)
            bcosθ = b_normalized * cosθ
            bsinθ = b_normalized * sinθ
            cosθ² = cosθ * cosθ
          end
          point_x2    = cosθ
          point_y2    = bsinθ

          Fx= point_x2
          Fy= point_y2

          begin
            fx=(Fx-Px)
            fy=(Fy-Py)
            f= fx*fx+fy*fy
          end
        end
        if abs(f-fold)<ϵ
          break
        end
        fold=f
      end
      # step block 4
      begin
        begin
          dpoint_x2dθ = -sinθ
          dpoint_y2dθ = bcosθ

          d²point_x2dθ² = -point_x2
          d²point_y2dθ² = -point_y2

          dpdθ_squared= 1- e²*cosθ*cosθ

          g = fx*dpoint_x2dθ+fy*dpoint_y2dθ
          h = dpdθ_squared+fx*d²point_x2dθ²+fy*d²point_y2dθ²
          p = -g/h
        end
        θ = mod2pi(θ+p)
        begin
          begin
            cosθ = cos(θ)
            sinθ = sin(θ)
            bcosθ = b_normalized * cosθ
            bsinθ = b_normalized * sinθ
            cosθ² = cosθ * cosθ
          end
          point_x2    = cosθ
          point_y2    = bsinθ

          Fx= point_x2
          Fy= point_y2

          begin
            fx=(Fx-Px)
            fy=(Fy-Py)
            f= fx*fx+fy*fy
          end
        end
        if abs(f-fold)<ϵ
          break
        end
        fold=f
      end
      # step block 5
      begin
        begin
          dpoint_x2dθ = -sinθ
          dpoint_y2dθ = bcosθ

          d²point_x2dθ² = -point_x2
          d²point_y2dθ² = -point_y2

          dpdθ_squared= 1- e²*cosθ*cosθ

          g = fx*dpoint_x2dθ+fy*dpoint_y2dθ
          h = dpdθ_squared+fx*d²point_x2dθ²+fy*d²point_y2dθ²
          p = -g/h
        end
        θ = mod2pi(θ+p)
        begin
          begin
            cosθ = cos(θ)
            sinθ = sin(θ)
            bcosθ = b_normalized * cosθ
            bsinθ = b_normalized * sinθ
            cosθ² = cosθ * cosθ
          end
          point_x2    = cosθ
          point_y2    = bsinθ

          Fx= point_x2
          Fy= point_y2

          begin
            fx=(Fx-Px)
            fy=(Fy-Py)
            f= fx*fx+fy*fy
          end
        end
        if abs(f-fold)<ϵ
          break
        end
        fold=f
      end
    end
    θ_out[i] = θ
  end
end

"""
    fast_minimization_distance!(t_out::A,θ_out::A,s_out::A,apoint_x::A,apoint_y::A,adirection_x::A,adirection_y::A,aθmin::A,aθmax::A,ascending::AB;δ=1e-10,kmax::Int=30, kwargs...) where {A<:AbstractArray{T},AB<:AbstractArray{Bool}} where T<:IEEEFloat

  Compute the minimum distance between a ray and an ellipse in 2D and modifies in place the optimization parameters `t_out`, `θ_out`, and `s_out` for each ray.
  It also updates the position apoint_x, apoint_y given the new t

  # Arguments
    - `t_out     <:Array{IEEEFloat}` : the optimization parameter t (length of the ray)
    - `θ_out     <:Array{IEEEFloat}` : the optimization parameter θ (astromonical angle of the Earth)
    - `s_out     <:Array{IEEEFloat}` : the optimization parameter s (tangent quote of the ray)
    - `apoint_x       <:Array{IEEEFloat}` : the x position of the ray
    - `apoint_y       <:Array{IEEEFloat}` : the y position of the ray
    - `adirection_x       <:Array{IEEEFloat}` : the x direction of the ray
    - `adirection_y       <:Array{IEEEFloat}` : the y direction of the ray
    - `aθmin     <:Array{IEEEFloat}` : the left angle of the wedge
    - `aθmax     <:Array{IEEEFloat}` : the right angle of the wedge
    - `ascending <:Array{Bool}`      : true defines the ascending part of the ray, finding the first true given the tangent quote

  # Optional Arguments
    - `δ::IEEEFloat=1e-10` : the tolerance for the minimum distance
    - `kmax::Int=30`       : the maximum number of iterations for the Newton method

  See also: [`fast_ray_bending`](@ref) [`fast_initization_theta`](@ref) [`fast_ray_tracing`](@ref)
"""
function fast_minimization_distance!(t_out::A,θ_out::A,s_out::A,apoint_x::A,apoint_y::A,adirection_x::A,adirection_y::A,aθmin::A,aθmax::A,ascending::AB;δ=1e-10,kmax::Int=30, kwargs...) where {A<:AbstractArray{T},AB<:AbstractArray{Bool}} where T<:IEEEFloat
  @assert size(adirection_x) == size(adirection_y) "direction coordinates has to be the same size"
  @assert size(apoint_x) == size(apoint_y) "position coordinates has to be the same size"
  @assert size(apoint_x) == size(adirection_x) "position and direction coordinates has to be the same size"
  @assert size(θ_out) == size(t_out) "the optimization parameters t and θ had to be the same size"
  @assert size(t_out) == size(adirection_x) "optimization parameter and problem size should be the same"
  @assert size(t_out) == size(s_out) "optimization parameter and problem size should be the same"
  @assert δ > 0 "δ has to be positive"

  b_normalized = get_minoraxis(T)
  e² = get_e²(T)
  kloops = div(kmax, 5,RoundUp)
  @inbounds for i = eachindex(t_out)
    local point_x1 = apoint_x[i]
    local point_y1 = apoint_y[i]
    local direction_x1 = adirection_x[i]
    local direction_y1 = adirection_y[i]
    local θ = θ_out[i]
    local s = s_out[i]        # aims for a level s, returns the effective level found
    local directional_sign = ascending[i] ? T(-1.0) :  T(1.0)  # we correct the level altitude with the real one
    local θmin = aθmin[i]
    local θmax = aθmax[i]
    # normalization of the direction
    begin
      norm = hypot(direction_x1, direction_y1)
      direction_x1 /= norm
      direction_y1 /= norm
    end
    origin_times_direction = point_x1 * direction_x1 + point_y1 * direction_y1

    begin
      # trigonometric function needed for computing the function
      begin
        cosθ = cos(θ)
        sinθ = sin(θ)
        bcosθ = b_normalized * cosθ
        bsinθ = b_normalized * sinθ
        cosθ² = cosθ * cosθ
      end

      # computation point on the ellipse at distance s
      begin
        R = 1 - e² * cosθ²
        N = 1 / sqrt(R)
        point_x2 = cosθ
        point_y2 = bsinθ
        direction_x2 = bcosθ
        direction_y2 = sinθ
        # point on the normal of the ellipse
        Fx = point_x2 + s * direction_x2 * N
        Fy = point_y2 + s * direction_y2 * N
      end
      # compute t minimum  for θ and the relative distance squared f
      # dropped 1/2 so that the square root is the displacement
      t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
      # point on the ray
      begin
        Px = point_x1 + t * direction_x1
        Py = point_y1 + t * direction_y1
      end
      begin
        fx = Fx - Px
        fy = Fy - Py
        f = fx * fx + fy * fy
      end

      # if the distance is alreadirection_y small enough or the function is nan
      # it means I am either at the end point or that the function has reached the minimum
      if abs(f)<δ || isnan(f)
        t_out[i] = NaN
        θ_out[i] = NaN
        s_out[i] = NaN
        return nothing
      end

      # set f-> fold
      fold = f
    end
    # Newton loop

    for k in 1:kloops
      # unroll block of 5 iterations
      # unroll block
      #   new function value and evaluation breaking criteria
      #   step block 1
      begin
        # update k ← k + 1
        # compute new θ
        begin
          # trigonometric function needed for computing the gradient and hessian
          begin
            sinθ² = sinθ * sinθ
            half_sin2θ = sinθ * cosθ
            cos2θ = cosθ² - sinθ²
          end
          # compute the gradient
          begin
            begin
              dpoint_x2dθ = -sinθ
              dpoint_y2dθ = bcosθ
              dpoint_x2dθ_0 = -bsinθ
              dpoint_y2dθ_0 = cosθ
            end
            begin
              half_dR = e² * half_sin2θ
              half_dR² = half_dR * half_dR
              half_d²R = e² * cos2θ
              N² = N * N
              dNdθ = -half_dR
              d²Ndθ² = -half_d²R + 3 * N² * half_dR²
              ddirection_x2dθ = (dpoint_x2dθ_0 + N² * dNdθ * direction_x2) * N
              ddirection_y2dθ = (dpoint_y2dθ_0 + N² * dNdθ * direction_y2) * N
            end
            begin
              dFxdθ = dpoint_x2dθ + s * ddirection_x2dθ
              dFydθ = dpoint_y2dθ + s * ddirection_y2dθ
              dtdθ = direction_x1 * dFxdθ + direction_y1 * dFydθ
              dPxdθ = direction_x1 * dtdθ
              dPydθ = direction_y1 * dtdθ
              dfxdθ = dFxdθ - dPxdθ
              dfydθ = dFydθ - dPydθ
            end
            g = fx * dfxdθ + fy * dfydθ
          end
          # compute the hessian
          begin
            begin
              d²point_x2dθ² = -point_x2
              d²point_y2dθ² = -point_y2
              d²direction_x2dθ²_0 = -direction_x2
              d²direction_y2dθ²_0 = -direction_y2
            end
            begin
              d²direction_x2dθ² = (d²direction_x2dθ²_0 + (2 * dpoint_x2dθ_0 * dNdθ + d²Ndθ² * N² * direction_x2) * N²) * N
              d²direction_y2dθ² = (d²direction_y2dθ²_0 + (2 * dpoint_y2dθ_0 * dNdθ + d²Ndθ² * N² * direction_y2) * N²) * N
              d²Fxdθ² = d²point_x2dθ² + d²direction_x2dθ²
              d²Fydθ² = d²point_y2dθ² + d²direction_y2dθ²
            end
            begin
              d²tdθ² = direction_x1 * d²Fxdθ² + direction_y1 * d²Fydθ²
              d²Pxdθ² = direction_x1 * d²tdθ²
              d²Pydθ² = direction_y1 * d²tdθ²
              dfxdθ_squared = dfxdθ * dfxdθ
              dfydθ_squared = dfydθ * dfydθ
              d²fxdθ² = d²Fxdθ² - d²Pxdθ²
              d²fydθ² = d²Fydθ² - d²Pydθ²
            end
            h = dfxdθ_squared + dfydθ_squared + fx * d²fxdθ² + fy * d²fydθ²
            # insure positive definess of the hessian by adding a const
            # similar to how LDLT works to ensure positive definiteness of matrix
            h = abs(h) > 10 ^ -5 ? abs(h) :  10 ^ -5
          end

          # compute the newton step
          p = -g / h
          # compute next θ
          θ = mod2pi(θ + p)
        end
        # compute new (f,t)
        begin
          begin
            cosθ = cos(θ)
            sinθ = sin(θ)
            bcosθ = b_normalized * cosθ
            bsinθ = b_normalized * sinθ
            cosθ² = cosθ * cosθ
          end
          begin
            R = 1 - e² * cosθ²
            N = 1 / sqrt(R)
            point_x2 = cosθ
            point_y2 = bsinθ
            direction_x2 = bcosθ
            direction_y2 = sinθ
            Fx = point_x2 + s * direction_x2 * N
            Fy = point_y2 + s * direction_y2 * N
          end
          t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
          begin
            Px = point_x1 + t * direction_x1
            Py = point_y1 + t * direction_y1
          end
          begin
            fx = Fx - Px
            fy = Fy - Py
            f = fx * fx + fy * fy
          end
        end
        # stopping criteria and update fold ← f
        begin
          if abs(f - fold) < δ
              break
          end
          fold = f
        end
      end
      #   step block 2
      begin
        # update k ← k + 1
        # compute new θ
        begin
          # trigonometric function needed for computing the gradient and hessian
          begin
            sinθ² = sinθ * sinθ
            half_sin2θ = sinθ * cosθ
            cos2θ = cosθ² - sinθ²
          end
          # compute the gradient
          begin
            half_dR = e² * half_sin2θ
            half_dR² = half_dR * half_dR
            half_d²R = e² * cos2θ
            N² = N * N
            dNdθ = -half_dR
            d²Ndθ² = -half_d²R + 3 * N² * half_dR²
            dpoint_x2dθ = -sinθ
            dpoint_y2dθ = bcosθ
            dpoint_x2dθ_0 = -bsinθ
            dpoint_y2dθ_0 = cosθ
            ddirection_x2dθ = (dpoint_x2dθ_0 + N² * dNdθ * direction_x2) * N
            ddirection_y2dθ = (dpoint_y2dθ_0 + N² * dNdθ * direction_y2) * N
            dFxdθ = dpoint_x2dθ + s * ddirection_x2dθ
            dFydθ = dpoint_y2dθ + s * ddirection_y2dθ
            dtdθ = direction_x1 * dFxdθ + direction_y1 * dFydθ
            dPxdθ = direction_x1 * dtdθ
            dPydθ = direction_y1 * dtdθ
            dfxdθ = dFxdθ - dPxdθ
            dfydθ = dFydθ - dPydθ
            g = fx * dfxdθ + fy * dfydθ
          end
          # compute the hessian
          begin
            d²point_x2dθ² = -point_x2
            d²point_y2dθ² = -point_y2
            d²direction_x2dθ²_0 = -direction_x2
            d²direction_y2dθ²_0 = -direction_y2
            d²direction_x2dθ² = (d²direction_x2dθ²_0 + (2 * dpoint_x2dθ_0 * dNdθ + d²Ndθ² * N² * direction_x2) * N²) * N
            d²direction_y2dθ² = (d²direction_y2dθ²_0 + (2 * dpoint_y2dθ_0 * dNdθ + d²Ndθ² * N² * direction_y2) * N²) * N
            d²Fxdθ² = d²point_x2dθ² + d²direction_x2dθ²
            d²Fydθ² = d²point_y2dθ² + d²direction_y2dθ²
            d²tdθ² = direction_x1 * d²Fxdθ² + direction_y1 * d²Fydθ²
            d²Pxdθ² = direction_x1 * d²tdθ²
            d²Pydθ² = direction_y1 * d²tdθ²
            dfxdθ_squared = dfxdθ * dfxdθ
            dfydθ_squared = dfydθ * dfydθ
            d²fxdθ² = d²Fxdθ² - d²Pxdθ²
            d²fydθ² = d²Fydθ² - d²Pydθ²
            h = dfxdθ_squared + dfydθ_squared + fx * d²fxdθ² + fy * d²fydθ²
            # insure positive definess of the hessian by adding a const
            # similar to how LDLT works to ensure positive definiteness of matrix
            h = abs(h) > 10 ^ -5 ? abs(h) :  10 ^ -5
          end

          # compute the newton step
          p = -g / h
          # compute next θ
          θ = mod2pi(θ + p)
        end
        # compute new (f,t)
        begin
          begin
            cosθ = cos(θ)
            sinθ = sin(θ)
            bcosθ = b_normalized * cosθ
            bsinθ = b_normalized * sinθ
            cosθ² = cosθ * cosθ
          end
          begin
            R = 1 - e² * cosθ²
            N = 1 / sqrt(R)
            point_x2 = cosθ
            point_y2 = bsinθ
            direction_x2 = bcosθ
            direction_y2 = sinθ
            Fx = point_x2 + s * direction_x2 * N
            Fy = point_y2 + s * direction_y2 * N
          end
          t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
          begin
            Px = point_x1 + t * direction_x1
            Py = point_y1 + t * direction_y1
          end
          begin
            fx = Fx - Px
            fy = Fy - Py
            f = fx * fx + fy * fy
          end
        end
        # stopping criteria and update fold ← f
        begin
          if abs(f - fold) < δ
              break
          end
          fold = f
        end
      end
      #   step block 3
      begin
        # update k ← k + 1
        # compute new θ
        begin
          # trigonometric function needed for computing the gradient and hessian
          begin
            sinθ² = sinθ * sinθ
            half_sin2θ = sinθ * cosθ
            cos2θ = cosθ² - sinθ²
          end
          # compute the gradient
          begin
            half_dR = e² * half_sin2θ
            half_dR² = half_dR * half_dR
            half_d²R = e² * cos2θ
            N² = N * N
            dNdθ = -half_dR
            d²Ndθ² = -half_d²R + 3 * N² * half_dR²
            dpoint_x2dθ = -sinθ
            dpoint_y2dθ = bcosθ
            dpoint_x2dθ_0 = -bsinθ
            dpoint_y2dθ_0 = cosθ
            ddirection_x2dθ = (dpoint_x2dθ_0 + N² * dNdθ * direction_x2) * N
            ddirection_y2dθ = (dpoint_y2dθ_0 + N² * dNdθ * direction_y2) * N
            dFxdθ = dpoint_x2dθ + s * ddirection_x2dθ
            dFydθ = dpoint_y2dθ + s * ddirection_y2dθ
            dtdθ = direction_x1 * dFxdθ + direction_y1 * dFydθ
            dPxdθ = direction_x1 * dtdθ
            dPydθ = direction_y1 * dtdθ
            dfxdθ = dFxdθ - dPxdθ
            dfydθ = dFydθ - dPydθ
            g = fx * dfxdθ + fy * dfydθ
          end
          # compute the hessian
          begin
            d²point_x2dθ² = -point_x2
            d²point_y2dθ² = -point_y2
            d²direction_x2dθ²_0 = -direction_x2
            d²direction_y2dθ²_0 = -direction_y2
            d²direction_x2dθ² = (d²direction_x2dθ²_0 + (2 * dpoint_x2dθ_0 * dNdθ + d²Ndθ² * N² * direction_x2) * N²) * N
            d²direction_y2dθ² = (d²direction_y2dθ²_0 + (2 * dpoint_y2dθ_0 * dNdθ + d²Ndθ² * N² * direction_y2) * N²) * N
            d²Fxdθ² = d²point_x2dθ² + d²direction_x2dθ²
            d²Fydθ² = d²point_y2dθ² + d²direction_y2dθ²
            d²tdθ² = direction_x1 * d²Fxdθ² + direction_y1 * d²Fydθ²
            d²Pxdθ² = direction_x1 * d²tdθ²
            d²Pydθ² = direction_y1 * d²tdθ²
            dfxdθ_squared = dfxdθ * dfxdθ
            dfydθ_squared = dfydθ * dfydθ
            d²fxdθ² = d²Fxdθ² - d²Pxdθ²
            d²fydθ² = d²Fydθ² - d²Pydθ²
            h = dfxdθ_squared + dfydθ_squared + fx * d²fxdθ² + fy * d²fydθ²
            # insure positive definess of the hessian by adding a const
            # similar to how LDLT works to ensure positive definiteness of matrix
            h = abs(h) > 10 ^ -5 ? abs(h) :  10 ^ -5
          end

          # compute the newton step
          p = -g / h
          # compute next θ
          θ = mod2pi(θ + p)
        end
        # compute new (f,t)
        begin
          begin
            cosθ = cos(θ)
            sinθ = sin(θ)
            bcosθ = b_normalized * cosθ
            bsinθ = b_normalized * sinθ
            cosθ² = cosθ * cosθ
          end
          begin
            R = 1 - e² * cosθ²
            N = 1 / sqrt(R)
            point_x2 = cosθ
            point_y2 = bsinθ
            direction_x2 = bcosθ
            direction_y2 = sinθ
            Fx = point_x2 + s * direction_x2 * N
            Fy = point_y2 + s * direction_y2 * N
          end
          t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
          begin
            Px = point_x1 + t * direction_x1
            Py = point_y1 + t * direction_y1
          end
          begin
            fx = Fx - Px
            fy = Fy - Py
            f = fx * fx + fy * fy
          end
        end
        # stopping criteria and update fold ← f
        begin
          if abs(f - fold) < δ
              break
          end
          fold = f
        end
      end
      #   step block 4
      begin
        # update k ← k + 1
        # compute new θ
        begin
          # trigonometric function needed for computing the gradient and hessian
          begin
            sinθ² = sinθ * sinθ
            half_sin2θ = sinθ * cosθ
            cos2θ = cosθ² - sinθ²
          end
          # compute the gradient
          begin
            half_dR = e² * half_sin2θ
            half_dR² = half_dR * half_dR
            half_d²R = e² * cos2θ
            N² = N * N
            dNdθ = -half_dR
            d²Ndθ² = -half_d²R + 3 * N² * half_dR²
            dpoint_x2dθ = -sinθ
            dpoint_y2dθ = bcosθ
            dpoint_x2dθ_0 = -bsinθ
            dpoint_y2dθ_0 = cosθ
            ddirection_x2dθ = (dpoint_x2dθ_0 + N² * dNdθ * direction_x2) * N
            ddirection_y2dθ = (dpoint_y2dθ_0 + N² * dNdθ * direction_y2) * N
            dFxdθ = dpoint_x2dθ + s * ddirection_x2dθ
            dFydθ = dpoint_y2dθ + s * ddirection_y2dθ
            dtdθ = direction_x1 * dFxdθ + direction_y1 * dFydθ
            dPxdθ = direction_x1 * dtdθ
            dPydθ = direction_y1 * dtdθ
            dfxdθ = dFxdθ - dPxdθ
            dfydθ = dFydθ - dPydθ
            g = fx * dfxdθ + fy * dfydθ
          end
          # compute the hessian
          begin
            d²point_x2dθ² = -point_x2
            d²point_y2dθ² = -point_y2
            d²direction_x2dθ²_0 = -direction_x2
            d²direction_y2dθ²_0 = -direction_y2
            d²direction_x2dθ² = (d²direction_x2dθ²_0 + (2 * dpoint_x2dθ_0 * dNdθ + d²Ndθ² * N² * direction_x2) * N²) * N
            d²direction_y2dθ² = (d²direction_y2dθ²_0 + (2 * dpoint_y2dθ_0 * dNdθ + d²Ndθ² * N² * direction_y2) * N²) * N
            d²Fxdθ² = d²point_x2dθ² + d²direction_x2dθ²
            d²Fydθ² = d²point_y2dθ² + d²direction_y2dθ²
            d²tdθ² = direction_x1 * d²Fxdθ² + direction_y1 * d²Fydθ²
            d²Pxdθ² = direction_x1 * d²tdθ²
            d²Pydθ² = direction_y1 * d²tdθ²
            dfxdθ_squared = dfxdθ * dfxdθ
            dfydθ_squared = dfydθ * dfydθ
            d²fxdθ² = d²Fxdθ² - d²Pxdθ²
            d²fydθ² = d²Fydθ² - d²Pydθ²
            h = dfxdθ_squared + dfydθ_squared + fx * d²fxdθ² + fy * d²fydθ²
            # insure positive definess of the hessian by adding a const
            # similar to how LDLT works to ensure positive definiteness of matrix
            h = abs(h) > 10 ^ -5 ? abs(h) :  10 ^ -5
          end

          # compute the newton step
          p = -g / h
          # compute next θ
          θ = mod2pi(θ + p)
        end
        # compute new (f,t)
        begin
          begin
            cosθ = cos(θ)
            sinθ = sin(θ)
            bcosθ = b_normalized * cosθ
            bsinθ = b_normalized * sinθ
            cosθ² = cosθ * cosθ
          end
          begin
            R = 1 - e² * cosθ²
            N = 1 / sqrt(R)
            point_x2 = cosθ
            point_y2 = bsinθ
            direction_x2 = bcosθ
            direction_y2 = sinθ
            Fx = point_x2 + s * direction_x2 * N
            Fy = point_y2 + s * direction_y2 * N
          end
          t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
          begin
            Px = point_x1 + t * direction_x1
            Py = point_y1 + t * direction_y1
          end
          begin
            fx = Fx - Px
            fy = Fy - Py
            f = fx * fx + fy * fy
          end
        end
        # stopping criteria and update fold ← f
        begin
          if abs(f - fold) < δ
              break
          end
          fold = f
        end
      end
      #   step block 5
      begin
        # update k ← k + 1
        # compute new θ
        begin
          # trigonometric function needed for computing the gradient and hessian
          begin
            sinθ² = sinθ * sinθ
            half_sin2θ = sinθ * cosθ
            cos2θ = cosθ² - sinθ²
          end
          # compute the gradient
          begin
            half_dR = e² * half_sin2θ
            half_dR² = half_dR * half_dR
            half_d²R = e² * cos2θ
            N² = N * N
            dNdθ = -half_dR
            d²Ndθ² = -half_d²R + 3 * N² * half_dR²
            dpoint_x2dθ = -sinθ
            dpoint_y2dθ = bcosθ
            dpoint_x2dθ_0 = -bsinθ
            dpoint_y2dθ_0 = cosθ
            ddirection_x2dθ = (dpoint_x2dθ_0 + N² * dNdθ * direction_x2) * N
            ddirection_y2dθ = (dpoint_y2dθ_0 + N² * dNdθ * direction_y2) * N
            dFxdθ = dpoint_x2dθ + s * ddirection_x2dθ
            dFydθ = dpoint_y2dθ + s * ddirection_y2dθ
            dtdθ = direction_x1 * dFxdθ + direction_y1 * dFydθ
            dPxdθ = direction_x1 * dtdθ
            dPydθ = direction_y1 * dtdθ
            dfxdθ = dFxdθ - dPxdθ
            dfydθ = dFydθ - dPydθ
            g = fx * dfxdθ + fy * dfydθ
          end
          # compute the hessian
          begin
            d²point_x2dθ² = -point_x2
            d²point_y2dθ² = -point_y2
            d²direction_x2dθ²_0 = -direction_x2
            d²direction_y2dθ²_0 = -direction_y2
            d²direction_x2dθ² = (d²direction_x2dθ²_0 + (2 * dpoint_x2dθ_0 * dNdθ + d²Ndθ² * N² * direction_x2) * N²) * N
            d²direction_y2dθ² = (d²direction_y2dθ²_0 + (2 * dpoint_y2dθ_0 * dNdθ + d²Ndθ² * N² * direction_y2) * N²) * N
            d²Fxdθ² = d²point_x2dθ² + d²direction_x2dθ²
            d²Fydθ² = d²point_y2dθ² + d²direction_y2dθ²
            d²tdθ² = direction_x1 * d²Fxdθ² + direction_y1 * d²Fydθ²
            d²Pxdθ² = direction_x1 * d²tdθ²
            d²Pydθ² = direction_y1 * d²tdθ²
            dfxdθ_squared = dfxdθ * dfxdθ
            dfydθ_squared = dfydθ * dfydθ
            d²fxdθ² = d²Fxdθ² - d²Pxdθ²
            d²fydθ² = d²Fydθ² - d²Pydθ²
            h = dfxdθ_squared + dfydθ_squared + fx * d²fxdθ² + fy * d²fydθ²
            # insure positive definess of the hessian by adding a const
            # similar to how LDLT works to ensure positive definiteness of matrix
            h = abs(h) > 10 ^ -5 ? abs(h) :  10 ^ -5
          end

          # compute the newton step
          p = -g / h
          # compute next θ
          θ = mod2pi(θ + p)
        end
        # compute new (f,t)
        begin
          begin
            cosθ = cos(θ)
            sinθ = sin(θ)
            bcosθ = b_normalized * cosθ
            bsinθ = b_normalized * sinθ
            cosθ² = cosθ * cosθ
          end
          begin
            R = 1 - e² * cosθ²
            N = 1 / sqrt(R)
            point_x2 = cosθ
            point_y2 = bsinθ
            direction_x2 = bcosθ
            direction_y2 = sinθ
            Fx = point_x2 + s * direction_x2 * N
            Fy = point_y2 + s * direction_y2 * N
          end
          t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
          begin
            Px = point_x1 + t * direction_x1
            Py = point_y1 + t * direction_y1
          end
          begin
            fx = Fx - Px
            fy = Fy - Py
            f = fx * fx + fy * fy
          end
        end
        # stopping criteria and update fold ← f
        begin
          if abs(f - fold) < δ
              break
          end
          fold = f
        end
      end

    end
    # clamping angle to the range
    # if θmin <= θmax  thene θ ∈ [θmin, θmax]
    begin
      if (θmin <= θmax && θmin < θ < θmax) || (θmin > θmax && (θ > θmin || θ < θmax))

          s += sqrt(f) * directional_sign

      else



        # angle clamping f(θ,θmin,θmax)
        begin
          if θmin<=θmax
            θ=clamp(θ,θmin,θmax)
          else
            dmin = mod(θ - θmin, 2π)
            dmax = mod(θmax - θ, 2π)
            θ = dmin < dmax ? θmin : θmax
          end
          θ
        end
        #update (point_x1,point_x2,θ,t,s)
        begin


          point_x2 = cos(θ)
          point_y2 = sin(θ) * b_normalized
          direction_x2 = b_normalized * cos(θ)
          direction_y2 = sin(θ)
          normF = hypot(direction_x2, direction_y2)
          direction_x2 /= normF
          direction_y2 /= normF

          det = direction_x1 * direction_y2 - direction_y1 * direction_x2

          ϵ = 1.0e-10        # TO DO: maybe add this value as a kwargs
          s = NaN        # initial s to NaN
          t = NaN        # initial t to NaN
          if det > ϵ
            t = ((point_x2 - point_x1) * direction_y2 - (point_y2 - point_y1) * direction_x2) / det
            s = ((point_x2 - point_x1) * direction_y1 - (point_y2 - point_y1) * direction_x1) / det
          end

          if (s < 0) # it starts from the ellipse surface, so a negative value is not possible
            s = NaN
          end
        end
      end
    end
    # check if the direction of the ray has changed and if it did, set the ascending to true
    # this helps to find easily the tangent quote
    # next iteration will look for the ascending s instead of the descending one
    if t < 0 && ascending[i] == false
      ascending[i] = true
    end
    #update all the output for the current iteration
    begin
      t_out[i] = t
      θ_out[i] = θ
      s_out[i] = s
      # update the position of the ray
      apoint_x[i] = point_x1 + t * direction_x1
      apoint_y[i] = point_y1 + t * direction_y1
    end
  end
end

"""
  fast_ray_bending!(θ_out::A,adirection_x::A,adirection_y::A,nᵢ::A,nₜ::A,aθmin::A,aθmax::A,ascending::AB;δ=1e-10, kwargs...)  where {AB<:AbstractArray{Bool},A<:AbstractArray{T}} where T<:IEEEFloat

  Compute the ray bending from passing through the atmosphere. This function is a fast version of the ray bending function that modifies in-place
  the arrays without any allocation.

  # Arguments
    - `θ_out     <:Array{IEEEFloat}` : the optimization parameter θ (astromonical angle of the Earth)
    - `adirection_x       <:Array{IEEEFloat}` : the x direction of the ray
    - `adirection_y       <:Array{IEEEFloat}` : the y direction of the ray
    - `nᵢ        <:Array{IEEEFloat}` : the incident index of refraction
    - `nₜ        <:Array{IEEEFloat}` :  the transmitted index of refraction
    - `aθmin     <:Array{IEEEFloat}` : the left angle of the wedge
    - `aθmax     <:Array{IEEEFloat}` : the right angle of the wedge
    - `ascending <:Array{Bool}`      : true defines the ascending part of the ray, finding the first true given the tangent quote

  # Optional arguments
    - `δ::Float64=1e-10` : the machine epsilon used to compute the minimum distance around the wedge ends to avoid numerical instability

  See also: [`fast_ray_bending`](@ref) [`fast_initization_theta`](@ref) [`fast_ray_tracing`](@ref)
"""
function fast_ray_bending!(θ_out::A,adirection_x::A,adirection_y::A,nᵢ::A,nₜ::A,aθmin::A,aθmax::A,ascending::AB; δ=1e-10, kwargs...) where {AB<:AbstractArray{Bool},A<:AbstractArray{T}} where T<:IEEEFloat
  @assert size(adirection_x) == size(adirection_y) "direction coordinates has to be the same size"
  @assert size(apoint_x) == size(apoint_y) "position coordinates has to be the same size"
  @assert size(apoint_x) == size(adirection_x) "position and direction coordinates has to be the same size"
  @assert size(θ_out) == size(adirection_x) "optimization parameter and problem size should be the same"

  @assert size(θ_out) == size(nᵢ) "optimization parameter and problem size should be the same"
  @assert size(θ_out) == size(nₜ) "optimization parameter and problem size should be the same"
  @assert size(θ_out) == size(aθmin) "optimization parameter and problem size should be the same"
  @assert size(θ_out) == size(aθmax) "optimization parameter and problem size should be the same"
  @assert δ > 0 "δ has to be positive"
  @assert size(ascending) == size(θ_out) "optimization parameter and problem size should be the same"
  b=get_minoraxis(T)

  @inbounds for i = eachindex(θ_out)
    local isAscending    = ascending[i]
    local n_incident     = nᵢ[i]
    local n_transmitted  = nₜ[i]
    local direction_x1 = adirection_x[i]
    local direction_y1 = adirection_y[i]
    local θ = θ_out[i]  # note: theta is the angle of the earth normal not the angle of the ray, this is important to remember

    norm_ray_direction =  hypot(direction_x1, direction_y1)
    direction_x1 /= norm_ray_direction
    direction_y1 /= norm_ray_direction
    # if the difference between the top and bottom s is less than the machine tollerance
    # it means that the ray is either inside the wedge or the next wedge has the same atmosphere
    n_incident ≈ n_transmitted && continue

    # the normal to an ellipse can be computed from the normal to the earth
    begin
      direction_x2 = b*cos(normal_earth)
      direction_y2 = sin(normal_earth)
      norm = hypot(direction_x2, direction_y2)
      direction_x2 /= norm
      direction_y2 /= norm
      begin
        begin
          # check if it is intersecting a level or a ray
          if (θ==θmax || θ==θmin)
            # the normal to a ray is defined as (y,-x)
            direction_x2 = direction_y2
            direction_y2 = -tmp
          elseif isAscending
            # the normal is inwards
            direction_x2 = -direction_x2
            direction_y2 = -direction_y2
          end

          # both directions are alreadirection_y normalized
          local Nx=direction_x2
          local Ny=direction_y2

          local direcion_ray_x=direction_x1
          local direction_ray_y=direction_y1
          local n01=n_incident/n_transmitted
          local n01²=n01*n01
          local cosθ_incident=-(Nx*direcion_ray_x+Ny*direction_ray_y)
          local sinθ²_transmitted =n01²*(1-cosθ_indicent*cosθ_incident)

          # check if the ray is internally reflected
          # this most likely happens if there is an issue with the atmosphere or if the tangent quote
          # happens to be at a level.
          if sinθ²_transmitted ≤ 1
            direction_x1= n01*direction_x1+(n01*cosθ_incident-sqrt(1-sinθ²_transmitted))*Nx
            direction_y1= n01*direction_y1+(n01*cosθ_incident-sqrt(1-sinθ²_transmitted))*Ny
          else
            direction_x1-=2*cosθ_incident*Nx
            direction_y1-=2*cosθ_incident*Ny
          end
        end

        # update the direction of the ray
        begin
          local norm_new_ray_direction = hypot(direction_x1, direction_y1)
          direction_x1 /= norm_new_ray_direction
          direction_y1 /= norm_new_ray_direction
          adirection_x[i] = direction_x1
          adirection_y[i] = direction_y1
        end
      end
    end
  end
end
##################################################################################### =#
"""
  fast_ray_tracing!(t_out::A,θ_out::A,s_out::A,apoint_x::A,apoint_y::A,adirection_x::A,adirection_y::A,incident_refractive_index::A,aθmin::A,aθmax::A,ascending::AB,
  atm_n::M,atm_θ::V1,atm_h::V2,retrieval_i::RETi,retrieval_j::RETi,retrieval_n::RETf,retrieval_θ::RETf,retrieval_t::RETf,retrieval_h::RETf,
  retrieval_point_x::RETf,retrieval_point_y::RETf,retrieval_direction_x::RETf,retrieval_direction_y::RETf,
  tangent_quote::A;δ=1e-10,kmax::Int=30,intersection_max::Int=140,initialized=false,free_space::T=T(1), kwargs...
  ) where {V1<:AbstractVector, V2<:AbstractVector, RETi, RETf, M, A<:AbstractArray{T}, AB<:AbstractArray{Bool}} where T<:IEEEFloat


Compute the ray tracing from passing through the atmosphere. This function is a fast version of the ray tracing function that modifies in-place
the arrays without any allocation. This function is a combination of the `fast_ray_bending!` and `fast_ray_tracing!` functions. The function
is optimized to avoid any allocation and to be as fast as possible.

# Arguments

- `t_out     <:Array{IEEEFloat}` : the optimization parameter t (distance from the origin)
- `θ_out     <:Array{IEEEFloat}` : the optimization parameter θ (astromonical angle of the Earth)
- `s_out     <:Array{IEEEFloat}` : the optimization parameter s (distance from the origin)
- `apoint_x       <:Array{IEEEFloat}` : the x position of the ray
- `apoint_y       <:Array{IEEEFloat}` : the y position of the ray
- `adirection_x       <:Array{IEEEFloat}` : the x direction of the ray
- `adirection_y       <:Array{IEEEFloat}` : the y direction of the ray
- `incident_refractive_index <:Array{IEEEFloat}` : the incident refractive index
- `aθmin     <:Array{IEEEFloat}` : the left angle of the wedge
- `aθmax     <:Array{IEEEFloat}` : the right angle of the wedge
- `ascending <:Array{Bool}`      : true defines the ascending part of the ray, finding the first true given the tangent quote
- `atm_n     <:Array{Array{IEEEFloat}}` : the refractive index of the atmosphere
- `atm_θ     <:Array{IEEEFloat}` : the angle of the atmosphere
- `atm_h     <:Array{IEEEFloat}` : the height of the atmosphere

# Retrievals
- `retrieval_i <:Array{Int}` : the i index of the retrieval
- `retrieval_j <:Array{Int}` : the j index of the retrieval
- `retrieval_n <:Array{IEEEFloat}` : the refractive index of the retrieval
- `retrieval_θ <:Array{IEEEFloat}` : the angle of the retrieval
- `retrieval_t <:Array{IEEEFloat}` : the distance of the retrieval
- `retrieval_h <:Array{IEEEFloat}` : the height of the retrieval
- `retrieval_point_x <:Array{IEEEFloat}` : the x position of the retrieval
- `retrieval_point_y <:Array{IEEEFloat}` : the y position of the retrieval
- `retrieval_direction_x <:Array{IEEEFloat}` : the x direction of the retrieval
- `retrieval_direction_y <:Array{IEEEFloat}` : the y direction of the retrieval


"""
function fast_ray_tracing!(t_out::A,θ_out::A,s_out::A,apoint_x::A,apoint_y::A,adirection_x::A,adirection_y::A,incident_refractive_index::A,aθmin::A,aθmax::A,ascending::AB,
  atm_n::M,atm_θ::V1,atm_h::V2,retrieval_i::RETi,retrieval_j::RETi,retrieval_n::RETf,retrieval_θ::RETf,retrieval_t::RETf,retrieval_h::RETf,
  retrieval_point_x::RETf,retrieval_point_y::RETf,retrieval_direction_x::RETf,retrieval_direction_y::RETf,
  tangent_quote::A;δ=1e-10,kmax::Int=30,intersection_max::Int=140,initialized=false,free_space::T=T(1), kwargs...
  ) where {V1<:AbstractVector, V2<:AbstractVector, RETi, RETf, M, A<:AbstractArray{T}, AB<:AbstractArray{Bool}} where T<:IEEEFloat
  @assert size(adirection_x) == size(adirection_y) "Direction coordinates (adirection_x, adirection_y) must have the same size"
  @assert size(apoint_x) == size(apoint_y) "Position coordinates (apoint_x, apoint_y) must have the same size"
  @assert size(apoint_x) == size(adirection_x) "Position (apoint_x, apoint_y) and direction (adirection_x, adirection_y) coordinates must have the same size"
  @assert size(θ_out) == size(t_out) "Optimization parameters θ_out and t_out must have the same size"
  @assert size(t_out) == size(adirection_x) "Optimization parameter t_out and direction coordinates (adirection_x, adirection_y) must have the same size"
  @assert size(t_out) == size(s_out) "Optimization parameters t_out and s_out must have the same size"
  @assert size(t_out) == size(ascending) "Optimization parameter t_out and ascending array must have the same size"
  @assert size(θ_out) == size(incident_refractive_index) "Optimization parameter θ_out and incident refractive index must have the same size"
  @assert size(θ_out) == size(aθmin) "Optimization parameter θ_out and aθmin must have the same size"
  @assert size(θ_out) == size(aθmax) "Optimization parameter θ_out and aθmax must have the same size"
  @assert size(tangent_quote) == size(θ_out) "Optimization parameter θ_out and tangent_quote must have the same size"
  @assert δ > 0 "δ must be positive"
  @assert intersection_max > 1 "intersection_max must be greater than 1"

  @assert isa(retrieval_i,VecOrMat) "retrieval_i has to be 2 dimensional"
  @assert size(retrieval_i)==size(retrieval_j) "retrieval_i and retrieval_j has to be the same size"
  @assert size(retrieval_i)==size(retrieval_n) "retrieval_i and retrieval_n has to be the same size"
  @assert size(retrieval_i)==size(retrieval_θ) "retrieval_i and retrieval_θ has to be the same size"
  @assert size(retrieval_i)==size(retrieval_h) "retrieval_i and retrieval_h has to be the same size"
  @assert size(retrieval_i)==size(retrieval_t)  "retrieval_i and retrieval_t has to be the same size"

  # I don't need it if using StaticArrays but I do not want to make assumptions
  Natm_n=size(atm_n,1)
  Matm_n=size(atm_n,2)
  Natm_θ=size(atm_θ,1)
  Matm_h=size(atm_h,1)

  @assert (Natm_θ-1)<=Natm_n<=Natm_θ "atm_n has $Natm_n rows but it has to have either the length of atm_θ-1, $(Natm_θ-1),  or  $(Natm_θ) (periodic radial atmosphere)"

  IsPeriodic = Natm_n == Natm_θ

  @assert Matm_n==Matm_h-1 "atm_n has $Matm_n columns but it has to have the length of atm_h-1, $(Matm_h-1)"


  NumRays=prod(size(adirection_x))  # size of the rays
  NumRaysRetrieval= size(retrieval_i,1) # size of the retrieval rays
  @assert size(retrieval_i,1) == NumRays "retrieval has to be have the number of rows $(NumRaysRetrieval) equivalent to the number of rays ($NumRays)"


  Miter=size(retrieval_i,2) # number of iterations

  iter_eff= min(Miter-1,intersection_max)

  # most of the code is generated in the same way as the previous ones
  # the only difference is that it all happens in a single loop

  b_normalized = get_minoraxis(T)
  e² = get_e²(T)
  ϵ  = 1.0e-10
  # early stop condition
  number_rays_stopped=0

  max_altitude=atm_h[1]

  @debug "Starting the ray tracing"
  @debug "max_altitude is $max_altitude"
  @debug " extrema of atm_h is $(minimum(atm_h)) and $(maximum(atm_h))"

  rho_max = T(50)
  kloops = div(kmax,5, RoundUp) #Internal loop of the Newton method
  #external loop over all the points
    @inbounds for iter in 1:iter_eff

    # internal loop of ray tracing
      #@batch for idirection_x_rays in eachindex(t_out)
      for idirection_x_rays in eachindex(t_out)
        @debug "Ray $idirection_x_rays and iteration $iter"
        # first iteration is diffent from the rest
        # I need to set the initial value of theta
        # and find the initial wedge, also all incident refractive index are 1
        # also s_top does not exist yet
        local θ = θ_out[idirection_x_rays]                # gibberish the first iteration
        local t = t_out[idirection_x_rays]                # gibberish always
        local directional_sign = ascending[idirection_x_rays] ? T(-1.0) :  T(1.0)
        local point_x1 = apoint_x[idirection_x_rays]
        local point_y1 = apoint_y[idirection_x_rays]
        local direction_x1 = adirection_x[idirection_x_rays]
        local direction_y1 = adirection_y[idirection_x_rays]
        local θmin = aθmin[idirection_x_rays]                                  # gibberish the first iteration not necessary but I prefer to have it for consistency and debugging
        local θmax = aθmax[idirection_x_rays]                                  # gibberish the first iteration not necessary but I prefer to have it for consistency and debugging
        local nᵢ = incident_refractive_index[idirection_x_rays]                   # gibberish the first iteration
        local nₜ  = T(0)
        local isAscending = ascending[idirection_x_rays]  # initially false
        local s=s_out[idirection_x_rays]                  # gibberish the first iteration
        local point_x2,point_y2,direction_x2,direction_y2 # normal to the ellipse
        local i_wedge= retrieval_i[idirection_x_rays,iter] # gibberish the first iteration
        local j_wedge= retrieval_j[idirection_x_rays,iter] # gibberish the first iteration
        local i_wedge_plus_1 = i_wedge+1
        local j_wedge_plus_1 = j_wedge
        local s_bottom  = T(-99)
        local s_top   = T(-99)

        if 0<j_wedge<=Matm_n
          s_top = atm_h[j_wedge]
        end
        if 0<j_wedge_plus_1<=Matm_n
          s_bottom = atm_h[j_wedge_plus_1]
        end


        retrieval_point_x[idirection_x_rays,iter]=point_x1
        retrieval_point_y[idirection_x_rays,iter]=point_y1
        retrieval_direction_x[idirection_x_rays,iter]=direction_x1
        retrieval_direction_y[idirection_x_rays,iter]=direction_y1


        # this value is computed when the wedge is found when the ray tracing is all in one loop
        #local nₜ = nₜ[idirection_x_rays]                   # gibberish the first iteration
        begin
          if iter==1
            if  initialized==false




              # set the initial value of s to be the top of the atmosphere
              s = max_altitude  # assumed to be in a descending order
              nᵢ= free_space    # assumed to be vacuum
              θ = atan(point_y1,point_x1) # find the initial value of θ to simplify the computation of the first iteration
              θmin = -Inf
              θmax =  Inf
              s_bottom = max_altitude
              s_top    = max_altitude + 100*ϵ  # just to be sure it is higher and it exists
              retrieval_i[idirection_x_rays,1]=-2 # value of θ
              retrieval_j[idirection_x_rays,1]= 0 # value of h
              retrieval_n[idirection_x_rays,1]=T(1)
              retrieval_θ[idirection_x_rays,1]=NaN
              retrieval_h[idirection_x_rays,1]=NaN




            else

              i_wedge = retrieval_i[idirection_x_rays,1]
              j_wedge = retrieval_j[idirection_x_rays,1]

              s_bottom = atm_h[j_wedge_plus_1]
              s_top    = atm_h[j_wedge]
              s        = isAscending ? s_top : s_bottom
              nᵢ       = retrieval_n[idirection_x_rays,1]
              θ        = retrieval_θ[idirection_x_rays,1]

              θmin     = atm_θ[i_wedge]
              θmax     = atm_θ[i_wedge_plus_1]
            end
          end
        end

        @debug "---------------------------------------------"
        @debug " try to minimize towards $s "
        @debug " top $(s_top) and bottom$(s_bottom))"
        @debug " isAscending is $isAscending"
        @debug "---------------------------------------------"

        # update indexes
        begin
          i_wedge        = retrieval_i[idirection_x_rays,iter]
          i_wedge_plus_1 = i_wedge+1
          j_wedge        = retrieval_j[idirection_x_rays,iter]
          j_wedge_plus_1 = j_wedge+1
          if IsPeriodic
            i_wedge = mod1(i_wedge,Natm_n)
            i_wedge_plus_1 = mod1(i_wedge+1,Natm_n)
          end
          i_wedge = i_wedge< Natm_n ? i_wedge : -1
          i_wedge_plus_1 = i_wedge_plus_1< Natm_n ? i_wedge_plus_1 : -1
          j_wedge = j_wedge< Matm_n ? j_wedge : -1
          j_wedge_plus_1 = j_wedge_plus_1< Matm_n ? j_wedge_plus_1 : -1
        end
        s_top = atm_h[j_wedge]
        s_bottom = atm_h[j_wedge_plus_1]

        # stop the iteration if the ray has reached the top of the atmosphere
        # needs to be at the top because the bottom has the stopping condition to stop when all rays have number_rays_stopped
        # if not
        #   1. The code would do some unnecessary computation
        #   2. The code would update incorrectly the number_rays_stopped, leading to an early stop
        if ((iter>1) &&  (j_wedge<1 ||  i_wedge<1 || i_wedge_plus_1<1 || j_wedge_plus_1<1))

          retrieval_h[idirection_x_rays,iter+1]=NaN
          retrieval_t[idirection_x_rays,iter+1]=NaN
          retrieval_point_x[idirection_x_rays,iter+1]=NaN
          retrieval_point_y[idirection_x_rays,iter+1]=NaN
          retrieval_direction_x[idirection_x_rays,iter+1]=NaN
          retrieval_direction_y[idirection_x_rays,iter+1]=NaN
          retrieval_i[idirection_x_rays,iter+1]=i_wedge
          retrieval_j[idirection_x_rays,iter+1]=j_wedge
          retrieval_n[idirection_x_rays,iter+1]=nᵢ
          retrieval_θ[idirection_x_rays,iter+1]=NaN
          continue
        end

        # begin the intersection loop

        begin
          # normalization of the direction
          local norm = hypot(direction_x1, direction_y1)
          begin
            direction_x1 /= norm
            direction_y1 /= norm
          end

          local origin_times_direction = point_x1 * direction_x1 + point_y1 * direction_y1
          local cosθ, sinθ, bcosθ, bsinθ, cosθ², sinθ², half_sin2θ, cos2θ
          local R, N, Fx, Fy, t, Px, Py, fx, fy, f, fold
          local dpoint_x2dθ, dpoint_y2dθ, dpoint_x2dθ_0, dpoint_y2dθ_0,  ddirection_x2dθ, ddirection_y2dθ
          local dFxdθ, dFydθ, dtdθ, dPxdθ,  dPydθ,  dfxdθ,  dfydθ
          local g,h,p
          local d²point_x2dθ²,d²point_y2dθ², d²direction_x2dθ²_0
          local d²direction_y2dθ²_0, d²direction_x2dθ², d²direction_y2dθ², d²Fxdθ², d²Fydθ²
          local d²tdθ², d²Pxdθ², d²Pydθ²
          local dfxdθ_squared, dfydθ_squared, d²fxdθ², d²fydθ²
          local penality_t = T(0)
          local f_with_penality
          local rho = 0.0 #penality factor
          local gamma = 1.5 #penality factor



          begin
            # trigonometric function needed for computing the function

            begin
              cosθ = cos(θ)
              sinθ = sin(θ)
              bcosθ = b_normalized * cosθ
              bsinθ = b_normalized * sinθ
              cosθ² = cosθ * cosθ
              #
            end


            # computation point on the ellipse at distance s
            begin
              R = 1 - e² * cosθ²
              N = 1 / sqrt(R)
              point_x2 = cosθ
              point_y2 = bsinθ
              direction_x2 = bcosθ
              direction_y2 = sinθ
              # point on the normal of the ellipse

              Fx = point_x2 + s * direction_x2 * N
              Fy = point_y2 + s * direction_y2 * N
            end
            # compute t minimum  for θ and the relative distance squared f
            # dropped 1/2 so that the square root is the displacement
            t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
            # point on the ray
            begin
              Px = point_x1 + t * direction_x1
              Py = point_y1 + t * direction_y1
            end
            begin
              fx = Fx - Px
              fy = Fy - Py
              f = fx * fx + fy * fy
            end

            # if the distance is alreadirection_y small enough or the function is nan
            # it means I am either at the end point or that the function has reached the minimum

            # set f-> fold
            fold = f
            @debug "f: $f"
          end
          # Newton loop

          for k in 1:kloops

            ############################
            # Updates every 5 iterations
            if penality_t>0
              rho=min(rho*gamma,rho_max)
            end
            ############################

            # unroll block of 5 iterations
            # unroll block
            #   new function value and evaluation breaking criteria
            #   step block 1
            begin
              # update k ← k + 1
              # compute new θ
              begin
                # trigonometric function needed for computing the gradient and hessian
                begin
                  sinθ² = sinθ * sinθ
                  half_sin2θ = sinθ * cosθ
                  cos2θ = cosθ² - sinθ²
                end
                # compute the gradient
                begin
                  begin
                    dpoint_x2dθ = -sinθ
                    dpoint_y2dθ = bcosθ
                    dpoint_x2dθ_0 = -bsinθ
                    dpoint_y2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddirection_x2dθ = (dpoint_x2dθ_0 + N² * dNdθ * direction_x2) * N
                    ddirection_y2dθ = (dpoint_y2dθ_0 + N² * dNdθ * direction_y2) * N
                  end
                  begin
                    dFxdθ = dpoint_x2dθ + s * ddirection_x2dθ
                    dFydθ = dpoint_y2dθ + s * ddirection_y2dθ
                    dtdθ = direction_x1 * dFxdθ + direction_y1 * dFydθ
                    dPxdθ = direction_x1 * dtdθ
                    dPydθ = direction_y1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
                  g = fx * dfxdθ + fy * dfydθ
                end
                # compute the hessian
                begin
                  begin
                    d²point_x2dθ² = -point_x2
                    d²point_y2dθ² = -point_y2
                    d²direction_x2dθ²_0 = -direction_x2
                    d²direction_y2dθ²_0 = -direction_y2
                  end
                  begin
                    d²direction_x2dθ² = (d²direction_x2dθ²_0 + (2 * dpoint_x2dθ_0 * dNdθ + d²Ndθ² * N² * direction_x2) * N²) * N
                    d²direction_y2dθ² = (d²direction_y2dθ²_0 + (2 * dpoint_y2dθ_0 * dNdθ + d²Ndθ² * N² * direction_y2) * N²) * N
                    d²Fxdθ² = d²point_x2dθ² + d²direction_x2dθ²
                    d²Fydθ² = d²point_y2dθ² + d²direction_y2dθ²
                  end
                  begin
                    d²tdθ² = direction_x1 * d²Fxdθ² + direction_y1 * d²Fydθ²
                    d²Pxdθ² = direction_x1 * d²tdθ²
                    d²Pydθ² = direction_y1 * d²tdθ²
                    dfxdθ_squared = dfxdθ * dfxdθ
                    dfydθ_squared = dfydθ * dfydθ
                    d²fxdθ² = d²Fxdθ² - d²Pxdθ²
                    d²fydθ² = d²Fydθ² - d²Pydθ²
                  end
                  h = dfxdθ_squared + dfydθ_squared + fx * d²fxdθ² + fy * d²fydθ²
                  # insure positive definess of the hessian by adding a const
                  # similar to how LDLT works to ensure positive definiteness of matrix
                  h = abs(h) > 10 ^ -5 ? abs(h) :  10 ^ -5
                end
                # penality for negative t
                begin
                  penality_t = 0
                  g_t = 0
                  h_t = 0
                  if t < 0
                    t² = t * t
                    penality_t = -rho * t*t²
                    g_t= -3t²*dtdθ
                    h_t = -t²*d²tdθ²-6t*dtdθ*dtdθ
                    g_t*=rho
                    h_t*=rho
                  end
                  f_with_penality = f + penality_t

                  g += g_t
                  h += h_t



                end
                # compute the newton step
                p = -g / h
                # compute next θ
                θ = mod2pi(θ + p)
              end
              # compute new (f,t)
              begin
                begin
                  cosθ = cos(θ)
                  sinθ = sin(θ)
                  bcosθ = b_normalized * cosθ
                  bsinθ = b_normalized * sinθ
                  cosθ² = cosθ * cosθ
                end
                begin
                  R = 1 - e² * cosθ²
                  N = 1 / sqrt(R)
                  point_x2 = cosθ
                  point_y2 = bsinθ
                  direction_x2 = bcosθ
                  direction_y2 = sinθ
                  Fx = point_x2 + s * direction_x2 * N
                  Fy = point_y2 + s * direction_y2 * N
                end
                t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
                begin
                  Px = point_x1 + t * direction_x1
                  Py = point_y1 + t * direction_y1
                end
                begin
                  fx = Fx - Px
                  fy = Fy - Py
                  f = fx * fx + fy * fy
                end
              end
              # stopping criteria and update fold ← f
              begin



                if abs(f-fold) < δ && (penality_t==0)

                    break
                end

                fold = f
              end
            end

            #   step block 2
            begin
              # update k ← k + 1
              # compute new θ
              begin
                # trigonometric function needed for computing the gradient and hessian
                begin
                  sinθ² = sinθ * sinθ
                  half_sin2θ = sinθ * cosθ
                  cos2θ = cosθ² - sinθ²
                end
                # compute the gradient
                begin
                  begin
                    dpoint_x2dθ = -sinθ
                    dpoint_y2dθ = bcosθ
                    dpoint_x2dθ_0 = -bsinθ
                    dpoint_y2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddirection_x2dθ = (dpoint_x2dθ_0 + N² * dNdθ * direction_x2) * N
                    ddirection_y2dθ = (dpoint_y2dθ_0 + N² * dNdθ * direction_y2) * N
                  end
                  begin
                    dFxdθ = dpoint_x2dθ + s * ddirection_x2dθ
                    dFydθ = dpoint_y2dθ + s * ddirection_y2dθ
                    dtdθ = direction_x1 * dFxdθ + direction_y1 * dFydθ
                    dPxdθ = direction_x1 * dtdθ
                    dPydθ = direction_y1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
                  g = fx * dfxdθ + fy * dfydθ
                end
                # compute the hessian
                begin
                  begin
                    d²point_x2dθ² = -point_x2
                    d²point_y2dθ² = -point_y2
                    d²direction_x2dθ²_0 = -direction_x2
                    d²direction_y2dθ²_0 = -direction_y2
                  end
                  begin
                    d²direction_x2dθ² = (d²direction_x2dθ²_0 + (2 * dpoint_x2dθ_0 * dNdθ + d²Ndθ² * N² * direction_x2) * N²) * N
                    d²direction_y2dθ² = (d²direction_y2dθ²_0 + (2 * dpoint_y2dθ_0 * dNdθ + d²Ndθ² * N² * direction_y2) * N²) * N
                    d²Fxdθ² = d²point_x2dθ² + d²direction_x2dθ²
                    d²Fydθ² = d²point_y2dθ² + d²direction_y2dθ²
                  end
                  begin
                    d²tdθ² = direction_x1 * d²Fxdθ² + direction_y1 * d²Fydθ²
                    d²Pxdθ² = direction_x1 * d²tdθ²
                    d²Pydθ² = direction_y1 * d²tdθ²
                    dfxdθ_squared = dfxdθ * dfxdθ
                    dfydθ_squared = dfydθ * dfydθ
                    d²fxdθ² = d²Fxdθ² - d²Pxdθ²
                    d²fydθ² = d²Fydθ² - d²Pydθ²
                  end
                  h = dfxdθ_squared + dfydθ_squared + fx * d²fxdθ² + fy * d²fydθ²
                  # insure positive definess of the hessian by adding a const
                  # similar to how LDLT works to ensure positive definiteness of matrix
                  h = abs(h) > 10 ^ -5 ? abs(h) :  10 ^ -5
                end
                # penality for negative t
                begin
                  penality_t = 0
                  g_t = 0
                  h_t = 0
                  if t < 0
                    t² = t * t
                    penality_t = -rho * t*t²
                    g_t= -3t²*dtdθ
                    h_t = -t²*d²tdθ²-6t*dtdθ*dtdθ
                    g_t*=rho
                    h_t*=rho
                  end
                  f_with_penality = f + penality_t

                  g += g_t
                  h += h_t
                end
                # compute the newton step
                p = -g / h
                # compute next θ
                θ = mod2pi(θ + p)
              end
              # compute new (f,t)
              begin
                begin
                  cosθ = cos(θ)
                  sinθ = sin(θ)
                  bcosθ = b_normalized * cosθ
                  bsinθ = b_normalized * sinθ
                  cosθ² = cosθ * cosθ
                end
                begin
                  R = 1 - e² * cosθ²
                  N = 1 / sqrt(R)
                  point_x2 = cosθ
                  point_y2 = bsinθ
                  direction_x2 = bcosθ
                  direction_y2 = sinθ
                  Fx = point_x2 + s * direction_x2 * N
                  Fy = point_y2 + s * direction_y2 * N
                end
                t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
                begin
                  Px = point_x1 + t * direction_x1
                  Py = point_y1 + t * direction_y1
                end
                begin
                  fx = Fx - Px
                  fy = Fy - Py
                  f = fx * fx + fy * fy
                end
              end
              # stopping criteria and update fold ← f
              begin


                if abs(f-fold) < δ && (penality_t==0)

                    break
                end

                fold = f
              end
            end

            #   step block 3
            begin
              # update k ← k + 1
              # compute new θ
              begin
                # trigonometric function needed for computing the gradient and hessian
                begin
                  sinθ² = sinθ * sinθ
                  half_sin2θ = sinθ * cosθ
                  cos2θ = cosθ² - sinθ²
                end
                # compute the gradient
                begin
                  begin
                    dpoint_x2dθ = -sinθ
                    dpoint_y2dθ = bcosθ
                    dpoint_x2dθ_0 = -bsinθ
                    dpoint_y2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddirection_x2dθ = (dpoint_x2dθ_0 + N² * dNdθ * direction_x2) * N
                    ddirection_y2dθ = (dpoint_y2dθ_0 + N² * dNdθ * direction_y2) * N
                  end
                  begin
                    dFxdθ = dpoint_x2dθ + s * ddirection_x2dθ
                    dFydθ = dpoint_y2dθ + s * ddirection_y2dθ
                    dtdθ = direction_x1 * dFxdθ + direction_y1 * dFydθ
                    dPxdθ = direction_x1 * dtdθ
                    dPydθ = direction_y1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
                  g = fx * dfxdθ + fy * dfydθ
                end
                # compute the hessian
                begin
                  begin
                    d²point_x2dθ² = -point_x2
                    d²point_y2dθ² = -point_y2
                    d²direction_x2dθ²_0 = -direction_x2
                    d²direction_y2dθ²_0 = -direction_y2
                  end
                  begin
                    d²direction_x2dθ² = (d²direction_x2dθ²_0 + (2 * dpoint_x2dθ_0 * dNdθ + d²Ndθ² * N² * direction_x2) * N²) * N
                    d²direction_y2dθ² = (d²direction_y2dθ²_0 + (2 * dpoint_y2dθ_0 * dNdθ + d²Ndθ² * N² * direction_y2) * N²) * N
                    d²Fxdθ² = d²point_x2dθ² + d²direction_x2dθ²
                    d²Fydθ² = d²point_y2dθ² + d²direction_y2dθ²
                  end
                  begin
                    d²tdθ² = direction_x1 * d²Fxdθ² + direction_y1 * d²Fydθ²
                    d²Pxdθ² = direction_x1 * d²tdθ²
                    d²Pydθ² = direction_y1 * d²tdθ²
                    dfxdθ_squared = dfxdθ * dfxdθ
                    dfydθ_squared = dfydθ * dfydθ
                    d²fxdθ² = d²Fxdθ² - d²Pxdθ²
                    d²fydθ² = d²Fydθ² - d²Pydθ²
                  end
                  h = dfxdθ_squared + dfydθ_squared + fx * d²fxdθ² + fy * d²fydθ²
                  # insure positive definess of the hessian by adding a const
                  # similar to how LDLT works to ensure positive definiteness of matrix
                  h = abs(h) > 10 ^ -5 ? abs(h) :  10 ^ -5
                end
                # penality for negative t
                begin
                  penality_t = 0
                  g_t = 0
                  h_t = 0
                  if t < 0
                    t² = t * t
                    penality_t = -rho * t*t²
                    g_t= -3t²*dtdθ
                    h_t = -t²*d²tdθ²-6t*dtdθ*dtdθ
                    g_t*=rho
                    h_t*=rho
                  end
                  f_with_penality = f + penality_t

                  g += g_t
                  h += h_t
                end
                # compute the newton step
                p = -g / h
                # compute next θ
                θ = mod2pi(θ + p)
              end
              # compute new (f,t)
              begin
                begin
                  cosθ = cos(θ)
                  sinθ = sin(θ)
                  bcosθ = b_normalized * cosθ
                  bsinθ = b_normalized * sinθ
                  cosθ² = cosθ * cosθ
                end
                begin
                  R = 1 - e² * cosθ²
                  N = 1 / sqrt(R)
                  point_x2 = cosθ
                  point_y2 = bsinθ
                  direction_x2 = bcosθ
                  direction_y2 = sinθ
                  Fx = point_x2 + s * direction_x2 * N
                  Fy = point_y2 + s * direction_y2 * N
                end
                t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
                begin
                  Px = point_x1 + t * direction_x1
                  Py = point_y1 + t * direction_y1
                end
                begin
                  fx = Fx - Px
                  fy = Fy - Py
                  f = fx * fx + fy * fy
                end
              end
              # stopping criteria and update fold ← f
              begin


                if abs(f-fold) < δ && (penality_t==0)

                    break
                end

                fold = f
              end
            end

            #   step block 4
            begin
              # update k ← k + 1
              # compute new θ
              begin
                # trigonometric function needed for computing the gradient and hessian
                begin
                  sinθ² = sinθ * sinθ
                  half_sin2θ = sinθ * cosθ
                  cos2θ = cosθ² - sinθ²
                end
                # compute the gradient
                begin
                  begin
                    dpoint_x2dθ = -sinθ
                    dpoint_y2dθ = bcosθ
                    dpoint_x2dθ_0 = -bsinθ
                    dpoint_y2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddirection_x2dθ = (dpoint_x2dθ_0 + N² * dNdθ * direction_x2) * N
                    ddirection_y2dθ = (dpoint_y2dθ_0 + N² * dNdθ * direction_y2) * N
                  end
                  begin
                    dFxdθ = dpoint_x2dθ + s * ddirection_x2dθ
                    dFydθ = dpoint_y2dθ + s * ddirection_y2dθ
                    dtdθ = direction_x1 * dFxdθ + direction_y1 * dFydθ
                    dPxdθ = direction_x1 * dtdθ
                    dPydθ = direction_y1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
                  g = fx * dfxdθ + fy * dfydθ
                end
                # compute the hessian
                begin
                  begin
                    d²point_x2dθ² = -point_x2
                    d²point_y2dθ² = -point_y2
                    d²direction_x2dθ²_0 = -direction_x2
                    d²direction_y2dθ²_0 = -direction_y2
                  end
                  begin
                    d²direction_x2dθ² = (d²direction_x2dθ²_0 + (2 * dpoint_x2dθ_0 * dNdθ + d²Ndθ² * N² * direction_x2) * N²) * N
                    d²direction_y2dθ² = (d²direction_y2dθ²_0 + (2 * dpoint_y2dθ_0 * dNdθ + d²Ndθ² * N² * direction_y2) * N²) * N
                    d²Fxdθ² = d²point_x2dθ² + d²direction_x2dθ²
                    d²Fydθ² = d²point_y2dθ² + d²direction_y2dθ²
                  end
                  begin
                    d²tdθ² = direction_x1 * d²Fxdθ² + direction_y1 * d²Fydθ²
                    d²Pxdθ² = direction_x1 * d²tdθ²
                    d²Pydθ² = direction_y1 * d²tdθ²
                    dfxdθ_squared = dfxdθ * dfxdθ
                    dfydθ_squared = dfydθ * dfydθ
                    d²fxdθ² = d²Fxdθ² - d²Pxdθ²
                    d²fydθ² = d²Fydθ² - d²Pydθ²
                  end
                  h = dfxdθ_squared + dfydθ_squared + fx * d²fxdθ² + fy * d²fydθ²
                  # insure positive definess of the hessian by adding a const
                  # similar to how LDLT works to ensure positive definiteness of matrix
                  h = abs(h) > 10 ^ -5 ? abs(h) :  10 ^ -5
                end
                # penality for negative t
                begin
                  penality_t = 0
                  g_t = 0
                  h_t = 0
                  if t < 0
                    t² = t * t
                    penality_t = -rho * t*t²
                    g_t= -3t²*dtdθ
                    h_t = -t²*d²tdθ²-6t*dtdθ*dtdθ
                    g_t*=rho
                    h_t*=rho
                  end
                  f_with_penality = f + penality_t

                  g += g_t
                  h += h_t
                end
                # compute the newton step
                p = -g / h
                # compute next θ
                θ = mod2pi(θ + p)
              end
              # compute new (f,t)
              begin
                begin
                  cosθ = cos(θ)
                  sinθ = sin(θ)
                  bcosθ = b_normalized * cosθ
                  bsinθ = b_normalized * sinθ
                  cosθ² = cosθ * cosθ
                end
                begin
                  R = 1 - e² * cosθ²
                  N = 1 / sqrt(R)
                  point_x2 = cosθ
                  point_y2 = bsinθ
                  direction_x2 = bcosθ
                  direction_y2 = sinθ
                  Fx = point_x2 + s * direction_x2 * N
                  Fy = point_y2 + s * direction_y2 * N
                end
                t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
                begin
                  Px = point_x1 + t * direction_x1
                  Py = point_y1 + t * direction_y1
                end
                begin
                  fx = Fx - Px
                  fy = Fy - Py
                  f = fx * fx + fy * fy
                end
              end
              # stopping criteria and update fold ← f
              begin


                if abs(f-fold) < δ && (penality_t==0)

                    break
                end

                fold = f
              end
            end
            #   step block 5

            begin
              # update k ← k + 1
              # compute new θ
              begin
                # trigonometric function needed for computing the gradient and hessian
                begin
                  sinθ² = sinθ * sinθ
                  half_sin2θ = sinθ * cosθ
                  cos2θ = cosθ² - sinθ²
                end
                # compute the gradient
                begin
                  begin
                    dpoint_x2dθ = -sinθ
                    dpoint_y2dθ = bcosθ
                    dpoint_x2dθ_0 = -bsinθ
                    dpoint_y2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddirection_x2dθ = (dpoint_x2dθ_0 + N² * dNdθ * direction_x2) * N
                    ddirection_y2dθ = (dpoint_y2dθ_0 + N² * dNdθ * direction_y2) * N
                  end
                  begin
                    dFxdθ = dpoint_x2dθ + s * ddirection_x2dθ
                    dFydθ = dpoint_y2dθ + s * ddirection_y2dθ
                    dtdθ = direction_x1 * dFxdθ + direction_y1 * dFydθ
                    dPxdθ = direction_x1 * dtdθ
                    dPydθ = direction_y1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
                  g = fx * dfxdθ + fy * dfydθ
                end
                # compute the hessian
                begin
                  begin
                    d²point_x2dθ² = -point_x2
                    d²point_y2dθ² = -point_y2
                    d²direction_x2dθ²_0 = -direction_x2
                    d²direction_y2dθ²_0 = -direction_y2
                  end
                  begin
                    d²direction_x2dθ² = (d²direction_x2dθ²_0 + (2 * dpoint_x2dθ_0 * dNdθ + d²Ndθ² * N² * direction_x2) * N²) * N
                    d²direction_y2dθ² = (d²direction_y2dθ²_0 + (2 * dpoint_y2dθ_0 * dNdθ + d²Ndθ² * N² * direction_y2) * N²) * N
                    d²Fxdθ² = d²point_x2dθ² + d²direction_x2dθ²
                    d²Fydθ² = d²point_y2dθ² + d²direction_y2dθ²
                  end
                  begin
                    d²tdθ² = direction_x1 * d²Fxdθ² + direction_y1 * d²Fydθ²
                    d²Pxdθ² = direction_x1 * d²tdθ²
                    d²Pydθ² = direction_y1 * d²tdθ²
                    dfxdθ_squared = dfxdθ * dfxdθ
                    dfydθ_squared = dfydθ * dfydθ
                    d²fxdθ² = d²Fxdθ² - d²Pxdθ²
                    d²fydθ² = d²Fydθ² - d²Pydθ²
                  end
                  h = dfxdθ_squared + dfydθ_squared + fx * d²fxdθ² + fy * d²fydθ²
                  # insure positive definess of the hessian by adding a const
                  # similar to how LDLT works to ensure positive definiteness of matrix
                  h = abs(h) > 10 ^ -5 ? abs(h) :  10 ^ -5
                end
                # penality for negative t
                begin
                  penality_t = 0
                  g_t = 0
                  h_t = 0
                  if t < 0
                    t² = t * t
                    penality_t = -rho * t*t²
                    g_t= -3t²*dtdθ
                    h_t = -t²*d²tdθ²-6t*dtdθ*dtdθ
                    g_t*=rho
                    h_t*=rho
                  end
                  f_with_penality = f + penality_t

                  g += g_t
                  h += h_t
                end
                # compute the newton step
                p = -g / h
                # compute next θ
                θ = mod2pi(θ + p)
              end
              # compute new (f,t)
              begin
                begin
                  cosθ = cos(θ)
                  sinθ = sin(θ)
                  bcosθ = b_normalized * cosθ
                  bsinθ = b_normalized * sinθ
                  cosθ² = cosθ * cosθ
                end
                begin
                  R = 1 - e² * cosθ²
                  N = 1 / sqrt(R)
                  point_x2 = cosθ
                  point_y2 = bsinθ
                  direction_x2 = bcosθ
                  direction_y2 = sinθ
                  Fx = point_x2 + s * direction_x2 * N
                  Fy = point_y2 + s * direction_y2 * N
                end
                t = -origin_times_direction + (direction_x1 * Fx + direction_y1 * Fy)
                begin
                  Px = point_x1 + t * direction_x1
                  Py = point_y1 + t * direction_y1
                end
                begin
                  fx = Fx - Px
                  fy = Fy - Py
                  f = fx * fx + fy * fy
                end
              end
              # stopping criteria and update fold ← f
              begin


                if abs(f-fold) < δ && (penality_t==0)

                    break
                end

                fold = f
              end
            end

          end
          # clamping angle to the range
          # if θmin <= θmax  thene θ ∈ [θmin, θmax]


          begin
            if (θmin <= θmax && θmin < θ < θmax) || (θmin > θmax && (θ > θmin || θ < θmax))


                s += sqrt(f) * directional_sign

            else



              # angle clamping f(θ,θmin,θmax)
              let
                if θmin<=θmax
                  θ=clamp(θ,θmin,θmax)

                else

                  dmin = mod(θ - θmin, 2π)
                  dmax = mod(θmax - θ, 2π)
                  θ = dmin < dmax ? θmin : θmax
                end
                θ
              end
              #update (point_x1,point_x2,θ,t,s)
              begin
                point_x2 = cos(θ)
                point_y2 = sin(θ) * b_normalized
                direction_x2 = b_normalized * cos(θ)
                direction_y2 = sin(θ)
                normF = hypot(direction_x2, direction_y2)
                direction_x2 /= normF
                direction_y2 /= normF

                det = direction_x1 * direction_y2 - direction_y1 * direction_x2
                ϵ = 1.0e-10        # TO DO: maybe add this value as a kwargs
                s = NaN        # initial s to NaN
                t = NaN        # initial t to NaN
                local Δp12x = point_x2 - point_x1
                local Δp12y = point_y2 - point_y1


                if abs(det) > ϵ
                  t = (Δp12x * direction_y2 - Δp12y * direction_x2) / det
                  s = (Δp12x * direction_y1 - Δp12y * direction_x1) / det
                end
                if (s < 0) # it starts from the ellipse surface, so a negative value is not possible
                  s = NaN
                end
              end
            end
          end
          # check if the direction of the ray has changed and if it did, set the ascending to true
          # this helps to find easily the tangent quote
          # next iteration will look for the ascending s instead of the descending one

          #update all the output for the current iteration
        end

        # find new index_i and index_j as well as the neighbors refractive index

        begin
          if iter==1 && initialized==false
            # find the index of the wedge using binary search
            # assume left orientation of the ray
            # TO DO: add the right orientation of the ray
            #@debug "t: $t, s: $s f: $f"
            #return (atm_θ,θ)  #debug
            i_wedge = findlast(atm_θ.<= θ)
            i_wedge_plus_1 = i_wedge+1
            j_wedge = 1
            j_wedge_plus_1 = 2
            #if s>max_altitude+ϵ  # the first ray never intersected the atmosphere
            #  j_wedge = -2
            #end
            # NOTE: the first intersection can only be downward so direction_x2 and direction_y2 are alreadirection_y the outward normal
            # to the wedge


          else
            @debug "--------------------------------------------------"
            @debug "  θ   : $(rad2deg(θ))°  "
            @debug "  θmin: $(rad2deg(θmin))°  "
            @debug "  θmax: $(rad2deg(θmax))°  "
            @debug "  Δh  : $(f)  "
            @debug "  s_top: $(s_top)  "
            @debug "  s_bottom: $(s_bottom)  "
            @debug "  s  : $(s)  "
            @debug "--------------------------------------------------"
            begin
              local tmp = direction_x2
              if  θmin<θ<θmax
                @debug "Between 2 wedges"
                if  (abs(f)>10^-5 && isAscending==false)
                    @debug " in the middle of the atmosphere"
                    tangent_quote[idirection_x_rays]=s
                    isAscending = true
                    # assume left handiness
                    # TO DO: implement the right handiness
                    # set to the next ray
                    θ=θmax
                    local t_old=t
                    local i_wedge_old=i_wedge
                    local j_wedge_old=j_wedge
                    # update to the next wedge
                    begin
                      i_wedge = i_wedge+1
                      i_wedge_plus_1 = i_wedge+1
                      direction_x2 =  direction_y2
                      direction_y2 = -tmp
                    end
                    # the ray is going up so I send it to the next ray after registering the
                    # tangent quote
                    #update (point_x1,point_x2,θ,t,s)
                    begin
                      point_x2 = cos(θ)
                      point_y2 = sin(θ) * b_normalized
                      direction_x2 = b_normalized * cos(θ)
                      direction_y2 = sin(θ)
                      normF = hypot(direction_x2, direction_y2)
                      direction_x2 /= normF
                      direction_y2 /= normF

                      det = direction_x1 * direction_y2 - direction_y1 * direction_x2
                      ϵ = 1.0e-10        # TO DO: maybe add this value as a kwargs
                      s = NaN        # initial s to NaN
                      t = NaN        # initial t to NaN
                      local Δp12x = point_x2 - point_x1
                      local Δp12y = point_y2 - point_y1


                      if abs(det) > ϵ
                        t = (Δp12x * direction_y2 - Δp12y * direction_x2) / det
                        s = (Δp12x * direction_y1 - Δp12y * direction_x1) / det
                      end
                      if (s < 0) # it starts from the ellipse surface, so a negative value is not possible
                        s = NaN
                      end
                    end


                elseif isAscending==true
                  @debug " going up"
                  j_wedge = j_wedge-1  # the direction of h is descending
                  direction_x2 = -direction_x2
                  direction_y2 = -direction_y2
                else
                  @debug " going down"
                  j_wedge = j_wedge+1  # the direction of h is ascending
                end
              elseif θ==θmax
                @debug "Touching left"
                begin
                  i_wedge = i_wedge+1
                  i_wedge_plus_1 = i_wedge+1
                  direction_x2 =  direction_y2
                  direction_y2 = -tmp
                end
              elseif θ==θmin
                @debug "Touching right"
                begin
                  i_wedge = i_wedge-1
                  direction_x2 = -direction_y2
                  direction_y2 =  tmp
                end
              # something went wrong
              else
                j_wedge = -2
              end
              # check i_index for periodic radial distribution
            end

            # update indexes
            begin
              i_wedge_plus_1 = i_wedge+1
              j_wedge_plus_1 = j_wedge+1
              if IsPeriodic
                i_wedge = mod1(i_wedge,Natm_n)
                i_wedge_plus_1 = mod1(i_wedge+1,Natm_n)
              end
              i_wedge = i_wedge< Natm_n ? i_wedge : -1
              i_wedge_plus_1 = i_wedge_plus_1< Natm_n ? i_wedge_plus_1 : -1
              j_wedge = j_wedge< Matm_n ? j_wedge : -1
              j_wedge_plus_1 = j_wedge_plus_1< Matm_n ? j_wedge_plus_1 : -1
            end





          end

          retrieval_i[idirection_x_rays,iter+1]=i_wedge
          retrieval_j[idirection_x_rays,iter+1]=j_wedge
        end


        # check if atmosphere has been reached or if the ray has failed to intersect
        begin
          # j_wedge is used as a flag in the code
          # -2 means an error
          # 0   means it left the atmosphere
          # -1  means it reached the ground
          if j_wedge<1 || i_wedge<1 || i_wedge_plus_1<1

            number_rays_stopped+=1

            continue
          end
        end
        # update the refractive index of the neighbor
        begin
          nₜ = atm_n[i_wedge,j_wedge]
          θmin = atm_θ[i_wedge]
          θmax = atm_θ[i_wedge_plus_1]
          s_top = atm_h[j_wedge]
          s_bottom = atm_h[j_wedge_plus_1]
        end
        ##################################################
        # DEBUG
        ##################################################
        @debug "--------------------------------------------------"
        @debug "Ray $(idirection_x_rays) at $(iter) iteration"
        @debug "--------------------------------------------------"
        @debug "i_wedge=$(i_wedge) j_wedge=$(j_wedge)"
        @debug "i_wedge_plus_1=$(i_wedge_plus_1) j_wedge_plus_1=$(j_wedge_plus_1)"
        @debug "n_incident=$(nᵢ) n_trasmitted=$(nₜ) n_01= $(n_incident/n_transmitted)"
        @debug "θmin=$(θmin) θmax=$(θmax)"
        @debug "s_top=$(s_top) s_bottom=$(s_bottom)"
        @debug "point_x1=$(point_x1) point_y1=$(point_y1)"
        @debug "point_x2=$(point_x2) point_y2=$(point_y2)"
        @debug "direction_x1=$(direction_x1) direction_y1=$(direction_y1)"
        @debug "direction_x2=$(direction_x2) direction_y2=$(direction_y2)"
        @debug "t=$(t) s=$(s)"
        @debug "θ=$(θ)"
        @debug "isAscenging? $(isAscending)"
        @debug "--------------------------------------------------"
        #return
        ##################################################
        # bending the ray
        if !(nᵢ==nₜ) # do bend only if the refractive index are different
          @debug " --------------------------------------------------"
          @debug "Bending the ray"
          @debug " --------------------------------------------------"
          let
            # check if it is intersecting a level or a ray
            # both directions are alreadirection_y normalized
            local n_incident = nᵢ
            local n_transmitted = nₜ
            local direcion_ray_x=direction_x1
            local direction_ray_y=direction_y1
            local Nx=direction_x2
            local Ny=direction_y2
            local n01=n_incident/n_transmitted
            local n01²=n01*n01
            local cosθ_incident=-(Nx*direcion_ray_x+Ny*direction_ray_y)
            local sinθ²_transmitted =n01²*(1-cosθ_incident*cosθ_incident)

            # check if the ray is internally reflected
            # this most likely happens if there is an issue with the atmosphere or if the tangent quote
            # happens to be at a level.
            if sinθ²_transmitted ≤ 1
              direction_x1= n01*direction_x1+(n01*cosθ_incident-sqrt(1-sinθ²_transmitted))*Nx
              direction_y1= n01*direction_y1+(n01*cosθ_incident-sqrt(1-sinθ²_transmitted))*Ny
            else
              direction_x1-=2*cosθ_incident*Nx
              direction_y1-=2*cosθ_incident*Ny
            end

          end
          @debug " direction_x_new=$(direction_x1) direction_y_new=$(direction_y1)"

          @debug " --------------------------------------------------"


        end
        # update the direction of the ray
        begin
          local norm_new_ray_direction = hypot(direction_x1, direction_y1)
          direction_x1 /= norm_new_ray_direction
          direction_y1 /= norm_new_ray_direction
          # first update the position using the array cause direction_x1,direction_y1 are alreadirection_y modified
          apoint_x[idirection_x_rays] = point_x1+t*adirection_x[idirection_x_rays]
          apoint_y[idirection_x_rays] = point_y1+t*adirection_y[idirection_x_rays]
          point_x1=apoint_x[idirection_x_rays]
          point_y1=apoint_y[idirection_x_rays]
          adirection_x[idirection_x_rays] = direction_x1
          adirection_y[idirection_x_rays] = direction_y1
          t_out[idirection_x_rays] = t
          θ_out[idirection_x_rays] = θ
          # set up next iteration
          s_out[idirection_x_rays] = isAscending ? s_top : s_bottom
          aθmin[idirection_x_rays] = θmin
          aθmax[idirection_x_rays] = θmax
          incident_refractive_index[idirection_x_rays]  = nₜ
          ascending[idirection_x_rays] = isAscending

        end

        # update output arrays
        begin
          # fill the retrieval output for the current iteration
          retrieval_i[idirection_x_rays,iter+1]=i_wedge
          retrieval_j[idirection_x_rays,iter+1]=j_wedge
          retrieval_θ[idirection_x_rays,iter+1]=θ
          retrieval_t[idirection_x_rays,iter+1]=t
          retrieval_h[idirection_x_rays,iter+1]=s

        end
    end
    # early stop condition
    if number_rays_stopped == NumRays
      break
    end
  end
end


@inline _rotation_matrix(θ)= [cosd(θ) sind(θ);-sind(θ) cosd(θ)]

function limb_angle(w,z,ang)
   θ = atan(z/w)
  (tx,ty)=(z,-w)|> x-> x./hypot(x...) .*-1.0
  #################################
  angle= ang*-1

  dir=_rotation_matrix(angle)*[tx,ty]
  return (dir[1],dir[2])
end

function nadir_angle(w,z,ang)
   θ = atan(z/w)
  (nx,ny)=(-w,-z)|> x-> x./hypot(x...) .*-1.0
  #################################
  angle= ang

  dir=_rotation_matrix(angle)*[nx,ny]
  return (dir[1],dir[2])
end

function nadir_angle_normal(nx,ny,ang;outward::Bool=true)
  inwardoutward = outward ? 1.0 : -1.0
  (nx,ny)=(nx,nx)|> x-> x./hypot(x...) .*inwardoutward
 #################################
 angle= ang

 dir=_rotation_matrix(angle)*[nx,ny]
 return (dir[1],dir[2])
end
