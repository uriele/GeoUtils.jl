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
  px::T
  py::T
  dx::T
  dy::T
  n ::T=1.0
  θmin::T=-Inf
  θmax::T=-Inf
  ascending::Bool=false
  InputRay{T}(px::T,py::T,dx::T,dy::T,n::T,θmin::T,θmax::T) where T= new{T}(px,py,dx,dy,n,θmin,θmax,false)
  InputRay(px::T,py::T,dx::T,dy::T,n::T,θmin::T,θmax::T) where T   = new{T}(px,py,dx,dy,n,θmin,θmax,false)
end

function Base.show(io::IO, p::InputRay{T}) where T
   print(io, "$((round(p.px;sigdigits=3),round(p.py;sigdigits=3))) ")
end

function Base.show(io::IO,  ::MIME"text/plain",p::InputRay{T}) where T
  px=round(p.px;sigdigits=3)
  py=round(p.py;sigdigits=3)
  dx=round(p.dx;sigdigits=3)
  dy=round(p.dy;sigdigits=3)
  n_minus_1 =round(p.n-1;sigdigits=3)
  θmin= p.θmin
  θmax= p.θmax

  if p.ascending
    print(io, "Ascending ")
  end
  println(io, "Ray{$(T)}")
  println(io, "  point      : $((px,py))")
  println(io, "  direction  : $((dx,dy))")
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
  px::T
  py::T
  dx::T
  dy::T
  ResultsRayTracing{T}(i::Int,j::Int,θ::T,t::T,h::T,n::T,px::T,py::T,dx::T,dy::T) where T= new{T}(i,j,θ,t,h,n,px,py,dx,dy)
  ResultsRayTracing(i::Int,j::Int,θ::T,t::T,h::T,n::T,px::T,py::T,dx::T,dy::T) where T   = new{T}(i,j,θ,t,h,n,px,py,dx,dy)
end

function Base.show(io::IO, p::ResultsRayTracing{T}) where T
  print(io, "$((round(p.px;sigdigits=3),round(p.py;sigdigits=3))) ")
end

function Base.show(io::IO,  ::MIME"text/plain",p::ResultsRayTracing{T}) where T
  px=round(p.px;sigdigits=3)
  py=round(p.py;sigdigits=3)
  dx=round(p.dx;sigdigits=3)
  dy=round(p.dy;sigdigits=3)
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
  println(io, "point     : $((px,py))")
  println(io, "direction : $((dx,dy))")
end


function initialize_theta(θin::T,px::T,py::T;ϵ=1e-2,kmax_init::Int=30)::T where T<:IEEEFloat
  @assert ϵ > 0 "δ has to be positive"
  @assert kmax_init > 0 "kmax_init has to be positive"

  local θ   = θin
  local Px  = px
  local Py  = py
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
    px2    = cosθ
    py2    = bsinθ

    Fx= px2
    Fy= py2

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
        dpx2dθ = -sinθ
        dpy2dθ = bcosθ

        d²px2dθ² = -px2
        d²py2dθ² = -py2

        dpdθ_squared= 1- e²*cosθ*cosθ

        g = fx*dpx2dθ+fy*dpy2dθ
        h = dpdθ_squared+fx*d²px2dθ²+fy*d²py2dθ²
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
        px2    = cosθ
        py2    = bsinθ

        Fx= px2
        Fy= py2

        begin
          fx=(Fx-Px)
          fy=(Fy-Py)
          f= fx*fx+fy*fy
        end
      end
      if abs(f-fold)<ϵ
        hasconverged = true
        #@info "Converged after $k iterations with Δf: $(abs(f-fold))"
        #return θ
      end
      fold=f
    end
  end
  #@info "Did not converged after $kmax_init iterations with Δf: $(abs(f-fold))"
  return θ

end

# Generated at: 2021-09-30T14:00:00.000
# High performance code obtained by interpolating the code from the quoted functions and the generated functions
# Manual inlining of functions and loops
# The code is optimized for the input parameters
 ##########################################################################################
"""
    initialize_theta!(θ_out::A,apx1::A,apy2::A;ϵ=1e-2,kmax_init::Int=30) where A<:AbstractArray{T} where T

  Initialize the initial angle theta as the approximate geocentric angle between the satellite and Earth

  # Arguments
    - `θ_out     <:Array{IEEEFloat}` : the optimization parameter θ (astromonical angle of the Earth)
    - `apx       <:Array{IEEEFloat}` : the x position of the ray
    - `apy       <:Array{IEEEFloat}` : the y position of the ray

  # Optional Arguments
    - `ϵ::IEEEFloat=1e-10` : the tolerance for the minimum distance
    - `kmax::Int=20`       : the maximum number of iterations for the Newton method

  See also: [`fast_ray_bending`](@ref) [`fast_initization_theta`](@ref) [`fast_ray_tracing`](@ref)
"""
function initialize_theta!(θ_out::A,apx1::A,apy2::A;ϵ=1e-2,kmax_init::Int=30) where A<:AbstractArray{T} where T
  @assert size(apx) == size(apy)   "position coordinates has to be the same size"
  @assert size(θ_out) == size(apx) "the optimization parameters θ has to be the same size"

  @assert δ > 0 "δ has to be positive"
  @assert kmax_init > 0 "kmax_init has to be positive"

  b_normalized = get_minoraxis(T)
  e² = get_e²(T)
  kinitloops = div(kmax, 5,RoundUp)
  @inbounds for i = eachindex(t_out)

    local fold=0.0
    local f   = 0.0
    local Px  = apx1[i]
    local Py  = apy1[i]
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
      px2    = cosθ
      py2    = bsinθ

      Fx= px2
      Fy= py2

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
          dpx2dθ = -sinθ
          dpy2dθ = bcosθ

          d²px2dθ² = -px2
          d²py2dθ² = -py2

          dpdθ_squared= 1- e²*cosθ*cosθ

          g = fx*dpx2dθ+fy*dpy2dθ
          h = dpdθ_squared+fx*d²px2dθ²+fy*d²py2dθ²
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
          px2    = cosθ
          py2    = bsinθ

          Fx= px2
          Fy= py2

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
          dpx2dθ = -sinθ
          dpy2dθ = bcosθ

          d²px2dθ² = -px2
          d²py2dθ² = -py2

          dpdθ_squared= 1- e²*cosθ*cosθ

          g = fx*dpx2dθ+fy*dpy2dθ
          h = dpdθ_squared+fx*d²px2dθ²+fy*d²py2dθ²
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
          px2    = cosθ
          py2    = bsinθ

          Fx= px2
          Fy= py2

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
          dpx2dθ = -sinθ
          dpy2dθ = bcosθ

          d²px2dθ² = -px2
          d²py2dθ² = -py2

          dpdθ_squared= 1- e²*cosθ*cosθ

          g = fx*dpx2dθ+fy*dpy2dθ
          h = dpdθ_squared+fx*d²px2dθ²+fy*d²py2dθ²
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
          px2    = cosθ
          py2    = bsinθ

          Fx= px2
          Fy= py2

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
          dpx2dθ = -sinθ
          dpy2dθ = bcosθ

          d²px2dθ² = -px2
          d²py2dθ² = -py2

          dpdθ_squared= 1- e²*cosθ*cosθ

          g = fx*dpx2dθ+fy*dpy2dθ
          h = dpdθ_squared+fx*d²px2dθ²+fy*d²py2dθ²
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
          px2    = cosθ
          py2    = bsinθ

          Fx= px2
          Fy= py2

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
          dpx2dθ = -sinθ
          dpy2dθ = bcosθ

          d²px2dθ² = -px2
          d²py2dθ² = -py2

          dpdθ_squared= 1- e²*cosθ*cosθ

          g = fx*dpx2dθ+fy*dpy2dθ
          h = dpdθ_squared+fx*d²px2dθ²+fy*d²py2dθ²
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
          px2    = cosθ
          py2    = bsinθ

          Fx= px2
          Fy= py2

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
    fast_minimization_distance!(t_out::A,θ_out::A,s_out::A,apx::A,apy::A,adx::A,ady::A,aθmin::A,aθmax::A,ascending::AB;δ=1e-10,kmax::Int=30, kwargs...) where {A<:AbstractArray{T},AB<:AbstractArray{Bool}} where T<:IEEEFloat

  Compute the minimum distance between a ray and an ellipse in 2D and modifies in place the optimization parameters `t_out`, `θ_out`, and `s_out` for each ray.
  It also updates the position apx, apy given the new t

  # Arguments
    - `t_out     <:Array{IEEEFloat}` : the optimization parameter t (length of the ray)
    - `θ_out     <:Array{IEEEFloat}` : the optimization parameter θ (astromonical angle of the Earth)
    - `s_out     <:Array{IEEEFloat}` : the optimization parameter s (tangent quote of the ray)
    - `apx       <:Array{IEEEFloat}` : the x position of the ray
    - `apy       <:Array{IEEEFloat}` : the y position of the ray
    - `adx       <:Array{IEEEFloat}` : the x direction of the ray
    - `ady       <:Array{IEEEFloat}` : the y direction of the ray
    - `aθmin     <:Array{IEEEFloat}` : the left angle of the wedge
    - `aθmax     <:Array{IEEEFloat}` : the right angle of the wedge
    - `ascending <:Array{Bool}`      : true defines the ascending part of the ray, finding the first true given the tangent quote

  # Optional Arguments
    - `δ::IEEEFloat=1e-10` : the tolerance for the minimum distance
    - `kmax::Int=30`       : the maximum number of iterations for the Newton method

  See also: [`fast_ray_bending`](@ref) [`fast_initization_theta`](@ref) [`fast_ray_tracing`](@ref)
"""
function fast_minimization_distance!(t_out::A,θ_out::A,s_out::A,apx::A,apy::A,adx::A,ady::A,aθmin::A,aθmax::A,ascending::AB;δ=1e-10,kmax::Int=30, kwargs...) where {A<:AbstractArray{T},AB<:AbstractArray{Bool}} where T<:IEEEFloat
  @assert size(adx) == size(ady) "direction coordinates has to be the same size"
  @assert size(apx) == size(apy) "position coordinates has to be the same size"
  @assert size(apx) == size(adx) "position and direction coordinates has to be the same size"
  @assert size(θ_out) == size(t_out) "the optimization parameters t and θ had to be the same size"
  @assert size(t_out) == size(adx) "optimization parameter and problem size should be the same"
  @assert size(t_out) == size(s_out) "optimization parameter and problem size should be the same"
  @assert δ > 0 "δ has to be positive"

  b_normalized = get_minoraxis(T)
  e² = get_e²(T)
  kloops = div(kmax, 5,RoundUp)
  @inbounds for i = eachindex(t_out)
    local px1 = apx[i]
    local py1 = apy[i]
    local dx1 = adx[i]
    local dy1 = ady[i]
    local θ = θ_out[i]
    local s = s_out[i]        # aims for a level s, returns the effective level found
    local directional_sign = ascending[i] ? T(-1.0) :  T(1.0)  # we correct the level altitude with the real one
    local θmin = aθmin[i]
    local θmax = aθmax[i]
    # normalization of the direction
    begin
      norm = hypot(dx1, dy1)
      dx1 /= norm
      dy1 /= norm
    end
    origin_times_direction = px1 * dx1 + py1 * dy1

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
        px2 = cosθ
        py2 = bsinθ
        dx2 = bcosθ
        dy2 = sinθ
        # point on the normal of the ellipse
        Fx = px2 + s * dx2 * N
        Fy = py2 + s * dy2 * N
      end
      # compute t minimum  for θ and the relative distance squared f
      # dropped 1/2 so that the square root is the displacement
      t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
      # point on the ray
      begin
        Px = px1 + t * dx1
        Py = py1 + t * dy1
      end
      begin
        fx = Fx - Px
        fy = Fy - Py
        f = fx * fx + fy * fy
      end

      # if the distance is already small enough or the function is nan
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
              dpx2dθ = -sinθ
              dpy2dθ = bcosθ
              dpx2dθ_0 = -bsinθ
              dpy2dθ_0 = cosθ
            end
            begin
              half_dR = e² * half_sin2θ
              half_dR² = half_dR * half_dR
              half_d²R = e² * cos2θ
              N² = N * N
              dNdθ = -half_dR
              d²Ndθ² = -half_d²R + 3 * N² * half_dR²
              ddx2dθ = (dpx2dθ_0 + N² * dNdθ * dx2) * N
              ddy2dθ = (dpy2dθ_0 + N² * dNdθ * dy2) * N
            end
            begin
              dFxdθ = dpx2dθ + s * ddx2dθ
              dFydθ = dpy2dθ + s * ddy2dθ
              dtdθ = dx1 * dFxdθ + dy1 * dFydθ
              dPxdθ = dx1 * dtdθ
              dPydθ = dy1 * dtdθ
              dfxdθ = dFxdθ - dPxdθ
              dfydθ = dFydθ - dPydθ
            end
            g = fx * dfxdθ + fy * dfydθ
          end
          # compute the hessian
          begin
            begin
              d²px2dθ² = -px2
              d²py2dθ² = -py2
              d²dx2dθ²_0 = -dx2
              d²dy2dθ²_0 = -dy2
            end
            begin
              d²dx2dθ² = (d²dx2dθ²_0 + (2 * dpx2dθ_0 * dNdθ + d²Ndθ² * N² * dx2) * N²) * N
              d²dy2dθ² = (d²dy2dθ²_0 + (2 * dpy2dθ_0 * dNdθ + d²Ndθ² * N² * dy2) * N²) * N
              d²Fxdθ² = d²px2dθ² + d²dx2dθ²
              d²Fydθ² = d²py2dθ² + d²dy2dθ²
            end
            begin
              d²tdθ² = dx1 * d²Fxdθ² + dy1 * d²Fydθ²
              d²Pxdθ² = dx1 * d²tdθ²
              d²Pydθ² = dy1 * d²tdθ²
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
                #@info "p_newton is $p  θ is $θ  t is $t"
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
            px2 = cosθ
            py2 = bsinθ
            dx2 = bcosθ
            dy2 = sinθ
            Fx = px2 + s * dx2 * N
            Fy = py2 + s * dy2 * N
          end
          t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
          begin
            Px = px1 + t * dx1
            Py = py1 + t * dy1
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
      return nothing
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
            dpx2dθ = -sinθ
            dpy2dθ = bcosθ
            dpx2dθ_0 = -bsinθ
            dpy2dθ_0 = cosθ
            ddx2dθ = (dpx2dθ_0 + N² * dNdθ * dx2) * N
            ddy2dθ = (dpy2dθ_0 + N² * dNdθ * dy2) * N
            dFxdθ = dpx2dθ + s * ddx2dθ
            dFydθ = dpy2dθ + s * ddy2dθ
            dtdθ = dx1 * dFxdθ + dy1 * dFydθ
            dPxdθ = dx1 * dtdθ
            dPydθ = dy1 * dtdθ
            dfxdθ = dFxdθ - dPxdθ
            dfydθ = dFydθ - dPydθ
            g = fx * dfxdθ + fy * dfydθ
          end
          # compute the hessian
          begin
            d²px2dθ² = -px2
            d²py2dθ² = -py2
            d²dx2dθ²_0 = -dx2
            d²dy2dθ²_0 = -dy2
            d²dx2dθ² = (d²dx2dθ²_0 + (2 * dpx2dθ_0 * dNdθ + d²Ndθ² * N² * dx2) * N²) * N
            d²dy2dθ² = (d²dy2dθ²_0 + (2 * dpy2dθ_0 * dNdθ + d²Ndθ² * N² * dy2) * N²) * N
            d²Fxdθ² = d²px2dθ² + d²dx2dθ²
            d²Fydθ² = d²py2dθ² + d²dy2dθ²
            d²tdθ² = dx1 * d²Fxdθ² + dy1 * d²Fydθ²
            d²Pxdθ² = dx1 * d²tdθ²
            d²Pydθ² = dy1 * d²tdθ²
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
                #@info "p_newton is $p  θ is $θ  t is $t"
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
            px2 = cosθ
            py2 = bsinθ
            dx2 = bcosθ
            dy2 = sinθ
            Fx = px2 + s * dx2 * N
            Fy = py2 + s * dy2 * N
          end
          t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
          begin
            Px = px1 + t * dx1
            Py = py1 + t * dy1
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
            dpx2dθ = -sinθ
            dpy2dθ = bcosθ
            dpx2dθ_0 = -bsinθ
            dpy2dθ_0 = cosθ
            ddx2dθ = (dpx2dθ_0 + N² * dNdθ * dx2) * N
            ddy2dθ = (dpy2dθ_0 + N² * dNdθ * dy2) * N
            dFxdθ = dpx2dθ + s * ddx2dθ
            dFydθ = dpy2dθ + s * ddy2dθ
            dtdθ = dx1 * dFxdθ + dy1 * dFydθ
            dPxdθ = dx1 * dtdθ
            dPydθ = dy1 * dtdθ
            dfxdθ = dFxdθ - dPxdθ
            dfydθ = dFydθ - dPydθ
            g = fx * dfxdθ + fy * dfydθ
          end
          # compute the hessian
          begin
            d²px2dθ² = -px2
            d²py2dθ² = -py2
            d²dx2dθ²_0 = -dx2
            d²dy2dθ²_0 = -dy2
            d²dx2dθ² = (d²dx2dθ²_0 + (2 * dpx2dθ_0 * dNdθ + d²Ndθ² * N² * dx2) * N²) * N
            d²dy2dθ² = (d²dy2dθ²_0 + (2 * dpy2dθ_0 * dNdθ + d²Ndθ² * N² * dy2) * N²) * N
            d²Fxdθ² = d²px2dθ² + d²dx2dθ²
            d²Fydθ² = d²py2dθ² + d²dy2dθ²
            d²tdθ² = dx1 * d²Fxdθ² + dy1 * d²Fydθ²
            d²Pxdθ² = dx1 * d²tdθ²
            d²Pydθ² = dy1 * d²tdθ²
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
                #@info "p_newton is $p  θ is $θ  t is $t"
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
            px2 = cosθ
            py2 = bsinθ
            dx2 = bcosθ
            dy2 = sinθ
            Fx = px2 + s * dx2 * N
            Fy = py2 + s * dy2 * N
          end
          t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
          begin
            Px = px1 + t * dx1
            Py = py1 + t * dy1
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
            dpx2dθ = -sinθ
            dpy2dθ = bcosθ
            dpx2dθ_0 = -bsinθ
            dpy2dθ_0 = cosθ
            ddx2dθ = (dpx2dθ_0 + N² * dNdθ * dx2) * N
            ddy2dθ = (dpy2dθ_0 + N² * dNdθ * dy2) * N
            dFxdθ = dpx2dθ + s * ddx2dθ
            dFydθ = dpy2dθ + s * ddy2dθ
            dtdθ = dx1 * dFxdθ + dy1 * dFydθ
            dPxdθ = dx1 * dtdθ
            dPydθ = dy1 * dtdθ
            dfxdθ = dFxdθ - dPxdθ
            dfydθ = dFydθ - dPydθ
            g = fx * dfxdθ + fy * dfydθ
          end
          # compute the hessian
          begin
            d²px2dθ² = -px2
            d²py2dθ² = -py2
            d²dx2dθ²_0 = -dx2
            d²dy2dθ²_0 = -dy2
            d²dx2dθ² = (d²dx2dθ²_0 + (2 * dpx2dθ_0 * dNdθ + d²Ndθ² * N² * dx2) * N²) * N
            d²dy2dθ² = (d²dy2dθ²_0 + (2 * dpy2dθ_0 * dNdθ + d²Ndθ² * N² * dy2) * N²) * N
            d²Fxdθ² = d²px2dθ² + d²dx2dθ²
            d²Fydθ² = d²py2dθ² + d²dy2dθ²
            d²tdθ² = dx1 * d²Fxdθ² + dy1 * d²Fydθ²
            d²Pxdθ² = dx1 * d²tdθ²
            d²Pydθ² = dy1 * d²tdθ²
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
                #@info "p_newton is $p  θ is $θ  t is $t"
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
            px2 = cosθ
            py2 = bsinθ
            dx2 = bcosθ
            dy2 = sinθ
            Fx = px2 + s * dx2 * N
            Fy = py2 + s * dy2 * N
          end
          t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
          begin
            Px = px1 + t * dx1
            Py = py1 + t * dy1
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
            dpx2dθ = -sinθ
            dpy2dθ = bcosθ
            dpx2dθ_0 = -bsinθ
            dpy2dθ_0 = cosθ
            ddx2dθ = (dpx2dθ_0 + N² * dNdθ * dx2) * N
            ddy2dθ = (dpy2dθ_0 + N² * dNdθ * dy2) * N
            dFxdθ = dpx2dθ + s * ddx2dθ
            dFydθ = dpy2dθ + s * ddy2dθ
            dtdθ = dx1 * dFxdθ + dy1 * dFydθ
            dPxdθ = dx1 * dtdθ
            dPydθ = dy1 * dtdθ
            dfxdθ = dFxdθ - dPxdθ
            dfydθ = dFydθ - dPydθ
            g = fx * dfxdθ + fy * dfydθ
          end
          # compute the hessian
          begin
            d²px2dθ² = -px2
            d²py2dθ² = -py2
            d²dx2dθ²_0 = -dx2
            d²dy2dθ²_0 = -dy2
            d²dx2dθ² = (d²dx2dθ²_0 + (2 * dpx2dθ_0 * dNdθ + d²Ndθ² * N² * dx2) * N²) * N
            d²dy2dθ² = (d²dy2dθ²_0 + (2 * dpy2dθ_0 * dNdθ + d²Ndθ² * N² * dy2) * N²) * N
            d²Fxdθ² = d²px2dθ² + d²dx2dθ²
            d²Fydθ² = d²py2dθ² + d²dy2dθ²
            d²tdθ² = dx1 * d²Fxdθ² + dy1 * d²Fydθ²
            d²Pxdθ² = dx1 * d²tdθ²
            d²Pydθ² = dy1 * d²tdθ²
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

          #@info "p_newton is $p  θ is $θ  t is $t"
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
            px2 = cosθ
            py2 = bsinθ
            dx2 = bcosθ
            dy2 = sinθ
            Fx = px2 + s * dx2 * N
            Fy = py2 + s * dy2 * N
          end
          t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
          begin
            Px = px1 + t * dx1
            Py = py1 + t * dy1
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
        #update (px1,px2,θ,t,s)
        begin


          px2 = cos(θ)
          py2 = sin(θ) * b_normalized
          dx2 = b_normalized * cos(θ)
          dy2 = sin(θ)
          normF = hypot(dx2, dy2)
          dx2 /= normF
          dy2 /= normF

          det = dx1 * dy2 - dy1 * dx2

          ϵ = 1.0e-10        # TO DO: maybe add this value as a kwargs
          s = NaN        # initial s to NaN
          t = NaN        # initial t to NaN
          if det > ϵ
            t = ((px2 - px1) * dy2 - (py2 - py1) * dx2) / det
            s = ((px2 - px1) * dy1 - (py2 - py1) * dx1) / det
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
      apx[i] = px1 + t * dx1
      apy[i] = py1 + t * dy1
    end
  end
end

"""
  fast_ray_bending!(θ_out::A,adx::A,ady::A,nᵢ::A,nₜ::A,aθmin::A,aθmax::A,ascending::AB;δ=1e-10, kwargs...)  where {AB<:AbstractArray{Bool},A<:AbstractArray{T}} where T<:IEEEFloat

  Compute the ray bending from passing through the atmosphere. This function is a fast version of the ray bending function that modifies in-place
  the arrays without any allocation.

  # Arguments
    - `θ_out     <:Array{IEEEFloat}` : the optimization parameter θ (astromonical angle of the Earth)
    - `adx       <:Array{IEEEFloat}` : the x direction of the ray
    - `ady       <:Array{IEEEFloat}` : the y direction of the ray
    - `nᵢ        <:Array{IEEEFloat}` : the incident index of refraction
    - `nₜ        <:Array{IEEEFloat}` :  the transmitted index of refraction
    - `aθmin     <:Array{IEEEFloat}` : the left angle of the wedge
    - `aθmax     <:Array{IEEEFloat}` : the right angle of the wedge
    - `ascending <:Array{Bool}`      : true defines the ascending part of the ray, finding the first true given the tangent quote

  # Optional arguments
    - `δ::Float64=1e-10` : the machine epsilon used to compute the minimum distance around the wedge ends to avoid numerical instability

  See also: [`fast_ray_bending`](@ref) [`fast_initization_theta`](@ref) [`fast_ray_tracing`](@ref)
"""
function fast_ray_bending!(θ_out::A,adx::A,ady::A,nᵢ::A,nₜ::A,aθmin::A,aθmax::A,ascending::AB; δ=1e-10, kwargs...) where {AB<:AbstractArray{Bool},A<:AbstractArray{T}} where T<:IEEEFloat
  @assert size(adx) == size(ady) "direction coordinates has to be the same size"
  @assert size(apx) == size(apy) "position coordinates has to be the same size"
  @assert size(apx) == size(adx) "position and direction coordinates has to be the same size"
  @assert size(θ_out) == size(adx) "optimization parameter and problem size should be the same"

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
    local dx1 = adx[i]
    local dy1 = ady[i]
    local θ = θ_out[i]  # note: theta is the angle of the earth normal not the angle of the ray, this is important to remember

    norm_ray_direction =  hypot(dx1, dy1)
    dx1 /= norm_ray_direction
    dy1 /= norm_ray_direction
    # if the difference between the top and bottom s is less than the machine tollerance
    # it means that the ray is either inside the wedge or the next wedge has the same atmosphere
    n_incident ≈ n_transmitted && continue

    # the normal to an ellipse can be computed from the normal to the earth
    begin
      dx2 = b*cos(normal_earth)
      dy2 = sin(normal_earth)
      norm = hypot(dx2, dy2)
      dx2 /= norm
      dy2 /= norm
      begin
        begin
          # check if it is intersecting a level or a ray
          if (θ==θmax || θ==θmin)
            # the normal to a ray is defined as (y,-x)
            dx2 = dy2
            dy2 = -tmp
          elseif isAscending
            # the normal is inwards
            dx2 = -dx2
            dy2 = -dy2
          end

          # both directions are already normalized
          local Nx=dx2
          local Ny=dy2

          local direcion_ray_x=dx1
          local direction_ray_y=dy1
          local n01=n_incident/n_transmitted
          local n01²=n01*n01
          local cosθ_incident=-(Nx*direcion_ray_x+Ny*direction_ray_y)
          local sinθ²_transmitted =n01²*(1-cosθ_indicent*cosθ_incident)

          # check if the ray is internally reflected
          # this most likely happens if there is an issue with the atmosphere or if the tangent quote
          # happens to be at a level.
          if sinθ²_transmitted ≤ 1
            dx1= n01*dx1+(n01*cosθ_incident-sqrt(1-sinθ²_transmitted))*Nx
            dy1= n01*dy1+(n01*cosθ_incident-sqrt(1-sinθ²_transmitted))*Ny
          else
            dx1-=2*cosθ_incident*Nx
            dy1-=2*cosθ_incident*Ny
          end
        end

        # update the direction of the ray
        begin
          local norm_new_ray_direction = hypot(dx1, dy1)
          dx1 /= norm_new_ray_direction
          dy1 /= norm_new_ray_direction
          adx[i] = dx1
          ady[i] = dy1
        end
      end
    end
  end
end
##################################################################################### =#
"""
  fast_ray_tracing!(t_out::A,θ_out::A,s_out::A,apx::A,apy::A,adx::A,ady::A,incident_refractive_index::A,aθmin::A,aθmax::A,ascending::AB,
  atm_n::M,atm_θ::V1,atm_h::V2,retrieval_i::RETi,retrieval_j::RETi,retrieval_n::RETf,retrieval_θ::RETf,retrieval_t::RETf,retrieval_h::RETf,
  retrieval_px::RETf,retrieval_py::RETf,retrieval_dx::RETf,retrieval_dy::RETf,
  tangent_quote::A;δ=1e-10,kmax::Int=30,intersection_max::Int=140,initialized=false,free_space::T=T(1), kwargs...
  ) where {V1<:AbstractVector, V2<:AbstractVector, RETi, RETf, M, A<:AbstractArray{T}, AB<:AbstractArray{Bool}} where T<:IEEEFloat


Compute the ray tracing from passing through the atmosphere. This function is a fast version of the ray tracing function that modifies in-place
the arrays without any allocation. This function is a combination of the `fast_ray_bending!` and `fast_ray_tracing!` functions. The function
is optimized to avoid any allocation and to be as fast as possible.

# Arguments

- `t_out     <:Array{IEEEFloat}` : the optimization parameter t (distance from the origin)
- `θ_out     <:Array{IEEEFloat}` : the optimization parameter θ (astromonical angle of the Earth)
- `s_out     <:Array{IEEEFloat}` : the optimization parameter s (distance from the origin)
- `apx       <:Array{IEEEFloat}` : the x position of the ray
- `apy       <:Array{IEEEFloat}` : the y position of the ray
- `adx       <:Array{IEEEFloat}` : the x direction of the ray
- `ady       <:Array{IEEEFloat}` : the y direction of the ray
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
- `retrieval_px <:Array{IEEEFloat}` : the x position of the retrieval
- `retrieval_py <:Array{IEEEFloat}` : the y position of the retrieval
- `retrieval_dx <:Array{IEEEFloat}` : the x direction of the retrieval
- `retrieval_dy <:Array{IEEEFloat}` : the y direction of the retrieval


"""
function fast_ray_tracing!(t_out::A,θ_out::A,s_out::A,apx::A,apy::A,adx::A,ady::A,incident_refractive_index::A,aθmin::A,aθmax::A,ascending::AB,
  atm_n::M,atm_θ::V1,atm_h::V2,retrieval_i::RETi,retrieval_j::RETi,retrieval_n::RETf,retrieval_θ::RETf,retrieval_t::RETf,retrieval_h::RETf,
  retrieval_px::RETf,retrieval_py::RETf,retrieval_dx::RETf,retrieval_dy::RETf,
  tangent_quote::A;δ=1e-10,kmax::Int=30,intersection_max::Int=140,initialized=false,free_space::T=T(1), kwargs...
  ) where {V1<:AbstractVector, V2<:AbstractVector, RETi, RETf, M, A<:AbstractArray{T}, AB<:AbstractArray{Bool}} where T<:IEEEFloat
  @assert size(adx) == size(ady) "Direction coordinates (adx, ady) must have the same size"
  @assert size(apx) == size(apy) "Position coordinates (apx, apy) must have the same size"
  @assert size(apx) == size(adx) "Position (apx, apy) and direction (adx, ady) coordinates must have the same size"
  @assert size(θ_out) == size(t_out) "Optimization parameters θ_out and t_out must have the same size"
  @assert size(t_out) == size(adx) "Optimization parameter t_out and direction coordinates (adx, ady) must have the same size"
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


  NumRays=prod(size(adx))  # size of the rays
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

  #@info "Starting the ray tracing"
  #@info "max_altitude is $max_altitude"
  #@info " extrema of atm_h is $(minimum(atm_h)) and $(maximum(atm_h))"

  rho_max = T(50)
  kloops = div(kmax,5, RoundUp) #Internal loop of the Newton method
  #external loop over all the points
    @inbounds for iter in 1:iter_eff

    # internal loop of ray tracing
      #@batch for idx_rays in eachindex(t_out)
      for idx_rays in eachindex(t_out)
        #@info "Ray $idx_rays and iteration $iter"
        # first iteration is diffent from the rest
        # I need to set the initial value of theta
        # and find the initial wedge, also all incident refractive index are 1
        # also s_top does not exist yet
        local θ = θ_out[idx_rays]                # gibberish the first iteration
        local t = t_out[idx_rays]                # gibberish always
        local directional_sign = ascending[idx_rays] ? T(-1.0) :  T(1.0)
        local px1 = apx[idx_rays]
        local py1 = apy[idx_rays]
        local dx1 = adx[idx_rays]
        local dy1 = ady[idx_rays]
        local θmin = aθmin[idx_rays]                                  # gibberish the first iteration not necessary but I prefer to have it for consistency and debugging
        local θmax = aθmax[idx_rays]                                  # gibberish the first iteration not necessary but I prefer to have it for consistency and debugging
        local nᵢ = incident_refractive_index[idx_rays]                   # gibberish the first iteration
        local nₜ  = T(0)
        local isAscending = ascending[idx_rays]  # initially false
        local s=s_out[idx_rays]                  # gibberish the first iteration
        local px2,py2,dx2,dy2 # normal to the ellipse
        local i_wedge= retrieval_i[idx_rays,iter] # gibberish the first iteration
        local j_wedge= retrieval_j[idx_rays,iter] # gibberish the first iteration
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

        #@info "---------------------------------------------"
        #@info " try to minimize towards $s "
        #@info " top $(s_top) and bottom$(s_bottom))"
        #@info " isAscending is $isAscending"
        #@info "---------------------------------------------"


        retrieval_px[idx_rays,iter]=px1
        retrieval_py[idx_rays,iter]=py1
        retrieval_dx[idx_rays,iter]=dx1
        retrieval_dy[idx_rays,iter]=dy1


        # this value is computed when the wedge is found when the ray tracing is all in one loop
        #local nₜ = nₜ[idx_rays]                   # gibberish the first iteration
        begin
          if iter==1
            if  initialized==false




              # set the initial value of s to be the top of the atmosphere
              s = max_altitude  # assumed to be in a descending order
              nᵢ= free_space    # assumed to be vacuum
              θ = atan(py1,px1) # find the initial value of θ to simplify the computation of the first iteration
              θmin = -Inf
              θmax =  Inf
              s_bottom = max_altitude
              s_top    = max_altitude + 100*ϵ  # just to be sure it is higher and it exists
              retrieval_i[idx_rays,1]=-2 # value of θ
              retrieval_j[idx_rays,1]= 0 # value of h
              retrieval_n[idx_rays,1]=T(1)
              retrieval_θ[idx_rays,1]=NaN
              retrieval_h[idx_rays,1]=NaN




            else

              i_wedge = retrieval_i[idx_rays,1]
              j_wedge = retrieval_j[idx_rays,1]

              s_bottom = atm_h[j_wedge_plus_1]
              s_top    = atm_h[j_wedge]
              s        = isAscending ? s_top : s_bottom
              nᵢ       = retrieval_n[idx_rays,1]
              θ        = retrieval_θ[idx_rays,1]

              θmin     = atm_θ[i_wedge]
              θmax     = atm_θ[i_wedge_plus_1]
            end
          end
        end

        # update indexes
        begin
          i_wedge        = retrieval_i[idx_rays,iter]
          i_wedge_plus_1 = i_wedge+1
          j_wedge        = retrieval_j[idx_rays,iter]
          j_wedge_plus_1 = j_wedge+1
          if IsPeriodic
            i_wedge = mod1(i_wedge,Natm_n)
            i_wedge_plus_1 = mod1(i_wedge+1,Natm_n)
          end
          i_wedge = i_wedge< N_atmn ? i_wedge : -1
          i_wedge_plus_1 = i_wedge_plus_1< N_atmn ? i_wedge_plus_1 : -1
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
          # never reaches here
          if j_wedge==0

          elseif j_wedge==-1

          elseif i_wedge==-2

          end

          if i_wedge==0

          elseif i_wedge==-1

          elseif i_wedge==-2

          end

          retrieval_h[idx_rays,iter+1]=NaN
          retrieval_t[idx_rays,iter+1]=NaN
          retrieval_px[idx_rays,iter+1]=NaN
          retrieval_py[idx_rays,iter+1]=NaN
          retrieval_dx[idx_rays,iter+1]=NaN
          retrieval_dy[idx_rays,iter+1]=NaN
          retrieval_i[idx_rays,iter+1]=i_wedge
          retrieval_j[idx_rays,iter+1]=j_wedge
          retrieval_n[idx_rays,iter+1]=nᵢ
          retrieval_θ[idx_rays,iter+1]=NaN
          continue
        end

        # begin the intersection loop

        begin
          # normalization of the direction
          local norm = hypot(dx1, dy1)
          begin
            dx1 /= norm
            dy1 /= norm
          end

          local origin_times_direction = px1 * dx1 + py1 * dy1
          local cosθ, sinθ, bcosθ, bsinθ, cosθ², sinθ², half_sin2θ, cos2θ
          local R, N, Fx, Fy, t, Px, Py, fx, fy, f, fold
          local dpx2dθ, dpy2dθ, dpx2dθ_0, dpy2dθ_0,  ddx2dθ, ddy2dθ
          local dFxdθ, dFydθ, dtdθ, dPxdθ,  dPydθ,  dfxdθ,  dfydθ
          local g,h,p
          local d²px2dθ²,d²py2dθ², d²dx2dθ²_0
          local d²dy2dθ²_0, d²dx2dθ², d²dy2dθ², d²Fxdθ², d²Fydθ²
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
              px2 = cosθ
              py2 = bsinθ
              dx2 = bcosθ
              dy2 = sinθ
              # point on the normal of the ellipse

              Fx = px2 + s * dx2 * N
              Fy = py2 + s * dy2 * N
            end
            # compute t minimum  for θ and the relative distance squared f
            # dropped 1/2 so that the square root is the displacement
            t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
            # point on the ray
            begin
              Px = px1 + t * dx1
              Py = py1 + t * dy1
            end
            begin
              fx = Fx - Px
              fy = Fy - Py
              f = fx * fx + fy * fy
            end

            # if the distance is already small enough or the function is nan
            # it means I am either at the end point or that the function has reached the minimum

            # set f-> fold
            fold = f
          end

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
                    dpx2dθ = -sinθ
                    dpy2dθ = bcosθ
                    ddx2dθ_0 = -bsinθ
                    ddy2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddx2dθ = (ddx2dθ_0 + N² * dNdθ * dx2) * N
                    ddy2dθ = (ddy2dθ_0 + N² * dNdθ * dy2) * N
                  end
                  begin
                    dFxdθ = dpx2dθ + s * ddx2dθ
                    dFydθ = dpy2dθ + s * ddy2dθ
                    dtdθ = dx1 * dFxdθ + dy1 * dFydθ
                    dPxdθ = dx1 * dtdθ
                    dPydθ = dy1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
                  g = fx * dfxdθ + fy * dfydθ
                end
                # compute the hessian
                begin
                  begin
                    d²px2dθ² = -px2
                    d²py2dθ² = -py2
                    d²dx2dθ²_0 = -dx2
                    d²dy2dθ²_0 = -dy2
                  end
                  begin
                    d²dx2dθ² = (d²dx2dθ²_0 + (2 * ddx2dθ_0 * dNdθ + d²Ndθ² * N² * dx2) * N²) * N
                    d²dy2dθ² = (d²dy2dθ²_0 + (2 * ddy2dθ_0 * dNdθ + d²Ndθ² * N² * dy2) * N²) * N
                    d²Fxdθ² = d²px2dθ² + d²dx2dθ²
                    d²Fydθ² = d²py2dθ² + d²dy2dθ²
                  end
                  begin
                    d²tdθ² = dx1 * d²Fxdθ² + dy1 * d²Fydθ²
                    d²Pxdθ² = dx1 * d²tdθ²
                    d²Pydθ² = dy1 * d²tdθ²
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

begin
  local N₀=N
  local N₀² =N²
local ∂N∂θ = dNdθ
local ∂²N∂θ² = d²Ndθ²
local ∂Px∂θ = dPxdθ
local ∂Py∂θ = dPydθ
local ∂Fx∂θ = dFxdθ
local ∂Fy∂θ = dFydθ
local ∂fx∂θ = dfxdθ
local ∂fy∂θ = dfydθ
local ∂²Px∂θ² = d²Pxdθ²
local ∂²Py∂θ² = d²Pydθ²
local ∂²Fx∂θ² = d²Fxdθ²
local ∂²Fy∂θ² = d²Fydθ²
local ∂²fx∂θ² = d²fxdθ²
local ∂²fy∂θ² = d²fydθ²
local g_step = g
local h_step = h
            #@info s
            #@info ddx2dθ_0
            #@info N² * dNdθ * dx2

            #@info "∂Fx∂θ $∂Fx∂θ"
            #@info "∂Fy∂θ $∂Fy∂θ"
            return
            #@info "∂fx∂θ $∂fx∂θ"
            #@info "∂fy∂θ $∂fy∂θ"
            #@info "∂²Px∂θ² $∂²Px∂θ²"
            #@info "∂²Py∂θ² $∂²Py∂θ²"
            #@info "∂²Fx∂θ² $∂²Fx∂θ²"
            #@info "∂²Fy∂θ² $∂²Fy∂θ²"
            #@info "∂²fx∂θ² $∂²fx∂θ²"
            #@info "∂²fy∂θ² $∂²fy∂θ²"
            #@info "g_step $g_step"
            #@info "h_step $h_step"

            #@info "p_step $(-g/h)"
            #return
end
                # compute the newton step
                p = -g / h
                # compute next θ
                θ = mod2pi(θ + p)

                #@info "p_newton is $p  θ is $θ  t is $t"
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
                  px2 = cosθ
                  py2 = bsinθ
                  dx2 = bcosθ
                  dy2 = sinθ
                  Fx = px2 + s * dx2 * N
                  Fy = py2 + s * dy2 * N
                end
                t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
                begin
                  Px = px1 + t * dx1
                  Py = py1 + t * dy1
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

          # Newton loop
          #@info "Starting the Newton loop θ0 is $θ"
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
                    dpx2dθ = -sinθ
                    dpy2dθ = bcosθ
                    dpx2dθ_0 = -bsinθ
                    dpy2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddx2dθ = (dpx2dθ_0 + N² * dNdθ * dx2) * N
                    ddy2dθ = (dpy2dθ_0 + N² * dNdθ * dy2) * N
                  end
                  begin
                    dFxdθ = dpx2dθ + s * ddx2dθ
                    dFydθ = dpy2dθ + s * ddy2dθ
                    dtdθ = dx1 * dFxdθ + dy1 * dFydθ
                    dPxdθ = dx1 * dtdθ
                    dPydθ = dy1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
                  g = fx * dfxdθ + fy * dfydθ
                end
                # compute the hessian
                begin
                  begin
                    d²px2dθ² = -px2
                    d²py2dθ² = -py2
                    d²dx2dθ²_0 = -dx2
                    d²dy2dθ²_0 = -dy2
                  end
                  begin
                    d²dx2dθ² = (d²dx2dθ²_0 + (2 * dpx2dθ_0 * dNdθ + d²Ndθ² * N² * dx2) * N²) * N
                    d²dy2dθ² = (d²dy2dθ²_0 + (2 * dpy2dθ_0 * dNdθ + d²Ndθ² * N² * dy2) * N²) * N
                    d²Fxdθ² = d²px2dθ² + d²dx2dθ²
                    d²Fydθ² = d²py2dθ² + d²dy2dθ²
                  end
                  begin
                    d²tdθ² = dx1 * d²Fxdθ² + dy1 * d²Fydθ²
                    d²Pxdθ² = dx1 * d²tdθ²
                    d²Pydθ² = dy1 * d²tdθ²
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


                #@info "e²: $e²  θ: $θ cosθ^2: $(cos(θ)^2) R₀: $(1-e² * cos(θ)^2)"
                #@info "N₀  $N"
                #@info "∂Nx∂θ $(dNdθ*N²)"
                #@info "∂²Nx∂θ² $(d²Ndθ² * N² * N²)"
                #@info "∂²Nx∂θ² $dPxdθ"
                #@info "∂Px∂θ $dPxdθ"
                #@info "∂Py∂θ $dPydθ"
                #@info "∂Fx∂θ $dFxdθ"
                #@info "∂Fy∂θ $dFydθ"
                #@info "∂fx∂θ $dfxdθ"
                #@info "∂fy∂θ $dfydθ"
                #@info "∂²Px∂θ² $d²Pxdθ²"
                #@info "∂²Py∂θ² $d²Pydθ²"
                #@info "∂²Fx∂θ² $d²Fxdθ²"
                #@info "∂²Fy∂θ² $d²Fydθ²"
                #@info "∂²fx∂θ² $d²fxdθ²"
                #@info "∂²fy∂θ² $d²fydθ²"
                #@info "g_step $g"
                #@info "h_step $h"
                #@info "p_step $(-g/h)"
                #return


                # compute the newton step
                p = -g / h
                # compute next θ
                θ = mod2pi(θ + p)

                #@info "p_newton is $p  θ is $θ  t is $t"
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
                  px2 = cosθ
                  py2 = bsinθ
                  dx2 = bcosθ
                  dy2 = sinθ
                  Fx = px2 + s * dx2 * N
                  Fy = py2 + s * dy2 * N
                end
                t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
                begin
                  Px = px1 + t * dx1
                  Py = py1 + t * dy1
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
                    dpx2dθ = -sinθ
                    dpy2dθ = bcosθ
                    ddx2dθ_0 = -bsinθ
                    ddy2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddx2dθ = (ddx2dθ_0 + N² * dNdθ * dx2) * N
                    ddy2dθ = (ddy2dθ_0 + N² * dNdθ * dy2) * N
                  end
                  begin
                    dFxdθ = dpx2dθ + s * ddx2dθ
                    dFydθ = dpy2dθ + s * ddy2dθ
                    dtdθ = dx1 * dFxdθ + dy1 * dFydθ
                    dPxdθ = dx1 * dtdθ
                    dPydθ = dy1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
                  g = fx * dfxdθ + fy * dfydθ
                end
                # compute the hessian
                begin
                  begin
                    d²px2dθ² = -px2
                    d²py2dθ² = -py2
                    d²dx2dθ²_0 = -dx2
                    d²dy2dθ²_0 = -dy2
                  end
                  begin
                    d²dx2dθ² = (d²dx2dθ²_0 + (2 * dpx2dθ_0 * dNdθ + d²Ndθ² * N² * dx2) * N²) * N
                    d²dy2dθ² = (d²dy2dθ²_0 + (2 * dpy2dθ_0 * dNdθ + d²Ndθ² * N² * dy2) * N²) * N
                    d²Fxdθ² = d²px2dθ² + d²dx2dθ²
                    d²Fydθ² = d²py2dθ² + d²dy2dθ²
                  end
                  begin
                    d²tdθ² = dx1 * d²Fxdθ² + dy1 * d²Fydθ²
                    d²Pxdθ² = dx1 * d²tdθ²
                    d²Pydθ² = dy1 * d²tdθ²
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
                #@info "p_newton is $p  θ is $θ  t is $t"
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
                  px2 = cosθ
                  py2 = bsinθ
                  dx2 = bcosθ
                  dy2 = sinθ
                  Fx = px2 + s * dx2 * N
                  Fy = py2 + s * dy2 * N
                end
                t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
                begin
                  Px = px1 + t * dx1
                  Py = py1 + t * dy1
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
                    dpx2dθ = -sinθ
                    dpy2dθ = bcosθ
                    dpx2dθ_0 = -bsinθ
                    dpy2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddx2dθ = (dpx2dθ_0 + N² * dNdθ * dx2) * N
                    ddy2dθ = (dpy2dθ_0 + N² * dNdθ * dy2) * N
                  end
                  begin
                    dFxdθ = dpx2dθ + s * ddx2dθ
                    dFydθ = dpy2dθ + s * ddy2dθ
                    dtdθ = dx1 * dFxdθ + dy1 * dFydθ
                    dPxdθ = dx1 * dtdθ
                    dPydθ = dy1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
                  g = fx * dfxdθ + fy * dfydθ
                end
                # compute the hessian
                begin
                  begin
                    d²px2dθ² = -px2
                    d²py2dθ² = -py2
                    d²dx2dθ²_0 = -dx2
                    d²dy2dθ²_0 = -dy2
                  end
                  begin
                    d²dx2dθ² = (d²dx2dθ²_0 + (2 * dpx2dθ_0 * dNdθ + d²Ndθ² * N² * dx2) * N²) * N
                    d²dy2dθ² = (d²dy2dθ²_0 + (2 * dpy2dθ_0 * dNdθ + d²Ndθ² * N² * dy2) * N²) * N
                    d²Fxdθ² = d²px2dθ² + d²dx2dθ²
                    d²Fydθ² = d²py2dθ² + d²dy2dθ²
                  end
                  begin
                    d²tdθ² = dx1 * d²Fxdθ² + dy1 * d²Fydθ²
                    d²Pxdθ² = dx1 * d²tdθ²
                    d²Pydθ² = dy1 * d²tdθ²
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
                #@info "p_newton is $p  θ is $θ  t is $t"
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
                  px2 = cosθ
                  py2 = bsinθ
                  dx2 = bcosθ
                  dy2 = sinθ
                  Fx = px2 + s * dx2 * N
                  Fy = py2 + s * dy2 * N
                end
                t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
                begin
                  Px = px1 + t * dx1
                  Py = py1 + t * dy1
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
                    dpx2dθ = -sinθ
                    dpy2dθ = bcosθ
                    dpx2dθ_0 = -bsinθ
                    dpy2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddx2dθ = (dpx2dθ_0 + N² * dNdθ * dx2) * N
                    ddy2dθ = (dpy2dθ_0 + N² * dNdθ * dy2) * N
                  end
                  begin
                    dFxdθ = dpx2dθ + s * ddx2dθ
                    dFydθ = dpy2dθ + s * ddy2dθ
                    dtdθ = dx1 * dFxdθ + dy1 * dFydθ
                    dPxdθ = dx1 * dtdθ
                    dPydθ = dy1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
                  g = fx * dfxdθ + fy * dfydθ
                end
                # compute the hessian
                begin
                  begin
                    d²px2dθ² = -px2
                    d²py2dθ² = -py2
                    d²dx2dθ²_0 = -dx2
                    d²dy2dθ²_0 = -dy2
                  end
                  begin
                    d²dx2dθ² = (d²dx2dθ²_0 + (2 * dpx2dθ_0 * dNdθ + d²Ndθ² * N² * dx2) * N²) * N
                    d²dy2dθ² = (d²dy2dθ²_0 + (2 * dpy2dθ_0 * dNdθ + d²Ndθ² * N² * dy2) * N²) * N
                    d²Fxdθ² = d²px2dθ² + d²dx2dθ²
                    d²Fydθ² = d²py2dθ² + d²dy2dθ²
                  end
                  begin
                    d²tdθ² = dx1 * d²Fxdθ² + dy1 * d²Fydθ²
                    d²Pxdθ² = dx1 * d²tdθ²
                    d²Pydθ² = dy1 * d²tdθ²
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
                #@info "p_newton is $p  θ is $θ  t is $t"
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
                  px2 = cosθ
                  py2 = bsinθ
                  dx2 = bcosθ
                  dy2 = sinθ
                  Fx = px2 + s * dx2 * N
                  Fy = py2 + s * dy2 * N
                end
                t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
                begin
                  Px = px1 + t * dx1
                  Py = py1 + t * dy1
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
                    dpx2dθ = -sinθ
                    dpy2dθ = bcosθ
                    dpx2dθ_0 = -bsinθ
                    dpy2dθ_0 = cosθ
                  end
                  begin
                    half_dR = e² * half_sin2θ
                    half_dR² = half_dR * half_dR
                    half_d²R = e² * cos2θ
                    N² = N * N
                    dNdθ = -half_dR
                    d²Ndθ² = -half_d²R + 3 * N² * half_dR²
                    ddx2dθ = (dpx2dθ_0 + N² * dNdθ * dx2) * N
                    ddy2dθ = (dpy2dθ_0 + N² * dNdθ * dy2) * N
                  end
                  begin
                    dFxdθ = dpx2dθ + s * ddx2dθ
                    dFydθ = dpy2dθ + s * ddy2dθ
                    dtdθ = dx1 * dFxdθ + dy1 * dFydθ
                    dPxdθ = dx1 * dtdθ
                    dPydθ = dy1 * dtdθ
                    dfxdθ = dFxdθ - dPxdθ
                    dfydθ = dFydθ - dPydθ
                  end
                  g = fx * dfxdθ + fy * dfydθ
                end
                # compute the hessian
                begin
                  begin
                    d²px2dθ² = -px2
                    d²py2dθ² = -py2
                    d²dx2dθ²_0 = -dx2
                    d²dy2dθ²_0 = -dy2
                  end
                  begin
                    d²dx2dθ² = (d²dx2dθ²_0 + (2 * dpx2dθ_0 * dNdθ + d²Ndθ² * N² * dx2) * N²) * N
                    d²dy2dθ² = (d²dy2dθ²_0 + (2 * dpy2dθ_0 * dNdθ + d²Ndθ² * N² * dy2) * N²) * N
                    d²Fxdθ² = d²px2dθ² + d²dx2dθ²
                    d²Fydθ² = d²py2dθ² + d²dy2dθ²
                  end
                  begin
                    d²tdθ² = dx1 * d²Fxdθ² + dy1 * d²Fydθ²
                    d²Pxdθ² = dx1 * d²tdθ²
                    d²Pydθ² = dy1 * d²tdθ²
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
                #@info "p_newton is $p  θ is $θ  t is $t"
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
                  px2 = cosθ
                  py2 = bsinθ
                  dx2 = bcosθ
                  dy2 = sinθ
                  Fx = px2 + s * dx2 * N
                  Fy = py2 + s * dy2 * N
                end
                t = -origin_times_direction + (dx1 * Fx + dy1 * Fy)
                begin
                  Px = px1 + t * dx1
                  Py = py1 + t * dy1
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
              #update (px1,px2,θ,t,s)
              begin
                px2 = cos(θ)
                py2 = sin(θ) * b_normalized
                dx2 = b_normalized * cos(θ)
                dy2 = sin(θ)
                normF = hypot(dx2, dy2)
                dx2 /= normF
                dy2 /= normF

                det = dx1 * dy2 - dy1 * dx2
                ϵ = 1.0e-10        # TO DO: maybe add this value as a kwargs
                s = NaN        # initial s to NaN
                t = NaN        # initial t to NaN
                local Δp12x = px2 - px1
                local Δp12y = py2 - py1


                if abs(det) > ϵ
                  t = (Δp12x * dy2 - Δp12y * dx2) / det
                  s = (Δp12x * dy1 - Δp12y * dx1) / det
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

        #@info "pre: t is $t  s: $s θ: $θ"
        #@info "px1: $px1, py1: $py1, px2: $px2, py2: $py2"
        #@info "dx1: $dx1, dy1: $dy1, dx2: $dx2, dy2: $dy2"
        #return (px1, py1, px2, py2, dx1, dy1, dx2, dy2, θ, t, s)

        ###############
        begin
          if iter==1 && initialized==false
            # find the index of the wedge using binary search
            # assume left orientation of the ray
            # TO DO: add the right orientation of the ray
            i_wedge = findlast(atm_θ.<= θ)
            i_wedge_plus_1 = i_wedge+1
            j_wedge = 1
            j_wedge_plus_1 = 2
            #if s>max_altitude+ϵ  # the first ray never intersected the atmosphere
            #  j_wedge = -2
            #end
            # NOTE: the first intersection can only be downward so dx2 and dy2 are already the outward normal
            # to the wedge


          else
            #@info "--------------------------------------------------"
            #@info "  θ   : $(rad2deg(θ))°  "
            #@info "  θmin: $(rad2deg(θmin))°  "
            #@info "  θmax: $(rad2deg(θmax))°  "
            #@info "  Δh  : $(f)  "
            #@info "  s_top: $(s_top)  "
            #@info "  s_bottom: $(s_bottom)  "
            #@info "  s  : $(s)  "
            #@info "--------------------------------------------------"
            begin
              local tmp = dx2
              if  θmin<θ<θmax
                #@info "Between 2 wedges"
                if  (abs(f)>10^-5 && isAscending==false)
                    #@info " in the middle of the atmosphere"
                    tangent_quote[idx_rays]=s
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
                      dx2 =  dy2
                      dy2 = -tmp
                    end
                    # the ray is going up so I send it to the next ray after registering the
                    # tangent quote
                    #update (px1,px2,θ,t,s)
                    begin
                      px2 = cos(θ)
                      py2 = sin(θ) * b_normalized
                      dx2 = b_normalized * cos(θ)
                      dy2 = sin(θ)
                      normF = hypot(dx2, dy2)
                      dx2 /= normF
                      dy2 /= normF

                      det = dx1 * dy2 - dy1 * dx2
                      ϵ = 1.0e-10        # TO DO: maybe add this value as a kwargs
                      s = NaN        # initial s to NaN
                      t = NaN        # initial t to NaN
                      local Δp12x = px2 - px1
                      local Δp12y = py2 - py1


                      if abs(det) > ϵ
                        t = (Δp12x * dy2 - Δp12y * dx2) / det
                        s = (Δp12x * dy1 - Δp12y * dx1) / det
                      end
                      if (s < 0) # it starts from the ellipse surface, so a negative value is not possible
                        s = NaN
                      end
                    end


                elseif isAscending==true
                  #@info " going up"
                  j_wedge = j_wedge-1  # the direction of h is descending
                  dx2 = -dx2
                  dy2 = -dy2
                else
                  #@info " going down"
                  j_wedge = j_wedge+1  # the direction of h is ascending
                end
              elseif θ==θmax
                #@info "Touching left"
                begin
                  i_wedge = i_wedge+1
                  i_wedge_plus_1 = i_wedge+1
                  dx2 =  dy2
                  dy2 = -tmp
                end
              elseif θ==θmin
                #@info "Touching right"
                begin
                  i_wedge = i_wedge-1
                  dx2 = -dy2
                  dy2 =  tmp
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
              i_wedge = i_wedge< N_atmn ? i_wedge : -1
              i_wedge_plus_1 = i_wedge_plus_1< N_atmn ? i_wedge_plus_1 : -1
              j_wedge = j_wedge< Matm_n ? j_wedge : -1
              j_wedge_plus_1 = j_wedge_plus_1< Matm_n ? j_wedge_plus_1 : -1
            end





          end

          retrieval_i[idx_rays,iter+1]=i_wedge
          retrieval_j[idx_rays,iter+1]=j_wedge
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
        #@info "--------------------------------------------------"
        #@info "Ray $(idx_rays) at $(iter) iteration"
        #@info "--------------------------------------------------"
        #@info "i_wedge=$(i_wedge) j_wedge=$(j_wedge)"
        #@info "i_wedge_plus_1=$(i_wedge_plus_1) j_wedge_plus_1=$(j_wedge_plus_1)"
        #@info "nₜ=$(nₜ) nₜ=$(nₜ)"
        #@info "θmin=$(θmin) θmax=$(θmax)"
        #@info "s_top=$(s_top) s_bottom=$(s_bottom)"
        #@info "px1=$(px1) py1=$(py1)"
        #@info "px2=$(px2) py2=$(py2)"
        #@info "dx1=$(dx1) dy1=$(dy1)"
        #@info "dx2=$(dx2) dy2=$(dy2)"
        #@info "t=$(t) s=$(s)"
        #@info "θ=$(θ)"
        #@info "isAscenging? $(isAscending)"
        #@info "--------------------------------------------------"
        #return
        ##################################################
        # bending the ray
        if !(nᵢ==nₜ) # do bend only if the refractive index are different
          #@info " --------------------------------------------------"
          #@info "Bending the ray"
          #@info " --------------------------------------------------"
          let
            # check if it is intersecting a level or a ray
            # both directions are already normalized
            local n_incident = nᵢ
            local n_transmitted = nₜ
            local direcion_ray_x=dx1
            local direction_ray_y=dy1
            local Nx=dx2
            local Ny=dy2
            local n01=n_incident/n_transmitted
            local n01²=n01*n01
            local cosθ_incident=-(Nx*direcion_ray_x+Ny*direction_ray_y)
            local sinθ²_transmitted =n01²*(1-cosθ_incident*cosθ_incident)

            # check if the ray is internally reflected
            # this most likely happens if there is an issue with the atmosphere or if the tangent quote
            # happens to be at a level.
            if sinθ²_transmitted ≤ 1
              dx1= n01*dx1+(n01*cosθ_incident-sqrt(1-sinθ²_transmitted))*Nx
              dy1= n01*dy1+(n01*cosθ_incident-sqrt(1-sinθ²_transmitted))*Ny
            else
              dx1-=2*cosθ_incident*Nx
              dy1-=2*cosθ_incident*Ny
            end

          end
          #@info " dx_new=$(dx1) dy_new=$(dy1)"

          #@info " --------------------------------------------------"


        end
        # update the direction of the ray
        begin
          local norm_new_ray_direction = hypot(dx1, dy1)
          dx1 /= norm_new_ray_direction
          dy1 /= norm_new_ray_direction
          # first update the position using the array cause dx1,dy1 are already modified
          apx[idx_rays] = px1+t*adx[idx_rays]
          apy[idx_rays] = py1+t*ady[idx_rays]
          px1=apx[idx_rays]
          py1=apy[idx_rays]
          adx[idx_rays] = dx1
          ady[idx_rays] = dy1
          t_out[idx_rays] = t
          θ_out[idx_rays] = θ
          # set up next iteration
          s_out[idx_rays] = isAscending ? s_top : s_bottom
          aθmin[idx_rays] = θmin
          aθmax[idx_rays] = θmax
          incident_refractive_index[idx_rays]  = nₜ
          ascending[idx_rays] = isAscending

        end

        # update output arrays
        begin
          # fill the retrieval output for the current iteration
          retrieval_i[idx_rays,iter+1]=i_wedge
          retrieval_j[idx_rays,iter+1]=j_wedge
          retrieval_θ[idx_rays,iter+1]=θ
          retrieval_t[idx_rays,iter+1]=t
          retrieval_h[idx_rays,iter+1]=s

        end
    end
    # early stop condition
    if number_rays_stopped == NumRays
      break
    end
  end
end
