

const OUTWARDS = true
const INWARDS = false
const CLOCKWISE = true
const COUNTCLOCKWISE = false
@inline _normal(nx::T,ny::T,  outwards::Bool)   where T =  outwards  ==OUTWARDS   ? (nx,ny) : (-ny,-nx)
@inline _tangent(nx::T,ny::T, clockwise::Bool) where T =   clockwise ==CLOCKWISE ? (ny,-nx) : (-ny,nx)

@inline function corrected_compute_newtop_step(s::T, dx1::T, dy1::T, cosθ::T, sinθ::T, Fx::T,Fy::T,Px::T,Py::T, b_normalized::T,e²::T,N₀::T)::T where T
  local cos2θ = cosθ^2 - sinθ^2
  local half_sin2θ = sinθ * cosθ
  local N₀² = N₀ * N₀
  # f=1/2*((Fx-Px)^2 + (Fy-Py)^2)
  # compute dNdθ
  local b  = b_normalized
  local px = cosθ
  local py = sinθ
  local ∂px∂θ = -sinθ
  local ∂py∂θ =  cosθ
  local ∂²px∂θ² = -px
  local ∂²py∂θ² = -py

  local s₀      = s
  local s₀b     = s₀ * b
  local ∂R₀∂θ   = e² * half_sin2θ
  local ∂²R₀∂θ² = e² * cos2θ
  local ∂N∂θ    = -  ∂R₀∂θ
  local ∂²N∂θ²  = (T(3) *N₀² * ∂R₀∂θ^2 -  ∂²R₀∂θ²)


  # First derivatives of Fx and Fy
  local  ∂Fx∂θ = ∂px∂θ   +  s₀b*N₀*(∂px∂θ+px*∂N∂θ*N₀²)
  local  ∂Fy∂θ = ∂py∂θ*b +  s₀ *N₀*(∂py∂θ+py*∂N∂θ*N₀²)

  local  ∂²Fx∂θ² = ∂²px∂θ²+s₀b*N₀*(∂²px∂θ²+(2*∂px∂θ*∂N∂θ + px*N₀² * ∂²N∂θ²)*N₀²)
  local  ∂²Fy∂θ² = ∂²py∂θ²+s₀ *N₀*(∂²py∂θ²+(2*∂py∂θ*∂N∂θ + py*N₀² * ∂²N∂θ²)*N₀²)

      # Compute t and its derivatives
  local ∂t∂θ = dx1 * ∂Fx∂θ + dy1 * ∂Fy∂θ
  local ∂²t∂θ² = dx1 * ∂²Fx∂θ² + dy1 * ∂²Fy∂θ²
  local ∂Px∂θ = dx1 * ∂t∂θ
  local ∂Py∂θ = dy1 * ∂t∂θ

  local ∂²Px∂θ² = dx1 * ∂²t∂θ²
  local ∂²Py∂θ² = dy1 * ∂²t∂θ²

  local fx= Fx - Px
  local fy= Fy - Py
  local ∂fx∂θ = ∂Fx∂θ - ∂Px∂θ
  local ∂fy∂θ = ∂Fy∂θ - ∂Py∂θ
  local ∂²fx∂θ² = ∂²Fx∂θ² - ∂²Px∂θ²
  local ∂²fy∂θ² = ∂²Fy∂θ² - ∂²Py∂θ²
  local g_step= (fx * ∂fx∂θ + fy * ∂fy∂θ)
  local h_step= (∂fx∂θ^2+∂fy∂θ^2 + fx * ∂²fx∂θ² + fy * ∂²fy∂θ²)

  local λ = 10^-10
#=
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
=#

  h_step = abs(h_step)>λ ? h_step : λ  # sufficient positiveness of the hessian

  return -g_step / h_step
end
@inline function _update_earth_ray(θ::T,b::T,e²::T,ϵ::T)::Tuple{T,T,T,T,T} where T
  return (cos(θ),b*sin(θ),b*cos(θ),sin(θ),_compute_N₀(θ,e²,ϵ))
end

@inline function _update_angle(θ::T,p_newton::T)::T where T
   return mod2pi(θ + p_newton)
end


@inline function _clampθ(θ::T,θmin::T,θmax::T)::T where T
  # simple case
  @debug "θmin: $(θmin) θ: $θ  θmax: $(θmax)"
  θmin<=θmax && return clamp(θ,θmin,θmax)
  @debug " [θmin,2pi) U [0,θmax]"
  θ>θmin || θ<θmax && return θ
  local distance_from_theta_min = mod(θ - θmin, 2π)
  local distance_from_theta_max = mod(θmax - θ, 2π)
  return distance_from_theta_min < distance_from_theta_max ? θmin : θmax
end


@inline function _bend_ray!(direction_ray_x::T,direction_ray_y::T,Nx::T,Ny::T,n_incident::T,n_transmitted::T)::Tuple{T,T} where T

  local norm_ray_direction = hypot(direction_ray_x,direction_ray_y)
  local norm_normal        = hypot(Nx,Ny)
  direction_ray_x/=norm_ray_direction
  direction_ray_y/=norm_ray_direction
  Nx/=norm_normal
  Ny/=norm_normal
  local n01=n_incident/n_transmitted
  local n01²=n01*n01
  local cosθ_incident=-(Nx*direction_ray_x+Ny*direction_ray_y)
  local sinθ²_transmitted =n01²*(1-cosθ_incident*cosθ_incident)

  if sinθ²_transmitted ≤ 1
    direction_ray_x= n01*direction_ray_x+(n01*cosθ_incident-sqrt(1-sinθ²_transmitted))*Nx
    direction_ray_y= n01*direction_ray_y+(n01*cosθ_incident-sqrt(1-sinθ²_transmitted))*Ny
  else
    direction_ray_x-=2*cosθ_incident*Nx
    direction_ray_y-=2*cosθ_incident*Ny
    #readline()
  end
  return (direction_ray_x,direction_ray_y)
end

@inline function  _interception_two_rays(x1::T,y1::T,x2::T,y2::T,dx1::T,dx2::T, dy1::T, dy2::T)::Tuple{T,T} where T
  local ϵ = 1.0e-10
  local d1_cross_d2 = dx1 * dy2 - dy1 * dx2 # determinant
  local Δx = x2 - x1
  local Δy = y2 - y1
  # check if the rays are parallel
  if abs(d1_cross_d2) < ϵ
    # note: s<0 is impossible in this context because it means the ray is inside the ellipse
    return T(-1),T(-1)
  end
  local t = (Δx * dy2 - Δy * dx2) / d1_cross_d2
  local s = (Δx * dy1 - Δy * dx1) / d1_cross_d2
  return t,s
end

@inline _check_if_stopped(i_wedge_right::Int,j_wedge_top)::Bool= j_wedge_top<1 || i_wedge_right<1


@inline function _compute_N₀(θ::T,e²::T,ϵ::T)::T where T
  local cosθ = cos(θ)
  local cosθ² = cosθ * cosθ
  local R₀ = max(T(1) - e² * cosθ²,ϵ)
  return T(1)/√(R₀)
end
@inline _increment_penality(ρ::T, γ::T, ρ_max::T) where T = min(ρ*γ,ρ_max)
@inline _distance²(cx::T,cy::T) where T = cx * cx + cy * cy
@inline _compute_minimum_t(origin_time_direction_i::T, dx_i::T, dy_j::T, pointx_j::T, pointy_j::T) where T = -origin_time_direction_i + (dx_i * pointx_j + dy_j * pointy_j)
@inline _compute_ray_at__t(point_i, direction_i, t_i) = point_i + t_i * direction_i
@inline function _top_skip_condition(i::Int, j::Int,iter::Int)
  iter==1    &&    return false # skip the first iteration
  @debug "iter greater than 1"
  (j>0 || i>0) && return false # if any of the index is negative, return false
  @debug "iter greater than 1 and j or i are negative"
  return true
end

@inline function _update_index_i(i::Int, N::Int, isPeriodic:: Bool)
  if isPeriodic
    return mod1(i,N)
  else
    return _update_index_j(i,N)
  end
end
@inline _update_index_j(j::Int, M::Int)= j<=M ? j : -1



"""
  marco_ray_tracing!(t_out::A,θ_out::A,s_out::A,apx::A,apy::A,adx::A,ady::A,incident_refractive_index::A,aθmin::A,aθmax::A,ascending::AB,
  atm_n::M,atm_θ::V1,atm_h::V2,retrieval_i::RETi,retrieval_j::RETi,retrieval_n::RETf,retrieval_θ::RETf,retrieval_t::RETf,retrieval_h::RETf,
  retrieval_px::RETf,retrieval_py::RETf,retrieval_dx::RETf,retrieval_dy::RETf,
  tangent_quote::A;δ=1e-10,kmax::Int=30,intersection_max::Int=140,initialized=false,free_space::T=T(1),rho_max::T=T(50.0) kwargs...
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
function marco_ray_tracing!(t_out::A,θ_out::A,s_out::A,apx::A,apy::A,adx::A,ady::A,incident_refractive_index::A,aθmin::A,aθmax::A,ascending::AB,
  atm_n::M,atm_θ::V1,atm_h::V2,retrieval_i::RETi,retrieval_j::RETi,retrieval_n::RETf,retrieval_θ::RETf,retrieval_t::RETf,retrieval_h::RETf,
  retrieval_px::RETf,retrieval_py::RETf,retrieval_dx::RETf,retrieval_dy::RETf,
  tangent_quote::A;δ=1e-10,kmax::Int=30,
  intersection_max::Int=140,initialized=false,free_space::T=T(1), kwargs...
  )::Nothing where {V1<:AbstractVector, V2<:AbstractVector, RETi, RETf, M, A<:AbstractArray{T}, AB<:AbstractArray{Bool}} where T<:IEEEFloat
  # consistency checks
  begin
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
  end
    # I don't need it if using StaticArrays but I do not want to make assumptions
  local Natm_n=size(atm_n,1)
  local Matm_n=size(atm_n,2)
  local Natm_θ=size(atm_θ,1)
  local Matm_h=size(atm_h,1)
  local IsPeriodic = Natm_n == Natm_θ

  local NumRays=prod(size(adx))  # size of the rays
  local NumRaysRetrieval= size(retrieval_i,1) # size of the retrieval rays

  local Miter=size(retrieval_i,2) # number of iterations
  local iter_eff= min(Miter-1,intersection_max)



  @assert (Natm_θ-1)<=Natm_n<=Natm_θ "atm_n has $Natm_n rows but it has to have either the length of atm_θ-1, $(Natm_θ-1),  or  $(Natm_θ) (periodic radial atmosphere)"
  @assert Matm_n==Matm_h-1 "atm_n has $Matm_n columns but it has to have the length of atm_h-1, $(Matm_h-1)"
  @assert size(retrieval_i,1) == NumRays "retrieval has to be have the number of rows $(NumRaysRetrieval) equivalent to the number of rays ($NumRays)"

  # ellipse parameters
  local  b_normalized = get_minoraxis(T)
  local e² = get_e²(T)
  local ϵ  = 1.0e-10
  # early stop condition
  local number_rays_stopped=0
  local max_altitude=atm_h[1]

  @debug "Starting the ray tracing"
  @debug "max_altitude is $max_altitude"
  @debug " extrema of atm_h is $(minimum(atm_h)) and $(maximum(atm_h))"

  @inbounds for iter in 1:iter_eff
    for idx_rays in eachindex(t_out)
      @debug "######################################"
      @debug "Ray $idx_rays and iteration $iter"

        #######################################
        # Initialize local variables
        #######################################
        ## Ray
        # length of the ray and its derivatives wrt θ
        local t = T(0)
        # direction and location
        local px1 = apx[idx_rays]
        local py1 = apy[idx_rays]
        local dx1 = adx[idx_rays]
        local dy1 = ady[idx_rays]

        # point 1 does not depend on theta only direction 1
        local Px = T(0)
        local Py = T(0)

        # normalization of the direction
        local norm_d1 = hypot(dx1, dy1)
        dx1 /= norm_d1
        dy1 /= norm_d1

        local origin_times_direction_1 = px1 * dx1 + py1 * dy1

        # function to minimize and its derivatives
        local f = T(0)
        local fold = T(0)


        ###########################################################################################
        # First initialization
        # I can assume this has been taken care during the end of the code after the first iteration
        local i_wedge_right  = retrieval_i[idx_rays,iter]
        local j_wedge_top    = retrieval_j[idx_rays,iter] # gibberish the first iteration
        ############################################################################################
        local i_wedge_left   = _update_index_i(i_wedge_right+1,Natm_n,IsPeriodic)
        local j_wedge_bottom = _update_index_j(j_wedge_top+1,Matm_n)
        local isAscending    = ascending[idx_rays]  # initially false
        #################################################
        # Initialize local variables
        # @turbo and @batch needs to know at compilation time
        #################################################

        local isAscending    = iter           > 1 ?  ascending[idx_rays] : false  # initially false

        local s_bottom  = j_wedge_top       > 0 ? atm_h[j_wedge_bottom]    : max_altitude
        local s_top     = j_wedge_bottom    > 0 ? atm_h[j_wedge_top] : T(-99)
        local nᵢ        = iter           > 1 ? incident_refractive_index[idx_rays] : free_space
        local s         = iter           > 1 ? s_out[idx_rays]   :  max_altitude
        local θ         = iter           > 1 ? θ_out[idx_rays]   : atan(py1,px1)  # initial guess
        local θmin      = iter           > 1 ? aθmin[idx_rays]   : -Inf
        local θmax      = iter           > 1 ? aθmax[idx_rays]   :  Inf
        local is_s_top  = s_top*isAscending
        local is_s_bottom = s_bottom*!isAscending

        @debug "s: $s s_top : $is_s_top is s bottom : $is_s_bottom   isAscending: $isAscending"
        local s_target  = iter           > 1 ? s_top*isAscending+ s_bottom*!isAscending : max_altitude
        @debug "s_target: $s_target  $(is_s_top+is_s_bottom)"

        local current_quote_tangent = iter >1 ? tangent_quote[idx_rays] : T(Inf)
        local θ_old = θ
        local θ_previous = θ
        local s_old = s


        if iter==1
          retrieval_i[idx_rays,iter]=i_wedge_right
          retrieval_j[idx_rays,iter]=j_wedge_top
          retrieval_θ[idx_rays,iter]=θ
          retrieval_t[idx_rays,iter]=T(0)
          retrieval_h[idx_rays,iter]=s_target
          retrieval_n[idx_rays,iter]=nᵢ
          retrieval_px[idx_rays,iter]=px1
          retrieval_py[idx_rays,iter]=py1
          retrieval_dx[idx_rays,iter]=dx1
          retrieval_dy[idx_rays,iter]=dy1
        end


        #################################################
        # DEBUG
        @debug "---------------------------------------------"
        @debug " try to minimize towards $s_target "
        @debug " top $(s_top) and bottom$(s_bottom))"
        @debug "j_wedge_top   : $j_wedge_top     s_top: $(s_top)"
        @debug "j_wedge_bottom: $j_wedge_bottom  s_bottom: $(s_bottom)"

        @debug " isAscending is $isAscending"
        @debug "---------------------------------------------"
        #################################################

        local px1_current = px1
        local py1_current = py1
        local dx1_current = dx1
        local dy1_current = dy1
        local i_previous = i_wedge_right
        local j_previous = j_wedge_top

        # if iter>1 and any index is negative


        if _top_skip_condition(min(i_wedge_left,i_wedge_right),min(j_wedge_top,j_wedge_bottom),iter)
            retrieval_h[idx_rays,iter+1]=T(-999)
            retrieval_t[idx_rays,iter+1]=T(-999)
            retrieval_px[idx_rays,iter+1]=px1
            retrieval_py[idx_rays,iter+1]=px1
            retrieval_dx[idx_rays,iter+1]=dx1
            retrieval_dy[idx_rays,iter+1]=dx1
            retrieval_i[idx_rays,iter+1]=i_wedge_right
            retrieval_j[idx_rays,iter+1]=j_wedge_top
            retrieval_n[idx_rays,iter+1]=nᵢ
            retrieval_θ[idx_rays,iter+1]=θ
            continue
        end

        local px2 = T(0)
        local py2 = T(0)
        local dx2 = T(0)
        local dy2 = T(0)

        begin
          local N₀ = T(0)

          (px2,py2,dx2,dy2,N₀) = _update_earth_ray(θ,b_normalized,e²,ϵ)
          local Fx   =  _compute_ray_at__t(px2, dx2*N₀,s_target)
          local Fy   =  _compute_ray_at__t(py2, dy2*N₀,s_target)
          t     = _compute_minimum_t(origin_times_direction_1, dx1, dy1, Fx, Fy)


          Px = _compute_ray_at__t(px1, dx1, t)
          Py = _compute_ray_at__t(py1, dy1, t)

          local fx= Fx - Px
          local fy= Fy - Py
          local k_iter = 0
          f = sqrt(_distance²(fx,fy))
          fold = f

          θ_old = θ
            # Newton loop
          @debug "Starting the Newton loop θ0 is $θ"
          for k_iter in 1:kmax
              p_newton =  corrected_compute_newtop_step(s_target, dx1, dy1, px2, dy2, Fx,Fy,Px,Py, b_normalized,e²,N₀ )

              θ = _update_angle(θ,p_newton)
              (px2,py2,dx2,dy2,N₀) = _update_earth_ray(θ,b_normalized,e²,ϵ)

              Fx     =  _compute_ray_at__t(px2, dx2*N₀,s_target)
              Fy     =  _compute_ray_at__t(py2, dy2*N₀,s_target)
              t      = _compute_minimum_t(origin_times_direction_1, dx1, dy1, Fx, Fy)
              Px     = _compute_ray_at__t(px1, dx1, t)
              Py     = _compute_ray_at__t(py1, dy1, t)
              fx= Fx - Px
              fy= Fy - Py
              f = sqrt(_distance²(fx,fy))

              #@debug "p_newton is $p_newton  θ is $θ  t is $t"



              if abs(f-fold)<δ  && abs(θ-θ_old)<δ
                @debug "Convergence reached in $k_iter iterations"
                @debug "f is $f and fold is $fold"
                @debug "θ is $θ and θold id $θ_old"
                @debug "θ is $θ and θold id $θ_previous"
                break
              end
              fold = f
              θ_old = θ
          end
          # compute the distance between the new point Px,Px and the ellipse with origin
          # at px2,py2 and the normal Fx,Fy
          # I know that the minimum distance between the ellipse and the point of interest is  the ⟂ distance
          # so even if it does not converge exactly to the ellipse, I can see if I am above below or stuck


          # I am moving towards the ellipse
          # 1. check if t>0
          if t<=0
            t=T(0)
            θ=θ_old
            s=s_old
          end

          # 2. check if it is in the bounds of θmin and θmax
          dx2=N₀*dx2
          dy2=N₀*dy2
          @debug "(t,s) $(t), $(s)"
          θ=_clampθ(θ,θmin,θmax)

          @debug "Pre (t,s) $(t), $(s)"
            @debug "s: $s s_target $s_target Δs: $(s-s_target)"
          (t,s) = _interception_two_rays(px1,py1,px2,py2,dx1,dx2,dy1,dy2)

          @debug "Post (t,s) $(t), $(s)"
            @debug "s: $s s_target $s_target Δs: $(s-s_target)"
          if (θ==θ_previous && t==0)
            @debug "θ is the same"

            isAscending = !isAscending

          else

            local nx = dx2
            local ny = dy2
            # 3. compute corrected (t,s)
            (px2,py2,dx2,dy2,N₀) = _update_earth_ray(θ,b_normalized,e²,ϵ)
            dx2*=N₀
            dy2*=N₀
            @debug "Pre: (t,s) is $t and $s"
            @debug "s: $s_target: $s_target Δs: $(s-s_target)"

            if θ==θmax
              @debug "θ=θmax"
              i_wedge_right = _update_index_i(i_wedge_right+1,Natm_n,IsPeriodic)
              i_wedge_left  = _update_index_i(i_wedge_right+1,Natm_n,IsPeriodic)

              @debug " nx: $dx2 ny: $dy2"
              (dx2,dy2) = _tangent(nx,ny,CLOCKWISE)

              @debug " tx: $dx2 ty: $dy2"
            elseif θ==θmin

              @debug "θ=θmin"

              @debug " nx: $dx2 ny: $dy2 tx: $(-dy2) ty: $(dx2)"
              i_wedge_right = _update_index_i(i_wedge_right-1,Natm_n,IsPeriodic)
              i_wedge_left = _update_index_i( i_wedge_right+1,Natm_n,IsPeriodic)

              (dx2,dy2) = _tangent(nx,ny,COUNTCLOCKWISE)
            elseif abs(s-s_target)<10*δ

              @debug "new_level"
              j_step = isAscending ? -1 : 1
              j_wedge_top = _update_index_j(j_wedge_top+j_step,Matm_n)
              j_wedge_bottom = _update_index_j(j_wedge_top+1,Matm_n)

              if isAscending
                (dx2,dy2) = _normal(nx,ny,INWARDS)

              else
                (dx2,dy2) = _normal(nx,ny,OUTWARDS)
              end
            else
              isAscending = !isAscending
              @debug "switch"
              @debug "Δs: $abs(s-s_target)"
              #readline()
            end
            #readline()
            @debug "Post: (t,s) is $t and $s"
          end
          @debug "isAscending is $isAscending"
          @debug "abs(s-s_target)<1e-3 is $(abs(s-s_target)<1e-3)"
          @debug "θ==θmax is $(θ==θmax)"
          @debug "θ==θmin is $(θ==θmin)"

          ## update position
          px1 = _compute_ray_at__t(px1, dx1, t)
          py1 = _compute_ray_at__t(py1, dy1, t)

        end

        if iter==1

          i_wedge_right = findlast(atm_θ.<= θ)
          j_wedge_top = 1
          i_wedge_left= _update_index_i(i_wedge_right+1,Natm_n,IsPeriodic)
          j_wedge_bottom = _update_index_j(j_wedge_top+1,Matm_n)
          @debug "i: $i_wedge_right  i+1: $i_wedge_left"
          @debug "j: $j_wedge_top   j+1: $j_wedge_bottom"
        end
        ###############################################################################
        #  DEBUG
        ###############################################################################
        @debug "---------------------------------------------"
        @debug "i is $i_wedge_right and j is $j_wedge_top"
        @debug "t is $t"
        @debug "theta is $θ"
        @debug "s_previous: $s_old"
        @debug "s_current: $s"
        @debug  "p_previuous: $px1_current, $py1_current"
        @debug "p_current: $px1, $py1"
        @debug "d_previuous: $dx1_current, $dy1_current"
        @debug "d_current: $dx1, $dy1"
        @debug "i: $i_wedge_right and j: $j_wedge_top"
        @debug "---------------------------------------------"
        px1_current = px1
        py1_current = py1
        dx1_current = dx1
        dy1_current = dy1

        ###############################################################################

        if _top_skip_condition(min(i_wedge_right,i_wedge_left),min(j_wedge_top,j_wedge_bottom),2) # avoid the first iteration
          @debug "Ray $idx_rays is stopped"
          @debug "i_wedge_right = $i_wedge_right and j_wedge_top = $j_wedge_top"

          number_rays_stopped+=1
          continue
        end


        begin
          @debug "------------- Bending the ray ---------------"
          @debug "i: $i_wedge_right  i+1: $i_wedge_left"
          @debug "j: $j_wedge_top   j+1: $j_wedge_bottom"
          @debug "i_previous: $i_previous "
          @debug "j_previous: $j_previous  j: $j_wedge_top"
          @debug "j: $j_wedge_top   j+1: $j_wedge_bottom"
          local nₜ = atm_n[i_wedge_right,j_wedge_top]
          @debug "nₜ is $nₜ and nᵢ is $nᵢ"
          if !(nᵢ==nₜ) # do bend only if the refractive index are different
            local direction_x = dx1
            local direction_y = dy1
            local normal_x = dx2
            local normal_y = dy2
            local n_incident = nᵢ
            local n_transmitted = nₜ
            (dx1,dy1) =_bend_ray!(direction_x,direction_y,normal_x,normal_y,n_incident,n_transmitted)
          end
        end
        begin
          local norm_new_ray_direction = hypot(dx1, dy1)
          dx1 /= norm_new_ray_direction
          dy1 /= norm_new_ray_direction
          # first update the position using the array cause dx1,dy1 are already modified
          apx[idx_rays] = px1
          apy[idx_rays] = py1
          adx[idx_rays] = dx1
          ady[idx_rays] = dy1
          t_out[idx_rays] = t
          θ_out[idx_rays] = θ
          # set up next iteration
          s_out[idx_rays] = s
          @debug "i is $i_wedge_right  i+1 is $i_wedge_left"
          @debug "j is $j_wedge_top  j+1 is $j_wedge_bottom"
          aθmin[idx_rays] = atm_θ[i_wedge_right]
          aθmax[idx_rays] = atm_θ[i_wedge_left]
          @debug "θmin is $(aθmin[idx_rays]) and θmax is $(aθmax[idx_rays])"
          incident_refractive_index[idx_rays]  = atm_n[i_wedge_right,j_wedge_top]
          ascending[idx_rays] = isAscending

        end

        # update output arrays
        begin
          # fill the retrieval output for the current iteration
          retrieval_i[idx_rays,iter+1]=i_wedge_right
          retrieval_j[idx_rays,iter+1]=j_wedge_top
          retrieval_θ[idx_rays,iter+1]=θ
          retrieval_t[idx_rays,iter+1]=t
          retrieval_h[idx_rays,iter+1]=s
          retrieval_n[idx_rays,iter+1]=incident_refractive_index[idx_rays]
          retrieval_px[idx_rays,iter+1]=px1
          retrieval_py[idx_rays,iter+1]=py1
          retrieval_dx[idx_rays,iter+1]=dx1
          retrieval_dy[idx_rays,iter+1]=dy1
          current_quote_tangent = max(0,min(current_quote_tangent,s))  #update tangent quote
          tangent_quote[idx_rays] = current_quote_tangent
          @debug "current_quote_tangent is $current_quote_tangent"
          @debug "s is $s"
          @debug "i is $i_wedge_right and j is $j_wedge_top"
          @debug "t is $t"
          @debug "---------------------------------------------"

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
  (nx,ny)=(nx,ny)|> x-> x./hypot(x...) .*inwardoutward
 #################################
 angle= ang

 dir=_rotation_matrix(angle)*[nx,ny]
 return (dir[1],dir[2])
end
