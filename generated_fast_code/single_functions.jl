@inline function _compute_newtop_step(s::T, dx1::T, dy1::T, cosθ::T, sinθ::T, Fx::T,Fy::T,Px::T,Py::T, b_normalized::T,e²::T,N₀::T )::T where T
  local cos2θ = cosθ * cosθ - sinθ * sinθ
  local sin2θ = 2 * sinθ * cosθ
  local N₀² = N₀ * N₀
  # f=1/2*((Fx-Px)^2 + (Fy-Py)^2)
  # compute dNdθ
  local b  = b_normalized
  local s₀ = s
  local s₀b = s₀ * b
  local ∂N∂θ   = - T(0.5) * e² * sin2θ * N₀² * N₀
  local ∂²N∂θ² = -e² *N₀² * N₀ *(cos2θ  + T(3) *N₀² * e² *e² * sin2θ^2)
  local parenthesis1 =( T(1)   + s₀b * N₀)
  local parenthesis2 =( s₀b    +       N₀)

  # First derivatives of Fx and Fy
  local  ∂Fx∂θ = -sinθ *parenthesis1+ cosθ * s₀b * ∂N∂θ
  local  ∂Fy∂θ =  cosθ *parenthesis2 + sinθ * ∂N∂θ

  local  ∂²Fx∂θ² = -cosθ * parenthesis1 -2 * sinθ *s₀b* ∂N∂θ + cosθ*s₀b* ∂²N∂θ²
  local  ∂²Fy∂θ² = -sinθ * parenthesis2 - 2 * cosθ * s₀ * b_normalized * ∂N∂θ + s₀ * b_normalized * ∂²N∂θ²

      # Compute t and its derivatives
  local ∂t∂θ = dx1 * ∂Fx∂θ + dy1 * ∂Fy∂θ
  local ∂²t∂θ² = dx1 * ∂²Fx∂θ² + dy1 * ∂²Fy∂θ²
  local ∂Px∂θ = dx1 * ∂t∂θ
  local ∂Py∂θ = dy1 * ∂t∂θ

  local ∂²Px∂θ² = dx1 * ∂²t∂θ²
  local ∂²Py∂θ² = dy1 * ∂²t∂θ²

  local Δfx= Fx - Px
  local Δfy= Fy - Py
  local ∂fx∂θ = ∂Fx∂θ - ∂Px∂θ
  local ∂fy∂θ = ∂Fy∂θ - ∂Py∂θ
  local ∂²fx∂θ² = ∂²Fx∂θ² - ∂²Px∂θ²
  local ∂²fy∂θ² = ∂²Fy∂θ² - ∂²Py∂θ²
  local g_step= 2*Δfx * ∂fx∂θ + 2*Δfy * ∂fy∂θ
  local h_step= 2 * (∂fy∂θ^2 + Δfx * ∂²fx∂θ² + Δfy * ∂²fy∂θ²)

  local λ = 5*10^-5
  local ϵ = 10^-5

  h_step = (h_step - λ) > ϵ ? h_step :  λ  # sufficient positiveness of the hessian

  return -g_step / h_step
end

@inline function _update_angle(θ::T,p_newton::T)::T where T
   return mod2pi(θ + p_newton)
end

@inline function _compute_N₀(cosθ::T,e²::T,ϵ::T)::T where T
  local cosθ² = cosθ * cosθ
  local R₀ = max(T(1) - e² * cosθ²,ϵ)
  return T(1)/√(R₀)
end
@inline _increment_penality(ρ::T, γ::T, ρ_max::T) = min(ρ*γ,ρ_max)
@inline _distance²(cx::T,cy::T) = cx * cx + cy * cy
@inline _compute_minimum_t(origin_time_direction_i::T, dx_i::T, dy_j::T, pointx_j::T, pointy_j::T) = -origin_time_direction_i + (dx_i * pointx_j + dy_j * pointy_j)
@inline _compute_ray_at__t(point_i, direction_i, t_i) = point_i + t_i * direction_i
@inline function _top_skip_condition(i::Int, j::Int,iter::Int)
  iter==1 &&    return true # skip the first iteration
  j>0 || i>0 || return false # if any of the index is negative, return false
  return true
end

@inline function _update_index_i(i::Int, N::Int, isPeriodic:: Bool)
  if isPeriodic
    return mod1(i,N)
  else
    return update_index_j(i,N)
  end
end
@inline _update_index_j(j::Int, M::Int)= j<=M ? j : -1

const Δθ_SHIFT = 1e-5

"""
  fast_ray_tracing!(t_out::A,θ_out::A,s_out::A,apx::A,apy::A,adx::A,ady::A,incident_refractive_index::A,aθmin::A,aθmax::A,ascending::AB,
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
function fast_ray_tracing!(t_out::A,θ_out::A,s_out::A,apx::A,apy::A,adx::A,ady::A,incident_refractive_index::A,aθmin::A,aθmax::A,ascending::AB,
  atm_n::M,atm_θ::V1,atm_h::V2,retrieval_i::RETi,retrieval_j::RETi,retrieval_n::RETf,retrieval_θ::RETf,retrieval_t::RETf,retrieval_h::RETf,
  retrieval_px::RETf,retrieval_py::RETf,retrieval_dx::RETf,retrieval_dy::RETf,
  tangent_quote::A;δ=1e-10,kmax::Int=30,rho_max::T=T(50.0),gamma::(T),
  intersection_max::Int=140,initialized=false,free_space::T=T(1), kwargs...
  ) where {V1<:AbstractVector, V2<:AbstractVector, RETi, RETf, M, A<:AbstractArray{T}, AB<:AbstractArray{Bool}} where T<:IEEEFloat
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

  @info "Starting the ray tracing"
  @info "max_altitude is $max_altitude"
  @info " extrema of atm_h is $(minimum(atm_h)) and $(maximum(atm_h))"

  @inbounds for iter in 1:iter_eff

      #@batch for idx_rays in eachindex(t_out)
      for idx_rays in eachindex(t_out)
        @info "Ray $idx_rays and iteration $iter"

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

        ## Earth Ray(s)
        local nₜ  = T(0)

        # function to minimize and its derivatives
        local f = T(0)
        local fold = T(0)


        ###########################################################################################
        # First initialization
        # I can assume this has been taken care during the end of the code after the first iteration
        local i_wedge_left= retrieval_i[idx_rays,iter] # gibberish the first iteration
        local j_wedge_top= retrieval_j[idx_rays,iter] # gibberish the first iteration
        ############################################################################################
        local i_wedge_right = _update_index_i(i_wedge_left+1,Natm_n,IsPeriodic)
        local j_wedge_bottom = _update_index_j(j_wedge_top+1,Matm_n)
        local isAscending = ascending[idx_rays]  # initially false
        #################################################
        # Initialize local variables
        # @turbo and @batch needs to know at compilation time
        #################################################

        local s_bottom  = j_wedge_top > 0 ? atm_h[j_wedge_bottom]    : max_altitude
        local s_top     = j_wedge_bottom    > 0 ? atm_h[j_wedge_top] : T(-99)
        local nᵢ        = iter           > 1 ? incident_refractive_index[idx_rays] : free_space
        local s         = iter           > 1 ? s_out[idx_rays]   :  max_altitude
        local θ         = iter           > 1 ? θ_out[idx_rays]   : atan(py1,px1)  # initial guess
        local θmin      = iter           > 1 ? aθmin[idx_rays]   : -Inf
        local θmax      = iter           > 1 ? aθmax[idx_rays]   :  Inf
        local s_target²  = s*s # to be used later
        local current_quote_tangent = iter >1 ? tangent_quote[idx_rays] : T(Inf)
        #################################################
        # DEBUG
        @info "---------------------------------------------"
        @info " try to minimize towards $s "
        @info " top $(s_top) and bottom$(s_bottom))"
        @info " isAscending is $isAscending"
        @info "---------------------------------------------"
        #################################################

        local px1_current = px1
        local py1_current = py1
        local dx1_current = dx1
        local dy1_current = dy1


        # if iter>1 and any index is negative

        if _top_skip_condition(min(i_wedge_left,i_wedge_right),min(j_wedge_top,j_wedge_bottom),iter)
            retrieval_h[idx_rays,iter+1]=T(-999)
            retrieval_t[idx_rays,iter+1]=T(-999)
            retrieval_px[idx_rays,iter+1]=px
            retrieval_py[idx_rays,iter+1]=px
            retrieval_dx[idx_rays,iter+1]=dx
            retrieval_dy[idx_rays,iter+1]=dx
            retrieval_i[idx_rays,iter+1]=i_wedge_left
            retrieval_j[idx_rays,iter+1]=j_wedge_top
            retrieval_n[idx_rays,iter+1]=nᵢ
            retrieval_θ[idx_rays,iter+1]=θ
            continue
        end

        begin

          ### Trugonometric functions
          local cosθ = cos(θ)
          local sinθ = sin(θ)
          local bcosθ = b_normalized * cosθ
          local bsinθ = b_normalized * sinθ


          local N₀ = _compute_N₀(θ,e²,ϵ)

          #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
          local penality_t = T(0)
          local f_with_penality
          local rho = 0.0 #penality factor
          local gamma = 1.5 #penality factor
          #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

          local px2 = cosθ
          local py2 = bsinθ
          local dx2 = bcosθ
          local dy2 = sinθ
              # point on the normal of the ellipse
          local Fx   =  _compute_ray_at__t(px2, dx2*N₀,s)
          local Fy   =  _compute_ray_at__t(py2, dy2*N₀,s)
          t     = _compute_minimum_t(origin_times_direction_1, dx1, dy1, Fx, Fy)


          Px = _compute_ray_at__t(px1, dx1, t)
          Px = _compute_ray_at__t(px1, dx1, t)
          Py = _compute_ray_at__t(py1, dy1, t)

          local fx= Fx - Px
          local fy= Fy - Py
          local k_iter = 0
          f = _distance²(fx,fy)
          fold = f

          # Newton loop

          for k_iter in 1:kmax
            p_newton =  _compute_newtop_step(s, dx1, dy1, cosθ, sinθ, Fx,Fy,Px,Py, b_normalized,e²,N₀ )
            θ = _update_angle(θ,p_newton) mod2pi(θ + p_newton)
            cosθ = cos(θ)
            sinθ = sin(θ)
            N₀ = _compute_N₀(cosθ,e²,ϵ)
            Fx   =  _compute_ray_at__t(px2, dx2*N₀,s)
            Fy   =  _compute_ray_at__t(py2, dy2*N₀,s)
            t     = _compute_minimum_t(origin_times_direction_1, dx1, dy1, Fx, Fy)
            Px    = _compute_ray_at__t(px1, dx1, t)
            Py    = _compute_ray_at__t(py1, dy1, t)
            fx= Fx - Px
            fy= Fy - Py
            f = _distance²(fx,fy)

            if abs(f-fold)<δ
              @info "Convergence reached in $k_iter iterations"
              break
            end
            fold = f
          end

          # compute the distance between the new point Px,Px and the ellipse with origin
          # at px2,py2 and the normal Fx,Fy
          # I know that the minimum distance between the ellipse and the point of interest is  the ⟂ distance
          # so even if it does not converge exactly to the ellipse, I can see if I am above below or stuck


          # I am moving towards the ellipse
          t = max(t,0)
          Px = _compute_ray_at__t(px1, dx1, t)
          Py = _compute_ray_at__t(py1, dy1, t)
          if t==0
            (s>=max_altitude) && j_wedge_top = -1
            (s<=0) && j_wedge_top = 0
            isAscending = !isAscending
            current_quote_tangent=min(s,current_quote_tangent)
          end

          local Δp12x = Px - px2  #earth to ray
          local Δp12y = Py - py2  #earth to ray

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

          @inline function  _interception_two_lines(x1::T,y1::T,x2::T,y2::T,d1::T, d) where T
            return (A-D)/(B-D)



          end


          if  (Δp12² <= s_target²) && (isAscending==false)

          elseif (Δp12² => s_target²) && (isAscending==true)
          end


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

        begin
          if iter==1 && initialized==false
            # find the index of the wedge using binary search
            # assume left orientation of the ray
            # TO DO: add the right orientation of the ray
            i_wedge_left = findlast(atm_θ.<= θ)
            i_wedge_right = i_wedge_left+1
            j_wedge_top = 1
            j_wedge_bottom = 2
            #if s>max_altitude+ϵ  # the first ray never intersected the atmosphere
            #  j_wedge_top = -2
            #end
            # NOTE: the first intersection can only be downward so dx2 and dy2 are already the outward normal
            # to the wedge


          else
            @info "--------------------------------------------------"
            @info "  θ   : $(rad2deg(θ))°  "
            @info "  θmin: $(rad2deg(θmin))°  "
            @info "  θmax: $(rad2deg(θmax))°  "
            @info "  Δh  : $(f)  "
            @info "  s_top: $(s_top)  "
            @info "  s_bottom: $(s_bottom)  "
            @info "  s  : $(s)  "
            @info "--------------------------------------------------"
            begin
              local tmp = dx2
              if  θmin<θ<θmax
                @info "Between 2 wedges"
                if  (abs(f)>10^-5 && isAscending==false)
                    @info " in the middle of the atmosphere"
                    tangent_quote[idx_rays]=s
                    isAscending = true
                    # assume left handiness
                    # TO DO: implement the right handiness
                    # set to the next ray
                    θ=θmax
                    local t_old=t
                    local i_wedge_old=i_wedge_left
                    local j_wedge_old=j_wedge_top
                    # update to the next wedge
                    begin
                      i_wedge_left = i_wedge_left+1
                      i_wedge_right = i_wedge_left+1
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


                      local Δp12² = _distance²(Δp12x,Δp12y)





                      if abs(det) > ϵ
                        t = (Δp12x * dy2 - Δp12y * dx2) / det
                        s = (Δp12x * dy1 - Δp12y * dx1) / det
                      end
                      if (s < 0) # it starts from the ellipse surface, so a negative value is not possible
                        s = NaN
                      end
                    end


                elseif isAscending==true
                  @info " going up"
                  j_wedge_top = j_wedge_top-1  # the direction of h is descending
                  dx2 = -dx2
                  dy2 = -dy2
                else
                  @info " going down"
                  j_wedge_top = j_wedge_top+1  # the direction of h is ascending
                end
              elseif θ==θmax
                @info "Touching left"
                begin
                  i_wedge_left = i_wedge_left+1
                  i_wedge_right = i_wedge_left+1
                  dx2 =  dy2
                  dy2 = -tmp
                end
              elseif θ==θmin
                @info "Touching right"
                begin
                  i_wedge_left = i_wedge_left-1
                  dx2 = -dy2
                  dy2 =  tmp
                end
              # something went wrong
              else
                j_wedge_top = -2
              end
              # check i_index for periodic radial distribution
            end

            # update indexes
            begin
              i_wedge_right = i_wedge_left+1
              j_wedge_bottom = j_wedge_top+1
              if IsPeriodic
                i_wedge_left = mod1(i_wedge_left,Natm_n)
                i_wedge_right = mod1(i_wedge_left+1,Natm_n)
              end
              i_wedge_left = i_wedge_left< N_atmn ? i_wedge_left : -1
              i_wedge_right = i_wedge_right< N_atmn ? i_wedge_right : -1
              j_wedge_top = j_wedge_top< Matm_n ? j_wedge_top : -1
              j_wedge_bottom = j_wedge_bottom< Matm_n ? j_wedge_bottom : -1
            end





          end

          retrieval_i[idx_rays,iter+1]=i_wedge_left
          retrieval_j[idx_rays,iter+1]=j_wedge_top
        end


        # check if atmosphere has been reached or if the ray has failed to intersect
        begin
          # j_wedge_top is used as a flag in the code
          # -2 means an error
          # 0   means it left the atmosphere
          # -1  means it reached the ground
          if j_wedge_top<1 || i_wedge_left<1 || i_wedge_right<1

            number_rays_stopped+=1

            continue
          end
        end
        # update the refractive index of the neighbor
        begin
          nₜ = atm_n[i_wedge_left,j_wedge_top]
          θmin = atm_θ[i_wedge_left]
          θmax = atm_θ[i_wedge_right]
          s_top = atm_h[j_wedge_top]
          s_bottom = atm_h[j_wedge_bottom]
        end
        ##################################################
        # DEBUG
        ##################################################
        @info "--------------------------------------------------"
        @info "Ray $(idx_rays) at $(iter) iteration"
        @info "--------------------------------------------------"
        @info "i_wedge_left=$(i_wedge_left) j_wedge_top=$(j_wedge_top)"
        @info "i_wedge_right=$(i_wedge_right) j_wedge_bottom=$(j_wedge_bottom)"
        @info "nₜ=$(nₜ) nₜ=$(nₜ)"
        @info "θmin=$(θmin) θmax=$(θmax)"
        @info "s_top=$(s_top) s_bottom=$(s_bottom)"
        @info "px1=$(px1) py1=$(py1)"
        @info "px2=$(px2) py2=$(py2)"
        @info "dx1=$(dx1) dy1=$(dy1)"
        @info "dx2=$(dx2) dy2=$(dy2)"
        @info "t=$(t) s=$(s)"
        @info "θ=$(θ)"
        @info "isAscenging? $(isAscending)"
        @info "--------------------------------------------------"
        #return
        ##################################################
        # bending the ray
        if !(nᵢ==nₜ) # do bend only if the refractive index are different
          @info " --------------------------------------------------"
          @info "Bending the ray"
          @info " --------------------------------------------------"
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
          @info " dx_new=$(dx1) dy_new=$(dy1)"

          @info " --------------------------------------------------"


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
          retrieval_i[idx_rays,iter+1]=i_wedge_left
          retrieval_j[idx_rays,iter+1]=j_wedge_top
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
