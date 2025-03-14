using Base: IEEEFloat

const DEFAULT_MINORAXIS= 0.5
const DEFAULT_MINORAXIS²=DEFAULT_MINORAXIS^2
const DEFAULT_E²=1-DEFAULT_MINORAXIS²

# Initialize the minor axis and its square by type
const MINORAXIS_F16=Ref{Float16}(DEFAULT_MINORAXIS)
const MINORAXIS²_F16=Ref{Float16}(DEFAULT_MINORAXIS²)
const E²_F16=Ref{Float16}(DEFAULT_E²)

const MINORAXIS_F32=Ref{Float32}(DEFAULT_MINORAXIS)
const MINORAXIS²_F32=Ref{Float32}(DEFAULT_MINORAXIS²)
const E²_F32=Ref{Float32}(DEFAULT_E²)

const MINORAXIS_F64=Ref{Float64}(DEFAULT_MINORAXIS)
const MINORAXIS²_F64=Ref{Float64}(DEFAULT_MINORAXIS²)
const E²_F64=Ref{Float64}(DEFAULT_E²)


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


for infoconstant in (:minoraxis,:minoraxis²,:e²)
  lower_const_name=String(infoconstant)
  upper_const_name=uppercase(lower_const_name)
  upper_const_symbol_F16=Symbol("$(upper_const_name)_F16")
  upper_const_symbol_F32=Symbol("$(upper_const_name)_F32")
  upper_const_symbol_F64=Symbol("$(upper_const_name)_F64")
  @show upper_const_symbol_F16
  fn=Symbol("get_$lower_const_name")

  @eval @generated function $fn(::Type{T}=Float64) where T
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


# Define the function to calculate the distance between a ray and the ellipsoid with major axis a=1

get_minoraxis(Float32)
set_minoraxis(0.3)





function min_distance_from_ellipse²(θ::T,s::T,Px0::T,Py0::T,Px1::T,Py1::T;δ=1e-9,kmax=10) where T

  origin_times_direction=Px0*Px1+Py0*Py1

  @inline function tmin(Fx,Fy)
    t = (Px1*Fx+Py1*Fy)-(origin_times_direction)
  end
  @inline function derivative_t_wrt_θ(dFx,dFy)
    dtdθ =Px1*dFx+Py1*dFy
  end


  @inline function point_ray(t)
    Px=(Px0+Px1*t)
    Py=(Py0+Py1*t)
    return (Px,Py)
  end

  @inline function fgmin(t,Fx,Fy,dFx,dFy,d²Fx,d²Fy)
    (Px,Py)=point_ray(t)
    Dx=Fx-Px
    Dy=Fy-Py

    # now t(θ) is a function of θ so we need to compute
    dt=derivative_t_wrt_θ(d²Fx,d²Fy)
    d²t=derivative_t_wrt_θ(d²Fx,d²Fy)

    dPx=Px1*dt
    dPy=Py1*dt

    d²Px=Px1*d²t
    d²Py=Py1*d²t



    dDx=dFx-dPx
    dDy=dFy-dPy


    dDx²=dDx*dDx
    dDy²=dDy*dDy

    d²Dx=d²Fx-d²Px
    d²Dy=d²Fy-d²Py

    g=Dx*dDx+Dy*dDy

    h=dDx²+dDy²+Dx*d²Dx+Dy*d²Dy

    h=abs(h)>eps(T) ? h : T(1)
    p=g/h

    f=Dx*Dx+Dy*Dy
    return f,p
  end

  @inline function precomputed_trigonometric_functions(θ)
    cosθ        = cos(θ)
    sinθ        = sin(θ)
    bcosθ       = get_minoraxis(T)*cosθ
    bsinθ       = get_minoraxis(T)*sinθ
    cosθ²       = cosθ*cosθ
    sinθ²       = sinθ*sinθ
    half_sin2θ  = sinθ*cosθ
    cos2θ       = cosθ²-sinθ²

    R=1-get_e²(T)*cosθ²

    N = 1/sqrt(R)
    N³=N*N*N
    N⁵=N³*N*N

    half_dR  = get_e²(T)*half_sin2θ
    half_dR² = half_dR*half_dR
    half_d²R = get_e²(T)*cos2θ
    dN=-N³*half_dR

    (Fx0,Fy0)=(cosθ,bsinθ)
    (Fx1,Fy1)=(bcosθ*R,sinθ*R)

    Fx=Fx0+s*Fx1
    Fy=Fy0+s*Fy1

    dFx0= -sinθ
    dFy0= bcosθ

    dFx1= -bsinθ
    dFy1=  cosθ

    dFx= dFx0+s*(dFx1*N+dN*Fx1)
    dFy= dFy0+s*(dFy1*N+dN*Fy1)


    d²Fx0=-Fx0
    d²Fy0=-Fy0

    d²Fx1= -Fx1
    d²Fy1= -Fy1

    d²N= -N³*half_d²R+3*N⁵*half_dR²

    d²Fx1= d²Fx1*N+2*dN*dFx1+d²N*Fx1
    d²Fy1= d²Fy1*N+2*dN*dFy1+d²N*Fy1


    d²Fx=d²Fx0+d²Fx1
    d²Fy=d²Fy0+d²Fy1

    return (Fx,Fy,
            dFx,dFy,
            d²Fx,d²Fy
            )
  end

  @inline function fgh!(θ)
    (Fx,Fy,dFx,dFy,d²Fx,d²Fy)=precomputed_trigonometric_functions(θ)
    t=tmin(Fx,Fy)
    (fold,dfdθ)=fgmin(t,Fx,Fy,dFx,dFy,d²Fx,d²Fy)
    return fold,dfdθ,t
  end
  #=
      while true
        k+=1

        if abs(fold-f)<δ
          break
        end

        if k>kmax
          @warn "The iteration did not converge"
          break
        end
  =#

    return fgh!(θ)
end

min_distance_from_ellipse²(θ,s,x0,y0,x1,y1)=_distance_from_ellipse²(promote(θ,s,x0,y0,x1,y1))



ang=LinRange(0,4pi,100)
using BenchmarkTools
using StructArrays
x0=pi
x=float(x0)
t=0.
for i in 1:100
  (f,p,t)=min_distance_from_ellipse²(x,sqrt(3).,2.,2.,-cos(pi/4),-sin(pi/4))
  x=-p
  @show i,t,f,p,x
end
ray(t)=(2.,2.).+(-cos(pi/4),-sin(pi/4)).*t

tt=LinRange(0,t,2)

lines!(ray.(tt))

_N1(θ)=1/sqrt(1-get_e²()*cos(0)^2)

lines!([((_N1(θ)+3)*cos(θ),(get_minoraxis²()*_N1(θ)+3)*sin(θ)) for θ in ang])
