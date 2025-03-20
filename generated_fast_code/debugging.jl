using Base: IEEEFloat
include("implementation_new_ray_tracing.jl")
include("single_functions.jl")
# Test Ray Traced
nrays= 1
niterations=5

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


limb_angle(5,2,65)
θ_limb=70
begin
  θ_out = zeros(Float64,nrays)
  t_out = similar(θ_out)
  s_out = similar(θ_out)
  aθmin = similar(θ_out)
  aθmax = similar(θ_out)
  apx   = ones(Float64,nrays).*1.8
  apy   = ones(Float64,nrays).*2
  if nrays>2
    dir_ang=[limb_angle(5,2,θ_limb+deg2rad(shift)) for shift in LinRange(-5,5,nrays)]
  else
    dir_ang=[limb_angle(5,2,θ_limb)]
  end
  adx   = [x for (x,_) in dir_ang]
  ady   = [y for (_,y) in dir_ang]

  incident_refractive_index = ones(Float64,nrays)
  adx
  ascending = [false for _ in 1:nrays]
  incident_refractive_index = ones(Float64,nrays)
  retrieval_i = zeros(Int,nrays,niterations+1)
  retrieval_j = similar(retrieval_i)
  retrieval_θ = zeros(Float64,nrays,niterations+1)
  retrieval_t = similar(retrieval_θ)
  retrieval_h = similar(retrieval_θ)
  retrieval_n = similar(retrieval_θ)
  retrieval_px = similar(retrieval_θ)
  retrieval_py = similar(retrieval_θ)
  retrieval_dx = similar(retrieval_θ)
  retrieval_dy = similar(retrieval_θ)
  tangent_quote = similar(θ_out)
  M=20
  N=180
  N_atmn= N
  M_atmn= M-1
  atm_h = [ exp(-x) for x in LinRange(0,3,M)]
  atm_θ = [ θ for θ in LinRange(0,2π,N+1)][1:end-1]
  atm_n = ones(Float64,N_atmn,M_atmn)
  M_atmn
  n_horizontal = 1.0.+0.00027.*(1.0.-atm_h[1:end-1])
  n_vertical = @. sin(atm_θ[1:end])*0.0000
  atm_n[:,:]=repeat(n_horizontal',N_atmn,1)
  atm_n[:,:]+=repeat(n_vertical,1,M_atmn)
end;

marco_ray_tracing!(t_out,θ_out,s_out,
  apx,apy,adx,ady,
  incident_refractive_index,
  aθmin,aθmax,ascending,
  atm_n,atm_θ,atm_h,
  retrieval_i,retrieval_j,
  retrieval_n,retrieval_θ,
  retrieval_t,retrieval_h,
  retrieval_px,retrieval_py,
  retrieval_dx,retrieval_dy,
  tangent_quote);

  _=fast_ray_tracing!(t_out,θ_out,s_out,
  apx,apy,adx,ady,
  incident_refractive_index,
  aθmin,aθmax,ascending,
  atm_n,atm_θ,atm_h,
  retrieval_i,retrieval_j,
  retrieval_n,retrieval_θ,
  retrieval_t,retrieval_h,
  retrieval_px,retrieval_py,
  retrieval_dx,retrieval_dy,
  tangent_quote);


  @benchmark fast_ray_tracing!($t_out,$θ_out,$s_out,
  $apx,$apy,$adx,$ady,
  $incident_refractive_index,
  $aθmin,$aθmax,$ascending,
  $atm_n,$atm_θ,$atm_h,
  $retrieval_i,$retrieval_j,
  $retrieval_n,$retrieval_θ,
  $retrieval_t,$retrieval_h,
  $retrieval_px,$retrieval_py,
  $retrieval_dx,$retrieval_dy,
  $tangent_quote)

begin
  b=get_minoraxis()
  e²=1-b*b
  ray(px,py,dx,dy,t)=(px,py).+(t*dx,t*dy)
  ellipse(θ,h)=(cos(θ),b*sin(θ)).+h.*(b*cos(θ),sin(θ))./sqrt(1-e²*cos(θ)*cos(θ))
  figure=Figure()
  ax=Axis(figure[1,1])

  for (px_step,py_step,dx_step,dy_step,t_step) in zip(retrieval_px[1:10,1:end-1],retrieval_py[1:10,1:end-1],retrieval_dx[1:10,1:end-1],retrieval_dy[1:end-1],  retrieval_t[1:10,2:end])
    tpass=0
    for (px,py,dx,dy,t) in zip(px_step,py_step,dx_step,dy_step,t_step)
        scatterlines!(ax,[ray(px,py,dx,dy,t) for t in (0,t)])
    end
  end
  lines!(ax,[ellipse(θ,h) for θ in LinRange(0,2pi+0.001,1000), h in atm_h][:],color="red")
  for θ in atm_θ
    lines!(ax,[ellipse(θ,h) for h in (0,1)][:],color="red")
  end
end
figure



pre_t= 1.8622978896214737
pre_s= 1.0
pre_θ= 1.6024482877051827
pre_px1= 1.8
pre_py1= 2.0
pre_px2= -0.031646676106628795
pre_py2= 0.996147987200637
pre_dx1= -0.9995057873530447
pre_dy1= -0.03143534710751977
pre_dx2= -0.03154067652592609
pre_dy2= 0.9995024690936419
post_t= 1.8623975572325817
post_s= 0.9457774526487894
post_θ= 1.6024482877051827
post_px1= 1.8
post_py1= 2.0
post_px2= -0.031646676106628795
post_py2= 0.996147987200637
post_dx1= -0.9995057873530447
post_dy1= -0.03143534710751977
post_dx2= -0.03154067652592609
dy2= 0.9995024690936419

figure=Figure()
ax=Axis(figure[1,1])

scatterlines!(ax,[(px,py) for (px,py) in zip(retrieval_px[2:end],retrieval_py[2:end])])
tangent_quote
retrieval_px
retrieval_dx
