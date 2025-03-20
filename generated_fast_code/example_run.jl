include("./implementation_new_ray_tracing.jl")
include("./single_functions.jl")
using NCDatasets
using DataFrames
using GeoUtils
using StructArrays
using WGLMakie
function nadir_angle(w,z,ang)
  θ = atan(z/w)
 (nx,ny)=(-w,-z)|> x-> x./hypot(x...) .*-1.0
 #################################
 angle= ang

 dir=_rotation_matrix(angle)*[nx,ny]
 return (dir[1],dir[2])
end
major_axis_earth = majoraxis(ellipsoid(WGS84Latest)) |> x-> uconvert(Unitful.km,x) |> ustrip
minor_axis_earth = minoraxis(ellipsoid(WGS84Latest)) |> x-> uconvert(Unitful.km,x) |> ustrip
eccentricity²_earth= eccentricity²(ellipsoid(WGS84Latest))
cd("generated_fast_code")
##################################################################
#  EARTH MODEL
##################################################################
b_wgs84 = minor_axis_earth/major_axis_earth
e_wgs84 = sqrt(1-b_wgs84^2)
ray(px,py,dx,dy,t)=(px,py).+(t*dx,t*dy)
ellipse(θ,h)=(cos(θ),b_wgs84*sin(θ)).+h.*(b_wgs84*cos(θ),sin(θ))./sqrt(1-e_wgs84*cos(θ)*cos(θ))
##################################################################
cairt_nc=Dataset("../generated_fast_code/cairt_trace_out_nuovo.nc")
cairt_nc.group
cairt_nc_orbit     = cairt_nc.group["orbit"]
cairt_nc_atm       = cairt_nc.group["atm"]
cairt_nc_scans     = cairt_nc.group["scans"]
cairt_nc_integrals = cairt_nc.group["integrals"]

# Lat and Lon
# Info scans
engineering_quotes = cairt_nc_scans["engineering_quotes"][:]
view_angles_deg = cairt_nc_scans["view_angles"][:,:]
tangent_points_z_km  = cairt_nc_scans["tangent_points"][1,:,:]
tangent_points_θ_deg = cairt_nc_scans["tangent_points"][2,:,:]
#######################################################
path_lengths_km          = cairt_nc_integrals["length"][:,:,:] # I assume it to be the path length  equivalent to t in my code
#######################################################
# Information about the atmosphere
#######################################################

phi_cloves=cairt_nc_atm["phi"][:]
z_cloves=cairt_nc_atm["z"][:]

altitudes_km = cairt_nc_atm["altitude"][:,:]

pressure_hPa = cairt_nc_atm["pressure"][:,:]
temperature_K = cairt_nc_atm["temperature"][:,:]

#########################################################
# Normalization
#########################################################
altitudes_km./=major_axis_earth
tangent_points_z_km./=major_axis_earth
path_lengths_km./=major_axis_earth
#########################################################

# get the average temperature and pressure in each clove
size(phi_cloves)
temperature_celsius =similar(temperature_K)

temperature_celsius = uconvert.(u"°C",temperature_K*u"K") |> fn-> @. ustrip(fn)
pressure_Pa = uconvert.(u"Pa",pressure_hPa*u"hPa") |> fn-> @. ustrip(fn)
reverse!(temperature_celsius,dims=1)
reverse!(pressure_Pa,dims=1)
reverse!(z_cloves)
#reverse!(phi_cloves)
z_cloves
@. mod2pi(deg2rad(phi_cloves))
#initialize atmosphere
atmosphere = StructArray(Matrix{AtmosphereProfile2D{Float64}}(undef, size(temperature_K,2),size(temperature_K,1)))

for j in eachindex(z_cloves)
    for i in eachindex(phi_cloves)
        atmosphere[i,j]=AtmosphereProfile2D( 1.0,
        let
          ismissing(temperature_celsius[j,i]) ? NaN : temperature_celsius[j,i]
        end,
        let
          ismissing(pressure_Pa[j,i]) ? NaN : pressure_Pa[j,i]
        end,
        mod2pi(deg2rad(phi_cloves[i])),
        z_cloves[j]./major_axis_earth
        )
    end
end

begin
  # average radially
  atmosphere.temperature_ave[1:end-1,:] = (atmosphere.temperature_ave[1:end-1,:]+atmosphere.temperature_ave[2:end,:])./2
  atmosphere.pressure_ave[1:end-1,:] = (atmosphere.pressure_ave[1:end-1,:]+atmosphere.pressure_ave[2:end,:])./2
  # average altitude
  atmosphere.temperature_ave[:,1:end-1] = (atmosphere.temperature_ave[:,1:end-1]+atmosphere.temperature_ave[:,2:end])./2
  atmosphere.pressure_ave[:,1:end-1] = ( atmosphere.pressure_ave[:,2:end]- atmosphere.pressure_ave[:,1:end-1])./(log.(atmosphere.pressure_ave[:,2:end]./atmosphere.pressure_ave[:,1:end-1]))
  for j in 1:size(atmosphere.temperature_ave,2)
    for i in 1:size(atmosphere.temperature_ave,1)
      atmosphere.refraction_index_ave[i,j]=refractive_index(temperature=atmosphere.temperature_ave[i,j],
        pressure=atmosphere.pressure_ave[i,j],CO2ppm=400)
    end
  end
  # average altitude
end
# get refractive index

@inline extrema_noNaN(x)=extrema(filter(x-> !( ismissing(x) || isnan(x)) ,x))


figure=Figure(size=(800/1.5,800))
with_theme(theme_latexfonts()) do
  ax1=Axis(figure[1,1][1,1],title="refraction_index (n-1)",xlabel="ϕ°",ylabel="h km")
  ax2=Axis(figure[2,1][1,1],title="temperature",xlabel="ϕ°",ylabel="h km")
  ax3=Axis(figure[3,1][1,1],title="pressure",xlabel="ϕ°",ylabel="h km")



  s1=surface!(ax1,atmosphere.θ_left,atmosphere.s_top,
  atmosphere.refraction_index_ave.-1.0)
  s2=surface!(ax2,atmosphere.θ_left,atmosphere.s_top,
  atmosphere.temperature_ave)
  s3=surface!(ax3,atmosphere.θ_left,atmosphere.s_top,
  atmosphere.pressure_ave./extrema_noNaN(atmosphere.pressure_ave)[2])
  Colorbar(figure[1,1][1,2],s1)
  Colorbar(figure[2,1][1,2],s2)
  Colorbar(figure[3,1][1,2],s3)
end
figure
save("./refraction_index_temperature_pressure.png",figure)

# Using equation from part 3 of Casia DEL006
r_satellite  = cairt_nc_orbit["radius"][1]/major_axis_earth
inclination  = cairt_nc_orbit["inclination"][1]

# angle of the current position in degrees
scans_phi = cairt_nc_scans["scans_phi"][:]
satx = @. r_satellite*cosd(scans_phi)
saty = @. r_satellite*sind(scans_phi)

total_scans = prod(size(tangent_points_z_km))
max_iterations = 140

##########################################################
#  INPUT ARRAYS
##########################################################
# Retrieval array
retrieval = StructArray(Matrix{ResultsRayTracing{Float64}}(undef, total_scans,max_iterations+1))
# Input array
inputray  = StructArray(Vector{InputRay{Float64}}(undef, total_scans))
θsat,hsat =  geocentric_xy_to_geodesic_θ(satx[2],saty[2])

# 1. Set the initial conditions
fov=5  # number of field of vies per scan
for i in eachindex(satx,saty)
  local θsat,hsat =  geocentric_xy_to_geodesic_θ(satx[i],saty[i])
  @info "satellite position" θsat,hsat
  local N = 1/sqrt(cosd(θsat)^2+b_wgs84^2*sind(θsat)^2)
  local sx1,sy1 = (N*cosd(θsat),b_wgs84^2*N*sind(θsat))
  local nx1,ny1 = (sx1,sy1/b_wgs84^2)
  local norm1 = hypot(nx1,ny1)
  nx1/=norm1
  ny1/=norm1
  @info " px: $(satx[i]), py: $(saty[i])"
  @info " px: $(sx1+hsat*nx1), px: $(sy1+hsat*ny1)"
  @info "------------------------------"
  @info " nx: $nx1, ny: $ny1"

  for j in 1:fov
    inputray.px[(i-1)*fov+j]=satx[i]
    inputray.py[(i-1)*fov+j]=saty[i]
    local angle= -view_angles_deg[(i-1)*fov+j]
    scan_direction=nadir_angle_normal(nx1,ny1,angle;outward=false)
    @info scan_direction
    inputray.dx[(i-1)*fov+j]=scan_direction[1]
    inputray.dy[(i-1)*fov+j]=scan_direction[2]
  end
end

marco_ray_tracing!(t_out::A,θ_out::A,s_out::A,
  inputray::IR,outputray::OR,atmosphere::ATM,
  tangent_quote::A;kwargs...) where {IR<:AbstractVector{InputRay{T}},
  OR<:AbstractMatrix{ResultsRayTracing{T}},
  ATM<:AbstractMatrix{AtmosphereProfile2D{T}},
  A<:AbstractVector{T}} where T=
marco_ray_tracing!(t_out,θ_out,s_out,
   inputray.px,inputray.py,
   inputray.dx,inputray.dy,
  inputray.n,inputray.θmin,inputray.θmax,inputray.ascending,
  atmosphere.refraction_index_ave[1:end-1,1:end-1],atmosphere.θ_left[1:end,1],atmosphere.s_top[1,1:end],
  outputray.i,outputray.j,outputray.n,outputray.θ,outputray.t,outputray.h,
  outputray.px,outputray.py,outputray.dx,outputray.dy,
  tangent_quote,kwargs...)

start_scan = 48
end_scan = 50
num_scans = 5
begin
test_retrieval=deepcopy(retrieval[start_scan:end_scan,1:num_scans+1])
test_inputray=deepcopy(inputray[start_scan:end_scan])
test_atmosphere=deepcopy(atmosphere)
tangent_quote =Array{Float64}(undef, size(test_inputray.px))
t_out= similar(tangent_quote)
θ_out= similar(tangent_quote)
s_out= similar(tangent_quote)
end


tangent_quote


#using BenchmarkTools
#@benchmark marco_ray_tracing!($t_out,$θ_out,$s_out, $test_inputray,$test_retrieval,$test_atmosphere,$tangent_quote)

marco_ray_tracing!(t_out,θ_out,s_out, test_inputray,test_retrieval,test_atmosphere,tangent_quote)


for j in 1:10
  fig=Figure(size=(600,1000))
  with_theme(theme_latexfonts()) do

    jjj=60
    ax1=Axis(fig[1,1][1,1],title="refractive index in clove",xlabel="iteration",ylabel="n-1")
    ax2=Axis(fig[2,1][1,1],title="length ray",xlabel="iteration",ylabel="km")
    ax3=Axis(fig[3,1][1,1],title="orbital coordinate",xlabel="iteration",ylabel="θ°")
    ax4=Axis(fig[4,1][1,1],title="height from surface",xlabel="iteration",ylabel="altitude km")


    for i in 1:5
    lines!(ax1,1:jjj-2,[(r-1) for r in test_retrieval.n[i+fov*(j-1),3:jjj]],label="LOS $i",linewidth=3)
      lines!(ax2,1:jjj-2,[r.*major_axis_earth for r in test_retrieval.t[i+fov*(j-1),3:jjj]],label="LOS $i",linewidth=3)
      lines!(ax3,1:jjj-2,[rad2deg(r) for r in test_retrieval.θ[i+fov*(j-1),3:jjj]],label="LOS $i",linewidth=3)
      lines!(ax4,1:jjj-2,[r.*major_axis_earth for r in test_retrieval.h[i+fov*(j-1),3:jjj]],label="LOS $i",linewidth=3)
      xlims!(ax1,1,jjj-2)
      xlims!(ax2,1,jjj-2)
      xlims!(ax3,1,jjj-2)
      xlims!(ax4,1,jjj-2)
    end
    Legend(fig[1,1][1,2],ax1)
    Legend(fig[2,1][1,2],ax2)
    Legend(fig[3,1][1,2],ax3)
    Legend(fig[4,1][1,2],ax4)

    Label(fig[0, :][1,1], "Scan $j";tellwidth=false, fontsize=30)

  end
  figure
  save("./scan$j.png",fig)
end
figure
reshape(tangent_quote,5,10).*major_axis_earth

test_retrieval.h[:,73:76].*major_axis_earth
reshape(tangent_quote .*major_axis_earth,5,10)
engineering_quotes
tangent_points_z_km' *major_axis_earth


fig=Figure()
ax=Axis(fig[1,1])
ww=48
bb=3
scatterlines!(ax,[(x,y) for (x,y) in zip(retrieval.px[ww,1:bb],retrieval.py)],label="engineering_quotes")
fig




figure=Figure()
ax=Axis(figure[1,1])

for inp in inputray
  lines!(ax,
  [ray(inp.px,inp.py,inp.dx,inp.dy,t) for t in (0,1)])
end
ext_theta=extrema(atmosphere.θ_left)
ext_h    =extrema(atmosphere.s_top)
θ_radii = atmosphere.θ_left[:,1]
h_levels= atmosphere.s_top[1,:]
for θ in θ_radii
  lines!(ax,[ellipse(θ,h) for h in ext_h],color=:black )
end
for h in h_levels
  lines!(ax,[ellipse(θ,h) for θ in LinRange(ext_theta...,1000)],color=:black )
end



scatter!(ax,[(px,py) for (px,py) in zip(test_retrieval.px[1,2:end],test_retrieval.py[1,2:end])])
TTT=10


scatter!(ax,[(px,py) for (px,py) in zip(test_retrieval.px[1,2:TTT],test_retrieval.dy[1,2:TTT])])

for (px,py,dx,dy,t) in zip(test_retrieval.px[1,2:end],test_retrieval.py[1,2:end],test_retrieval.dx[1,2:end],test_retrieval.dy[1,2:end],test_retrieval.t[1,2:end])
  lines!(ax,[ray(px,py,dx,dy,t) for t in (0,1)])
end

test_retrieval.dx'
