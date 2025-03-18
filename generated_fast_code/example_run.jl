include("../generated_fast_code/implementation_new_ray_tracing.jl")
using NCDatasets
using DataFrames
using GeoUtils
using StructArrays

major_axis_earth = majoraxis(ellipsoid(WGS84Latest)) |> x-> uconvert(Unitful.km,x) |> ustrip
minor_axis_earth = minoraxis(ellipsoid(WGS84Latest)) |> x-> uconvert(Unitful.km,x) |> ustrip
eccentricity²_earth= eccentricity²(ellipsoid(WGS84Latest))
cd("generated_fast_code")
##################################################################
#  EARTH MODEL
##################################################################
b_wgs84 = minor_axis_earth/major_axis_earth
set_minoraxis(b_wgs84)
b_mine  = get_minoraxis()
e²_mine = get_e²()
eccentricity²_earth≈get_e²()
b_wgs84==b_mine
##################################################################

cairt_nc=Dataset("../generated_fast_code/cairt_trace_out_test.nc")
cairt_nc.group
cairt_nc_orbit     = cairt_nc.group["orbit"]
cairt_nc_atm       = cairt_nc.group["atm"]
cairt_nc_scans     = cairt_nc.group["scans"]
cairt_nc_integrals = cairt_nc.group["integrals"]

# Lat and Lon
lat_degrees=cairt_nc_atm["latitude"][:]
lon_degrees=cairt_nc_atm["longitude"][:]

# Info scans
polar_coordinates_along_orbit_deg = cairt_nc_scans["scans_phi"][:]
# Don't know what is geometry
# engineering quotes in km
engineering_quotes = cairt_nc_scans["engineering_quotes"][:]
ccd_angles_deg = cairt_nc_scans["ccd_angles"][:]
view_angles_deg = cairt_nc_scans["view_angles"][:,:]
tangent_points_z_km  = cairt_nc_scans["tangent_points"][1,:,:]
tangent_points_θ_deg = cairt_nc_scans["tangent_points"][2,:,:]
#######################################################
cairt_nc_orbit
#######################################################
npaths_scans = cairt_nc_integrals["npaths"][:,:]  # I assume it to be the number of interceptions

phi_index_cloves_per_los = cairt_nc_integrals["phi_idx"][:,:,:] # I assume it to be the index of the phi clove
z_index_per_los          = cairt_nc_integrals["z_idx"][:,:,:] # I assume it to be the index of the z clove
path_lengths_km          = cairt_nc_integrals["length"][:,:,:] # I assume it to be the path length  equivalent to t in my code

cairt_nc_integrals
#######################################################
# Information about the atmosphere
#######################################################

phi_cloves=cairt_nc_atm["phi"][:]
z_cloves=cairt_nc_atm["z"][:]

ngrids= cairt_nc_atm["ngrid"][:]

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
extrema(ngrids)
extrema(phi_index_cloves_per_los)
extrema(z_index_per_los)  # i-1
# point of entrance is largest clove index (h-1)
min_altitude=extrema(abs.(z_index_per_los[:,:,:]))[1]
max_altitude=extrema(z_index_per_los[:,:,:])[2]

# get the average temperature and pressure in each clove
size(phi_cloves)
temperature_celsius =similar(temperature_K)

temperature_celsius = uconvert.(u"°C",temperature_K*u"K") |> fn-> @. ustrip(fn)
pressure_Pa = uconvert.(u"Pa",pressure_hPa*u"hPa") |> fn-> @. ustrip(fn)


#initialize atmosphere
atmosphere = StructArray(Matrix{AtmosphereProfile2D{Float64}}(undef, size(temperature_K,2),size(temperature_K,1)))

for j in eachindex(z_cloves)
    for i in eachindex(phi_cloves)
        atmosphere[i,j]=AtmosphereProfile2D( 1.0,
        let
          ismissing(temperature_K[end-j+1,i]) ? NaN : temperature_celsius[end-j+1,i]
        end,
        let
          ismissing(pressure_hPa[end-j+1,i]) ? NaN : pressure_Pa[end-j+1,i]
        end,
        deg2rad(lat_degrees[i]),
        z_cloves[end-j+1]./major_axis_earth
        )
    end
end

z_cloves

atm_h

begin
  # average radially
  atmosphere.temperature_ave[1:end-1,:] = (atmosphere.temperature_ave[1:end-1,:]+atmosphere.temperature_ave[2:end,:])./2
  atmosphere.pressure_ave[1:end-1,:] = (atmosphere.pressure_ave[1:end-1,:]+atmosphere.pressure_ave[2:end,:])./2
  # average altitude
  atmosphere.temperature_ave[:,1:end-1] = (atmosphere.temperature_ave[:,1:end-1]+atmosphere.temperature_ave[:,2:end])./2
  atmosphere.pressure_ave[:,1:end-1] = ( atmosphere.pressure_ave[:,1:end-1]- atmosphere.pressure_ave[:,2:end])./(log.(atmosphere.pressure_ave[:,1:end-1]./atmosphere.pressure_ave[:,2:end]))
  for j in 1:size(atmosphere.temperature_ave,2)
    for i in 1:size(atmosphere.temperature_ave,1)
      atmosphere.refraction_index_ave[i,j]=refractive_index(temperature=atmosphere.temperature_ave[i,j],
        pressure=atmosphere.pressure_ave[i,j],CO2ppm=400)
    end
  end
  # average altitude
end
# get refractive index



figure=Figure(size=(800/2,800))
begin
ax1=Axis(figure[1,1][1,1],title="refraction_index (n-1)",xlabel="ϕ°",ylabel="h km")
ax2=Axis(figure[2,1][1,1],title="temperature",xlabel="ϕ°",ylabel="h km")
ax3=Axis(figure[3,1][1,1],title="pressure",xlabel="ϕ°",ylabel="h km")

s1=surface!(ax1,phi_cloves,z_cloves,atmosphere.refraction_index_ave.-1.0)
s2=surface!(ax2,phi_cloves,z_cloves,atmosphere.temperature_ave)
s3=surface!(ax3,phi_cloves,z_cloves,atmosphere.pressure_ave)
Colorbar(figure[1,1][1,2],s1)
Colorbar(figure[2,1][1,2],s2)
Colorbar(figure[3,1][1,2],s3)
end
figure
save("./refraction_index_temperature_pressure.png",figure)


# Using equation from part 3 of Casia DEL006
r_satellite  = cairt_nc_orbit["radius"][1]/major_axis_earth
inclination  = cairt_nc_orbit["inclination"][1]

inclination  = cairt_nc_orbit["inclination"]
cairt_nc_orbit
cairt_nc_atm

lon0= cairt_nc_orbit["lon0"][1]
# reduced coordinates
px = r_satellite.*cosd(inclination)
py = r_satellite.*sind(inclination)
#px = r_satellite.*cosd(lon0)
#py = r_satellite.*sind(lon0)

cairt_nc_atm
total_scans = prod(size(tangent_points_z_km))
max_iterations = 140

##########################################################
#  INPUT ARRAYS
##########################################################
# Retrieval array
retrieval = StructArray(Matrix{ResultsRayTracing{Float64}}(undef, total_scans,max_iterations+1))
# Input array
inputray  = StructArray(Vector{InputRay{Float64}}(undef, total_scans))

# need to compute the normal to the surface I can use the gradient to the surface with
# Find the x,y given the θ_geodesic
# 1. find angle of the satellite in geodesic coordinates
#   1.1 Initial guess geodesic angle
    (θ0_geodesic,h_sat)= geocentric_xy_to_geodesic_θ(px,py)
#   1.2 correct the angle
    θ0_geodesic_rad=deg2rad(θ0_geodesic)
    θ_geodesic_rad = initialize_theta(θ0_geodesic_rad,px,py)
# 2. Compute the cartesian cordinates
    earth_x,earth_y = (cos(θ_geodesic_rad),sin(θ_geodesic_rad)*b_mine)
# 3. Get the normal to the surface
    normal_earth_x,normal_earth_y = earth_x,earth_y/b_mine^2
    (normal_earth_x,normal_earth_y) = (normal_earth_x,normal_earth_y)./hypot(normal_earth_x,normal_earth_y)
# 4. Confirm the result
    p′_x,p′_y = (earth_x,earth_y) .+ (normal_earth_x,normal_earth_y).*h_sat
    @assert hypot(p′_x-px,p′_y-py)<1e10  # true
# 5. Set the inputray with the origin of the satellite
    inputray.px[:].=px
    inputray.py[:].=py
# 6. Compute the direction of the scan angles
# The scan is give in the geodesic coordinates AND as a nadir angle
# so I need to use the nadir angle normal to compute the direction
    for i in eachindex(view_angles_deg)
      local angle= view_angles_deg[i]
      scan_direction=nadir_angle_normal(normal_earth_y,normal_earth_x,angle;outward=true)
      inputray.dx[i]=scan_direction[1]
      inputray.dy[i]=scan_direction[2]
    end
# 7. Visually check if the direction makes sense

    # SHIFT THE ORIGIN OF THE GRID BY THE ANGLE OF THE SATELLITE
    atmosphere.θ_left[:,:]=repeat(deg2rad.(lat_degrees[:]),1,size(atmosphere,2)).+θ_geodesic_rad

    figure=Figure()
    ax=Axis(figure[1,1][1,1])
    let
      hh= extrema(atmosphere.s_top)[2]
      atmosphere.θ_left[:,1]
      for θ in atmosphere.θ_left[:,1]
        lines!(ax,[ellipse(θ,h) for h in (0,hh) ],color=:black)
      end
      (thetamin,thetamax)= extrema(atmosphere.θ_left)

      for h in atmosphere.s_top[1,:]
        lines!(ax,[ellipse(θ,h) for θ in LinRange(thetamin,thetamax,1000) ],color=:black)
      end
      xlims!(ax,0.05,0.55)
      ylims!(ax,0.85,1.018)
      for i  in 1:total_scans
         _px=inputray.px[i]
         _py=inputray.py[i]
         _dx=inputray.dx[i]
         _dy=inputray.dy[i]
       lines!(ax,[(_px,_py).+(_dx,_dy).*t for t in (0,5)] ,color=:red)
      end
    end

    # Atmosphere array
atmosphere = atmosphere
# tangent quote
tangent_quote =Array{Float64}(undef, total_scans)
t_out= similar(tangent_quote)
θ_out= similar(tangent_quote)
s_out= similar(tangent_quote)

fast_ray_tracing!(t_out::A,θ_out::A,s_out::A,
  inputray::IR,outputray::OR,atmosphere::ATM,
  tangent_quote::A;kwargs...) where {IR<:AbstractVector{InputRay{T}},
  OR<:AbstractMatrix{ResultsRayTracing{T}},
  ATM<:AbstractMatrix{AtmosphereProfile2D{T}},
  A<:AbstractVector{T}} where T=
fast_ray_tracing!(t_out,θ_out,s_out,
   inputray.px,inputray.py,
   inputray.dx,inputray.dy,
  inputray.n,inputray.θmin,inputray.θmax,inputray.ascending,
  atmosphere.refraction_index_ave[1:end-1,1:end-1],atmosphere.θ_left[1:end,1],atmosphere.s_top[1,1:end],
  outputray.i,outputray.j,outputray.n,outputray.θ,outputray.t,outputray.h,
  outputray.px,outputray.py,outputray.dx,outputray.dy,
  tangent_quote,kwargs...)

scan_testing = 1
num_scans = 30
begin
test_retrieval=deepcopy(retrieval[1:1,1:num_scans+1])
test_inputray=deepcopy(inputray[1:1])
test_atmosphere=deepcopy(atmosphere)
tangent_quote =Array{Float64}(undef, scan_testing)
t_out= similar(tangent_quote)
θ_out= similar(tangent_quote)
s_out= similar(tangent_quote)
end




fast_ray_tracing!(t_out,θ_out,s_out, test_inputray,test_retrieval,test_atmosphere,tangent_quote)
t_out
test_retrieval[1,2].n
test_retrieval[1,3].n

test_retrieval[1,:].py
tangent_quote
scatter!(ax,test_retrieval[1,:].px,test_retrieval[1,:].py,color=:black)

inputray.n
tangent_quote
atmosphere.refraction_index,ave[54,1]

reft_out
tangent_quote

tangent_quote


atmosphere
typeof(atmosphere)
typeof(inputray.px)

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
end

fast_ray_tracing!(t_out,θ_out,s_out,
apx,apy,adx,ady,
incident_refractive_index,
aθmin,aθmax,ascending,
atm_n,atm_θ,atm_h,
retrieval_i,retrieval_j,
retrieval_n,retrieval_θ,
retrieval_t,retrieval_h,
retrieval_px,retrieval_py,
retrieval_dx,retrieval_dy,
tangent_quote)





fig = Figure()
ax  = Axis(fig[1,1][1,1])
lines!(ax,[(cosd(θ),b_mine*sind(θ)) for θ in range(0,stop=365,length=100)],color=:black)
scatter!(ax,px,py,markersize=51)
