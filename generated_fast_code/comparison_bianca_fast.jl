include("./implementation_new_ray_tracing.jl")
using GeoUtils
using StructArrays

using Unitful:°,hPa
const MODEL=Ref{AirModel}(Carlotti())
const INTERPOLATION=Ref{AbstractPressureInterpolation}(LinearPressure())

setModel(A::AirModel)= MODEL[]=A
getModel()=MODEL[]
setModel()=MODEL[]=(Ciddor())
setPressureInterpolation(A::AbstractPressureInterpolation)= INTERPOLATION[]=A
setPressureInterpolation()= INTERPOLATION[]=(LinearPressure())
getPressureInterpolation()=INTERPOLATION[]
getPressureInterpolation()
setModel()
include("./initialization_script_fast.jl")
major_axis_earth = majoraxis(ellipsoid(WGS84Latest)) |> x-> uconvert(Unitful.km,x) |> ustrip
minor_axis_earth = minoraxis(ellipsoid(WGS84Latest)) |> x-> uconvert(Unitful.km,x) |> ustrip
eccentricity²_earth= eccentricity²(ellipsoid(WGS84Latest))

##################################################################
#  EARTH MODEL
##################################################################
b_wgs84 = minor_axis_earth/major_axis_earth

#########################################################
# Normalization
#########################################################

#########################################################
temperature
#initialize atmosphere
atmosphere = StructArray(Matrix{AtmosphereProfile2D{Float64}}(undef, length(θᵢ),length(hᵢ)-1))


for j in axes(atmosphere,2)
    for i in axes(atmosphere,1)
        atmosphere[i,j]=AtmosphereProfile2D(
        let
          ismissing(refractive[j,i]) ? NaN : refractive[j,i]
        end,
        let
          ismissing(temperature[j,i]) ? NaN : temperature[j,i]
        end,
        let
          ismissing(pressure[j,i]) ? NaN : pressure[j,i]
        end,
        deg2rad(θᵢ[i]),
        hᵢ[j]
        )
    end
end

phi_cloves=θᵢ
z_cloves  =hᵢ

phi_cloves_ave=phi_cloves
z_cloves_ave  =(z_cloves[1:end-1]+z_cloves[2:end])./2
figure=Figure(size=(800/2,800))
begin
ax1=Axis(figure[1,1][1,1],title="refraction_index",xlabel="ϕ°",ylabel="h km")
ax2=Axis(figure[2,1][1,1],title="temperature",xlabel="ϕ°",ylabel="h km")
ax3=Axis(figure[3,1][1,1],title="pressure",xlabel="ϕ°",ylabel="h km")

s1=surface!(ax1,phi_cloves,z_cloves_ave,atmosphere.refraction_index_ave.-1.0)
s2=surface!(ax2,phi_cloves,z_cloves_ave,atmosphere.temperature_ave.+273.15)
s3=surface!(ax3,phi_cloves,z_cloves_ave,atmosphere.pressure_ave.*1e-2)
Colorbar(figure[1,1][1,2],s1;label="(n-1)")
Colorbar(figure[2,1][1,2],s2,label="°K")
Colorbar(figure[3,1][1,2],s3,label="hPa")
end
figure
save("./refraction_index_temperature_pressure_bianca.png",figure)

# Using equation from part 3 of Casia DEL006
point_ray_x    = Array{Float64}(undef,size(rays))
point_ray_y    = similar(point_ray_x)
direction_ray_x= similar(point_ray_x)
direction_ray_y= similar(point_ray_y)


total_scans = prod(size(point_ray_x))
max_iterations = 140

##########################################################
#  INPUT ARRAYS
##########################################################
# Retrieval array
retrieval = StructArray(Matrix{ResultsRayTracing{Float64}}(undef, total_scans,max_iterations+1));
# Input array
inputray  = StructArray(Vector{InputRay{Float64}}(undef, total_scans));
for (i,ray) in enumerate(rays)
  inputray.point_x[i] = ray.x
  inputray.point_y[i] = ray.y
  inputray.direction_x[i] = ray.direction_x
  inputray.direction_y[i] = ray.direction_y
end

e²_wgs84 = eccentricity²(ellipsoid(WGS84Latest))

ellipse(θ,h)= (cos(θ),sin(θ)*b_wgs84).+(b_wgs84*cos(θ),sin(θ)).*h./(sqrt(1-e²_wgs84*cos(θ)^2))

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
end
# tangent quote
tangent_quote =Array{Float64}(undef, total_scans)
t_out= similar(tangent_quote)
θ_out= similar(tangent_quote)
s_out= similar(tangent_quote)

fast_ray_tracing!(t_out::A,θ_out::A,s_out::A,
  inputray::IR,outputray::OR,atmosphere::ATM,θᵢ,hᵢ,
  tangent_quote::A;kwargs...) where {IR<:AbstractVector{InputRay{T}},
  OR<:AbstractMatrix{ResultsRayTracing{T}},
  ATM<:AbstractMatrix{AtmosphereProfile2D{T}},
  A<:AbstractVector{T}} where T=
fast_ray_tracing!(t_out,θ_out,s_out,
   inputray.point_x,inputray.point_y,
   inputray.direction_x,inputray.direction_y,
  inputray.n,inputray.θmin,inputray.θmax,inputray.ascending,
  atmosphere.refraction_index_ave,θᵢ,hᵢ,
  outputray.i,outputray.j,outputray.n,outputray.θ,outputray.t,outputray.h,
  outputray.point_x,outputray.point_y,outputray.direction_x,outputray.direction_y,
  tangent_quote,kwargs...)

scan_testing = 1 #total_scans
num_scans = 10
begin
test_retrieval=deepcopy(retrieval[1:scan_testing,1:num_scans+1])
test_inputray=deepcopy(inputray[1:scan_testing])
test_atmosphere=deepcopy(atmosphere)
tangent_quote =Array{Float64}(undef, scan_testing)
t_out= similar(tangent_quote)
θ_out= similar(tangent_quote)
s_out= similar(tangent_quote)
end

θᵢ_rad=@. deg2rad(θᵢ)

lines!(ax,[Point(
  test_inputray.point_x[1]+t*test_inputray.direction_x[1],
  test_inputray.point_y[1]+t*test_inputray.direction_y[1]
) for t in (0,1)])

fast_ray_tracing!(t_out,θ_out,s_out, test_inputray,test_retrieval,test_atmosphere,θᵢ_rad,hᵢ,tangent_quote)

scatter!(ax,[Point(x,y) for (x,y) in zip(test_retrieval.point_x[:],test_retrieval.point_y[:])])

test_retrieval[1,2].n
test_retrieval[1,3].n
test_retrieval
test_retrieval[1,:].point_y
tangent_quote
scatter!(ax,test_retrieval[1,:].point_x,test_retrieval[1,:].point_y,color=:black)

inputray.n
tangent_quote
atmosphere.refraction_index,ave[54,1]


save("./ray_tracing.png",figure)
reft_out
tangent_quote

tangent_quote


atmosphere
typeof(atmosphere)
typeof(inputray.point_x)

incident_refractive_index = ones(Float64,nrays)
adirection_x
ascending = [false for _ in 1:nrays]
incident_refractive_index = ones(Float64,nrays)
retrieval_i = zeros(Int,nrays,niterations+1)
retrieval_j = similar(retrieval_i)
retrieval_θ = zeros(Float64,nrays,niterations+1)
retrieval_t = similar(retrieval_θ)
retrieval_h = similar(retrieval_θ)
retrieval_n = similar(retrieval_θ)
retrieval_point_x = similar(retrieval_θ)
retrieval_point_y = similar(retrieval_θ)
retrieval_direction_x = similar(retrieval_θ)
retrieval_direction_y = similar(retrieval_θ)
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
apoint_x,apoint_y,adirection_x,adirection_y,
incident_refractive_index,
aθmin,aθmax,ascending,
atm_n,atm_θ,atm_h,
retrieval_i,retrieval_j,
retrieval_n,retrieval_θ,
retrieval_t,retrieval_h,
retrieval_point_x,retrieval_point_y,
retrieval_direction_x,retrieval_direction_y,
tangent_quote)





fig = Figure()
ax  = Axis(fig[1,1][1,1])
lines!(ax,[(cosd(θ),b_mine*sind(θ)) for θ in range(0,stop=365,length=100)],color=:black)
scatter!(ax,point_x,point_y,markersize=51)
