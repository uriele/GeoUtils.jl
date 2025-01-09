using Test
using Interpolations
function mse(x,y)
  return sum((x-y).^2)/length(x)
end
const ACCETTABLE_ERROR_LINEAR_SHARP=1e-3

function xave(x)
   _x=(x[1:end-1]+x[2:end])
   push!(_x,x[1]+x[end])
   return _x./2
end
function generate_random_knots(xmin,xmax,n,idx)
  x_points=Vector{Float64}(undef,n)
  # create random knots vector
  [x_points[i]=xmin+(xmax-xmin)*sin(pi*(i-1)/(n-1)) for i in 1:n]
  return (xave(x_points),x_points[idx])
end


function xaveexp_linear(x)
  x[1:end-1].*exp.((log.(x[2:end]).-log.(x[1:end-1]))./2.0)
end


function generate_random_knots_exp(xmin,xmax,points,idx)
  n=length(points)
  h=points[end]-points[1]
  x_points=Vector{Float64}(undef,n)
  # create random knots vector
  lambda=log(xmax/xmin)/h
  [x_points[i]=xmin*exp(lambda*point) for (i,point) in enumerate(points)]
  return (x_points,x_points[idx])
end

function xaveexp_logaritmic(x)
  (x[1:end-1]-x[2:end])./log.(x[1:end-1]./x[2:end])
end

@testset "Interpolations" begin
  n=120;
  np=div(n,2);
  idx=vcat(1,unique(sort(rand(2:n-1,np))),n);
  points=LinRange(0,10,n);
  knots=points[idx];
  (p0,p1,t0,t1)=(1013.25,1000.0,15,30);
  (t_ave,t_knots)=generate_random_knots(t0,t1,length(points),idx);
  struct UnimplementedPressure<:GeoUtils.AbstractPressureInterpolation end

  @testset "interpolation with theta" begin
    (p_ave,p_knots)=generate_random_knots(p0,p1,length(points),idx)

    p_int,t_int=GeoUtils._interpolate_pt_theta(knots,p_knots,t_knots,points,Periodic())

    mse_p=mse(p_int,p_ave)
    mse_t=mse(t_int,t_ave)

    @test mse_p<ACCETTABLE_ERROR_LINEAR_SHARP
    @test mse_t<ACCETTABLE_ERROR_LINEAR_SHARP
  end
  @testset "interpolation with h" begin
    (p_points,p_knots)=generate_random_knots_exp(p0,p1,points,idx)
    p_ave1=xaveexp_linear(p_points)
    p_ave2=xaveexp_logaritmic(p_points)
    (p_int1,t_int1)=GeoUtils._interpolate_pt_h(LinearPressure(),knots,p_knots,t_knots,points,Flat())
    (p_int2,t_int2)=GeoUtils._interpolate_pt_h(LogarithmicPressure(),knots,p_knots,t_knots,points,Flat())

    mse_plinear    =mse(p_int1,p_ave1)
    mse_plogaritmic=mse(p_int2,p_ave2)
    mse_t1=mse(t_int1,t_ave[1:end-1])
    mse_t2=mse(t_int2,t_ave[1:end-1])

    @test mse_plinear<GeoUtils.atol(Float32)
    @test mse_plogaritmic<GeoUtils.atol(Float64)
    @test mse_t1<ACCETTABLE_ERROR_LINEAR_SHARP
    @test mse_t2<ACCETTABLE_ERROR_LINEAR_SHARP
    @test mse_t1==mse_t2
    @test_throws ArgumentError GeoUtils._interpolate_pt_h(UnimplementedPressure(),knots,p_knots,t_knots,points,Flat())
  end
end

@testset "Conversion LL2D to ECEF2D" begin
  🌎 = ellipsoid(NormalizedEarth)
  majoraxis_earth= majoraxis(🌎)
  minoraxis_earth= minoraxis(🌎)
  ecc²= eccentricity²(🌎)
  w=1.0
  z=1.0

  # Compare the results with what we know works the CoordRefSystems package

  lat,alt=convert(LatLonAlt{NormalizedEarth},Cartesian{NormalizedEarth}(w,0.0,z)) |> x->(latitude(x),altitude(x)) |> x-> ustrip.(x)
  lat2,alt2=convert(LLA2D{NormalizedEarth},ECEF2D{NormalizedEarth}(w,z)) |> x->(x.θ,x.h) |> x-> ustrip.(x)

  ww,zz=convert(Cartesian{NormalizedEarth},LatLonAlt{NormalizedEarth}(lat,0.0,alt)) |> x->(x.x,x.z) |> x-> ustrip.(x)
  ww2,zz2=convert(ECEF2D{NormalizedEarth},LLA2D{NormalizedEarth}(alt,lat)) |> x->(x.w,x.z) |> x-> ustrip.(x)
  h1,ang1=GeoUtils._extract_h_from_position(promote(w,z,majoraxis_earth,minoraxis_earth,ecc²)...)
  @test lat==lat2
  @test alt==alt2
  @test ww≈ww2
  @test zz≈zz2
  @test alt≈h1
  @test lat≈ang1
end
