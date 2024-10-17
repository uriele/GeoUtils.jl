using GeoUtils,GeoMakie,WGLMakie
using ScatteredInterpolation
using Unitful
str=["vash_"*string(i)*"_prof" for i in 1:10];
df=GeoUtils.get_data("data",str,11);

df_filter= df[(30°.<df.lat.<45°) .& (20km.<df.alt.<100km),:];
_lon_length=122;
_lat_length=length(unique(df_filter.lat));
_alt_length=length(unique(df_filter.alt));

_lon=df_filter.lon;
_lat=df_filter.lat;
_alt=df_filter.alt;
_val=df_filter.vash_7_prof;



colormap = to_colormap(:plasma)
colormap[1] = RGBAf(0,0,0,0)

_val_m=reshape(_val,_lat_length,_lon_length,_alt_length);
_lon_m=ustrip.(reshape(_lon,_lat_length,_lon_length,_alt_length));
_lat_m=ustrip.(reshape(_lat,_lat_length,_lon_length,_alt_length));
_alt_m=ustrip.(reshape(_alt,_lat_length,_lon_length,_alt_length));
(lat_min,lat_max)=ustrip.(extrema(_lat));
(lon_min,lon_max)=ustrip.(extrema(_lon));


fig=Figure(resolution = (800, 800))
ax=GeoAxis(fig[1,1])
xlims!(ax, lon_min, lon_max)
ylims!(ax, lat_min, lat_max)

i=50
lines!(ax,GeoMakie.coastlines())
surface!(ax, Matrix(_lon_m[:,:,i]), Matrix(_lat_m[:,:,i]), Matrix(_val_m[:,:,i]), colormap = colormap)
surface!(ax, Matrix(_lon_m[:,:,60]), Matrix(_lat_m[:,:,60]), Matrix(_val_m[:,:,60]), colormap = colormap)

_t=convert.(Cartesian,LatLonAlt.(_lat,_lon,_alt))
_x=map(x->x.x,_t);
_y=map(x->x.y,_t);
_z=map(x->x.z,_t);

points_xyz =ustrip.( @. uconvert(km, [_x _y _z]))
val=df_filter.vash_7_prof;

full=Float32.([points_xyz val]);
unixmin,xmax=uconvert.(km,(minimum(_x),maximum(_x)))
ymin,ymax=uconvert.(km,(minimum(_y),maximum(_y)))
zmin,zmax=uconvert.(km,(minimum(_z),maximum(_z)))

using ScatteredInterpolation

itp=interpolate(NearestNeighbor() ,(points_xyz),df_filter.vash_7_prof);
