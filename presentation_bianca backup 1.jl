### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ f1dacc60-8be3-11ef-2407-d5dedc774620
begin
    import Pkg
    # careful: this is _not_ a reproducible environment
    # activate the global environment
    Pkg.activate(".")
	Pkg.add("GeoMakie")
	Pkg.add("PlutoUI")
	Pkg.add("GADM")
	Pkg.add("GeometryOps")
	Pkg.instantiate()
	using GeoUtils
	using Makie
	using GLMakie
	using ScatteredInterpolation
	using Unitful
	using PlutoUI
	using ImageFiltering
	using Printf
	using GADM
	Pkg.add("GeoMakie")
	using GeoMakie
	import GeometryOps as GO
	str=["vash_"*string(i)*"_prof" for i in 1:10];
	str_skip=repeat([11],length(str))
	str=[str;"aerden_17_prof";"aerden_prof";"ice_prof";"iden_prof";"pres";"temp";"vmr_prof"]
	str_skip=[str_skip;11;11;11;11;13;13;5]
	@info length(str)
	@info length(str_skip)
	
	df=unique(GeoUtils.get_data("data",str,str_skip;initial_lat=9));
	
	sort!(df,[:alt,:lon_folder,:lat])
	_lat,_lon,_alt=unique.((df.lat,df.lon,df.alt));
	(_lon_min,_lon_max)=extrema(df.lon)
	df
	
end

# ╔═╡ 8fb3b1b8-d0ed-45f9-883b-db15e62596bb
Pkg.add("GeoMakie")

# ╔═╡ e1de2c8e-7aa2-4043-9ee1-629aab06414e
Pkg.add("ImageFiltering")

# ╔═╡ 8d049e9e-08fe-4422-9921-1a0f7759d18c
val_plot=names(df[:,Not(:lon,:lat,:alt,:lon_folder)])


# ╔═╡ 0da23ae1-b5b9-4a06-913a-eba588452bb1
md"""
1. Filter the data in the range of Latitude 35° to 45°
2. Filter the altitude between 20km to 100km (all other values are of the order of 10^-6)
"""

# ╔═╡ d4be1fab-5a6a-44c0-8348-d19cf0dabb8c
md"""
# Range of Interest
1. Latitude Minimum $(@bind _lat_min PlutoUI.Select(_lat,default=_lat[1]))
1. Latitude Maximum $(@bind _lat_max PlutoUI.Select(_lat,default=_lat[end-1]))
1. Altitude Minimum $(@bind _alt_min PlutoUI.Select(_alt,default=_alt[21]))
1. Altitude Maximum $(@bind _alt_max PlutoUI.Select(_alt,default=_alt[end]))
1. File $(@bind _info PlutoUI.Select(val_plot))
"""

# ╔═╡ 9a95bae7-ede0-4766-8097-b63632d639cd
begin
	df_filtered=df[(_lat_min.<=df.lat.<_lat_max) .& (_alt_min.<=df.alt.<=_alt_max),:]
	#sort!(df_filtered,[:alt,:lat]);
	df_filtered[:,[:lat,:lon,:alt,Symbol(_info)]]
end

# ╔═╡ c79e1274-9de4-4102-9dd7-94e0422e9187
begin
	__lat=ustrip.(df_filtered.lat)
	__lon=ustrip.(df_filtered.lon)
	__alt=df_filtered.alt
	__val=df_filtered[!,_info]
	__lat_length=length(unique(__lat));
	__lon_length=122;
	__alt_length=length(unique(__alt))
	
	__lat,__lon,__alt,__val=reshape.([__lat,__lon,__alt,__val],__lat_length,__lon_length,__alt_length)
	__alt1=copy(__alt)
	__lat1=copy(__lat)
	__alt=__alt[1,1,:]
	__lat1=unique(__lat1)	
	__val1=copy(__val);
	__val1.=(__val1.-minimum(__val1))./(maximum(__val1)-minimum(__val1));
end

# ╔═╡ 9fbd0b1a-5f5a-406b-b7ef-73633ecbbc1e
begin
	points=Matrix([hcat(__lat[:]) hcat(__lon[:])  ustrip.(hcat(__alt1[:]))]');
	samples=__val1[:];
end

# ╔═╡ 48eb9f7a-c062-4b38-a7bc-1ee2572d5283
begin
	itp=interpolate(NearestNeighbor(),points,samples)
end

# ╔═╡ ec07b2c4-4882-4318-b1e3-2986497e9e8f
begin
	range_lat=extrema(__lat)
	range_lon=extrema(__lon)
	range_lat=(35.0,45.0)
	range_lon=(13.0,21.0)
	range_alt=extrema(ustrip.(__alt1))
	
	n=200;
	m=600;
	Lat=range(range_lat[1],stop=range_lat[2],length=n)
	Lon=range(range_lon[1],stop=range_lon[2],length=n)
	Alt=range(range_alt[1],stop=range_alt[2],length=m)
	len_lat,len_lon,len_alt=(length(Lat),
		length(Lon),length(Alt));
	len_tot=len_lat*len_lon*len_alt
	gridPoints_lat=Array{Float32,3}(undef,len_lat,len_lon,len_alt)
	gridPoints_lon=similar(gridPoints_lat)
	gridPoints_alt=similar(gridPoints_lat)
	[(gridPoints_lat[:,j,k].=Lat,gridPoints_lon[i,:,k].=Lon,gridPoints_alt[i,j,:].=Alt) for i in 1:n, j in 1:n,k in 1:m]
	kernel=Kernel.gaussian((2,2,2))
	val_itp=evaluate(itp,[gridPoints_lat[:] gridPoints_lon[:] gridPoints_alt[:]]') |>
	x-> reshape(x,len_lat,len_lon,len_alt) 
	#|>
	#x-> imfilter(x,kernel)
	
end

# ╔═╡ 7d9f5e20-9c3e-4df0-8703-a5092df4f36b
begin
	ita_df = GADM.get("ITA",("Sicily")) |> DataFrame
	ita_centroid = GO.centroid(ita_df.geom)
	for reg in ["Campania","Basilicata","Apulia","Calabria"]
		r=GADM.get("ITA",reg) |> DataFrame
		global ita_df=vcat(ita_df,r)
	end

	function make_fig(_fig,_title)
		ga = GeoAxis(
		    _fig[1, 1];
			title=_title,
		    dest = "+proj=ortho +lon_0=$(ita_centroid[1]) +lat_0=$(ita_centroid[2])"
		)	
		poly!(
		    ga, ita_df.geom;
		    color = 1:size(ita_df, 1),
		    strokecolor = :blue, strokewidth = 1, shading = NoShading
		    )
		ylims!(ga,[36,41])
		xlims!(ga,[11,21])
		return ga	
	end

	colormap=to_colormap(:plasma)
	colormap[1] = RGBAf(0,0,0,0)
	_val_itp=val_itp./maximum(val_itp);
	
	fig1=Figure()
	coord_lon=Observable(Matrix(gridPoints_lon[:,:,1]));
	coord_lat=Observable(Matrix(gridPoints_lat[:,:,1]));
	field    =Observable(Matrix(_val_itp[:,:,1]));

	num_altitudes=length(gridPoints_alt[1,1,:])
	
	frames =1:num_altitudes;
	
	frame_rate= 60;
	
	title_obs = Observable(@sprintf("Altitude: %.2f ",gridPoints_alt[1,1,1]))
	
	
	gx = make_fig(fig1,title_obs);
	surface!(gx,coord_lon,coord_lat,field,colormap=colormap,colorrange=(0,1))
	
	xlims!(gx,11,21)
	ylims!(gx,35,41)
	Colorbar(fig1[1,2],colormap=colormap,limits=extrema(__val),label="μg/kg-dryair")
	
	record(fig1,"unfiltered_"*_info*"_overmap.mp4",frames;framerate=90) do frame
		idx=frame
    	title_obs[] = @sprintf("Altitude: %.2f ",gridPoints_alt[1,1,idx])
		coord_lon[]=Matrix(gridPoints_lon[:,:,idx]);
		coord_lat[]=Matrix(gridPoints_lat[:,:,idx]);
		vv=Matrix(_val_itp[:,:,idx]);
		vv[end]=1.e-6;
		field[]    =vv;
	end
	
end

# ╔═╡ 2027d64b-8140-4155-878b-f4ae0aebbcc7
let
	fig2=Figure(size=(1200,600),title=_info)
	
	ax1=Axis(fig2[1,1],
	xlabel="Latitude °",
	ylabel="Altitude km")
	
	ax2=Axis3(fig2[1,2],
	xlabel="Latitude °",
	ylabel="Longitude °",
	zlabel="Altitude km",
	azimuth=-1/6*pi,
	aspect=:data)
	


	
	#xlims!(ax1,range_lat)
	#ylims!(ax1,0,20)
	#ylims!(ax2,14,20.5)
	current_p=Point3(range_lat[1],Lon[1],range_alt[1])
	current_rect=Observable(Rect3f(current_p,(41-37,0.01,20)))
	
	current_range=Observable((Lon[1],range_lon[2]))
	current_val  =Observable(_val_itp[:,:,:])
	val2d=Observable(Matrix(_val_itp[:,1,:]))

	current_x=Observable(Matrix(gridPoints_lat[:,1,:]))
	current_y=Observable(Matrix(gridPoints_lon[:,1,:]))
	current_z=Observable(Matrix(gridPoints_alt[:,1,:]))
	


	
	
	volume!(ax2,range_lat,current_range,range_alt,current_val;shading=NoShading,colormap=colormap,colorrange=(0,1))
	
	mesh!(ax2,current_rect,color=:black,alpha=0.8)
	surface!(ax2,current_x,current_y,current_z;color=val2d,shading=NoShading,colormap=colormap,colorrange=(0,1))
	
	surface!(ax1,Matrix(gridPoints_lat[:,1,:]),Matrix(gridPoints_alt[:,1,:]);color=val2d,shading=NoShading,colormap=colormap,colorrange=(0,1))
	frames=1:length(Lon)
	
	
	record(fig2,"unfiltered_"*_info*"cross_nube_Lon.mp4",frames;framerate=20) do frame
			
		current_p=Point3(range_lat[1],Lon[frame],range_alt[1])
		current_rect[]=Rect3f(current_p,(range_lat[2]-range_lat[1],0.01,20))
		current_x[]=Matrix(gridPoints_lat[:,frame,:])
		current_y[]=Matrix(gridPoints_lon[:,frame,:])
		current_z[]=Matrix(gridPoints_alt[:,frame,:])
		current_range[]=(Lon[frame],range_lon[2])
		current_val[]=_val_itp[:,frame:end,:]
		val2d[]=Matrix(_val_itp[:,frame,:])
	end

	



end

# ╔═╡ 10eab0ea-4fbb-4253-85c7-76c6ec9753c7
let
	fig2r=Figure(size=(1200,1200),title=_info)
	
	ax2r=Axis3(fig2r[1,1],
	xlabel="Latitude °",
	ylabel="Longitude °",
	zlabel="Altitude km",
	aspect=:data)
	
	current_p=Point3(range_lat[1],Lon[1],range_alt[1])
	current_rect=Observable(Rect3f(current_p,(41-37,0.01,20)))
	
	current_range=(Lon[1],range_lon[2])
	current_val  =_val_itp[:,:,:]
	val2d        =Matrix(_val_itp[:,1,:])


	volume!(ax2r,range_lat,current_range,range_alt,current_val;shading=NoShading,colormap=colormap,colorrange=(0,1))
	
	record(fig2r, "unfiltered_"*_info*"_rotate_deg.mp4", 1:120) do frame
    	ax2r.azimuth[] = 1.7pi + 0.6 * sin(2pi * frame / 120)
	end

end

# ╔═╡ Cell order:
# ╠═f1dacc60-8be3-11ef-2407-d5dedc774620
# ╟─8fb3b1b8-d0ed-45f9-883b-db15e62596bb
# ╟─8d049e9e-08fe-4422-9921-1a0f7759d18c
# ╠═0da23ae1-b5b9-4a06-913a-eba588452bb1
# ╟─d4be1fab-5a6a-44c0-8348-d19cf0dabb8c
# ╠═9a95bae7-ede0-4766-8097-b63632d639cd
# ╟─ebee726a-f5ed-4bf9-9226-22923f9e28dd
# ╟─c79e1274-9de4-4102-9dd7-94e0422e9187
# ╠═9fbd0b1a-5f5a-406b-b7ef-73633ecbbc1e
# ╠═48eb9f7a-c062-4b38-a7bc-1ee2572d5283
# ╠═ec07b2c4-4882-4318-b1e3-2986497e9e8f
# ╠═7d9f5e20-9c3e-4df0-8703-a5092df4f36b
# ╠═2027d64b-8140-4155-878b-f4ae0aebbcc7
# ╠═10eab0ea-4fbb-4253-85c7-76c6ec9753c7
