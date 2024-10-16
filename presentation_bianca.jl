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
	str=["vash_"*string(i)*"_prof" for i in 1:10];
	str_skip=repeat([11],length(str))
	str=[str;"aerden_17_prof";"aerden_prof";"ice_prof";"iden_prof";"pres";"temp";"vmr_prof"]
	str_skip=[str_skip;11;11;11;11;13;13;5]
	@info length(str)
	@info length(str_skip)
	
	df=unique(GeoUtils.get_data("data",str,str_skip;initial_lat=9));
	
	sort!(df,[:alt,:lat])
	_lat,_lon,_alt=unique.((df.lat,df.lon,df.alt));
	(_lon_min,_lon_max)=extrema(df.lon)
	df
	
end

# ╔═╡ 8fb3b1b8-d0ed-45f9-883b-db15e62596bb
Pkg.add("GeoMakie")

# ╔═╡ c3475450-4381-45a0-8837-e6891eaa8db4
using GeoMakie

# ╔═╡ f678896f-0826-438c-97dd-400b99387111


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
1. File $(@bind _info PlutoUI.Select(str))
"""

# ╔═╡ 9a95bae7-ede0-4766-8097-b63632d639cd
begin
	df_filtered=df[(_lat_min.<=df.lat.<_lat_max) .& (_alt_min.<=df.alt.<=_alt_max),:]
	sort!(df_filtered,[:alt,:lat]);
	
end

# ╔═╡ ee572a46-c410-4981-b40d-1ef490accc93
let
	str=names(df_filtered[:,Not(:lon,:lat,:alt)])
end

# ╔═╡ 9f9ee0e2-285b-435c-8a14-88b4da2ecb7a
names(df_filtered[:,Not(:lon,:lat,:alt)])

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
	
	__alt=__alt[1,1,:]
	md"Altitude $(@bind _selected_height_ PlutoUI.Slider(__alt))"	
end

# ╔═╡ 96e93603-098b-4a80-84c6-32220c208b85
begin
	using GADM
	import GeometryOps as GO
	ita_df = GADM.get("ITA",("Sicily")) |> DataFrame
	ita_centroid = GO.centroid(ita_df.geom)
	for reg in ["Campania","Basilicata","Apulia","Calabria"]
		r=GADM.get("ITA",reg) |> DataFrame
		global ita_df=vcat(ita_df,r)
	end

	fig = Figure()		
	
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

	ga=make_fig(fig,string(_selected_height_))
	idx=findfirst(x->x==_selected_height_,__alt)
	@info idx
	colormap=to_colormap(:plasma)
	colormap[1] = RGBAf(0,0,0,0)
	#zlims!(ga,extrema[__val])
	__val_lim=extrema(__val)
	__val1=(__val.-__val_lim[1])./(__val_lim[2]-__val_lim[1])
	scatter!(ga,Matrix(__lon[:,:,idx])[:],Matrix(__lat[:,:,idx])[:];color=Matrix(__val1[:,:,idx])[:],
	colormap=colormap)
	Colorbar(fig[1,2],colormap=colormap,limits=extrema(__val))
	#=
	surface!(ga,Matrix(__lon[:,:,idx]),Matrix(__lat[:,:,idx]),Matrix(__val[:,:,idx]);
	colormap=colormap)
	=#
	

	fig1=Figure()
	coord_lon=Observable(Matrix(__lon[:,:,idx])[:]);
	coord_lat=Observable(Matrix(__lat[:,:,idx])[:]);
	field    =Observable(Matrix(__val1[:,:,idx])[:]);
	
	frames =1:(length(unique(df_filtered.alt))-1);
	
	frame_rate= 10;
	
	title_obs = Observable("Altitude: "*string(__alt[1]))
	
	
	gx = make_fig(fig1,title_obs);
	scatter!(gx,coord_lon,coord_lat;color=field,colormap=colormap)
	
		
	Colorbar(fig1[1,2],colormap=colormap,limits=extrema(__val))
	fig1
	
	record(fig1,"fig1_bianca.mp4",frames;framerate=frame_rate) do frame
		
		
    	title_obs[] = "Altitude: "*string(__alt[frame])
		coord_lon[]=Matrix(__lon[:,:,frame])[:];
		coord_lat[]=Matrix(__lat[:,:,frame])[:];
		vv=Matrix(__val1[:,:,frame])[:];
		vv[end]=1.e-16;
		field[]    =vv;
	end
	
end

# ╔═╡ 960fc084-75e9-4557-8694-2831611ec665
extrema(__val1[:,:,end])

# ╔═╡ 0116c652-f918-4757-8328-accb91294cc3
ga.title.val

# ╔═╡ 60b3f5e1-3b2e-4fed-8b33-35c27aa94626
let
cc=GADM.get("ITA")
end

# ╔═╡ ea6a43c5-dc4f-4bdb-92d5-0a14fdbe0a82
extrema(__lon[:,:,10])

# ╔═╡ 10eab0ea-4fbb-4253-85c7-76c6ec9753c7
let
	
	fig2=Figure()
	ax2=Axis3(fig2[1,1],title=_info,
	xlabel="X km",
	ylabel="Y km",
	zlabel="Z km")

	
	__lat=df_filtered.lat
	__lon=df_filtered.lon
	__alt=df_filtered.alt
	__val=df_filtered[!,_info]
	__lat_length=length(unique(__lat));
	__lon_length=122;
	__alt_length=length(unique(__alt))
	
	__lat,__lon,__alt,__val=reshape.([__lat,__lon,__alt,__val],__lat_length,__lon_length,__alt_length)


	__x,__y,__z= LatLonAlt.(__lat,__lon,__alt) |>
	x-> convert.(Cartesian,x) |>
	x-> (uconvert.(km,[x.x for x in x]),uconvert.(km,[x.y for x in x]),uconvert.(km,[x.z for x in x]))

	@info size(__x)
	@info size(__y)
	@info size(__z)
	@info size(__val1)
	
	meshscatter!(ax2,ustrip.(__x)[:],ustrip.(__y)[:],ustrip.(__z)[:];color=__val1[:],colormap=colormap,alpha=0.5)
	fig2

	xlims!(ax2,4500,5000)
	ylims!(ax2,1000,1700)
	zlims!(ax2,3800,4400)
	record(fig2, "rotate.mp4", 1:120) do frame
    	ax2.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120)
    
	end
end

# ╔═╡ 72ec1dbd-54ad-4237-8fc4-2686461ff0e0
let
	
	fig2=Figure()
	ax2=Axis3(fig2[1,1],title=_info,
	xlabel="Latitude °",
	ylabel="Longitude °",
	zlabel="Altitude km")

	
	__lat=df_filtered.lat
	__lon=df_filtered.lon
	__alt=df_filtered.alt
	__val=df_filtered[!,_info]
	__lat_length=length(unique(__lat));
	__lon_length=122;
	__alt_length=length(unique(__alt))
	
	__lat,__lon,__alt,__val=reshape.([__lat,__lon,__alt,__val],__lat_length,__lon_length,__alt_length)


	
	meshscatter!(ax2,ustrip.(__lat)[:],ustrip.(__lon)[:],ustrip.(__alt)[:];color=__val1[:],colormap=colormap,alpha=0.8)
	fig2

	record(fig2, "rotate_deg.mp4", 1:120) do frame
    
    	ax2.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120)
    
	end
end

# ╔═╡ Cell order:
# ╠═f1dacc60-8be3-11ef-2407-d5dedc774620
# ╠═f678896f-0826-438c-97dd-400b99387111
# ╠═ee572a46-c410-4981-b40d-1ef490accc93
# ╟─8fb3b1b8-d0ed-45f9-883b-db15e62596bb
# ╟─c3475450-4381-45a0-8837-e6891eaa8db4
# ╠═0da23ae1-b5b9-4a06-913a-eba588452bb1
# ╠═d4be1fab-5a6a-44c0-8348-d19cf0dabb8c
# ╠═9f9ee0e2-285b-435c-8a14-88b4da2ecb7a
# ╠═9a95bae7-ede0-4766-8097-b63632d639cd
# ╠═c79e1274-9de4-4102-9dd7-94e0422e9187
# ╠═96e93603-098b-4a80-84c6-32220c208b85
# ╠═960fc084-75e9-4557-8694-2831611ec665
# ╠═0116c652-f918-4757-8328-accb91294cc3
# ╠═60b3f5e1-3b2e-4fed-8b33-35c27aa94626
# ╟─ea6a43c5-dc4f-4bdb-92d5-0a14fdbe0a82
# ╠═10eab0ea-4fbb-4253-85c7-76c6ec9753c7
# ╠═72ec1dbd-54ad-4237-8fc4-2686461ff0e0
