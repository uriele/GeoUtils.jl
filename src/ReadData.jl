
"""
get_data(path::String,info::String="temp",skip=13;initial_lat=0,initial_alt=0,sink=DataFrame,orbital_coordinates=true)

This function reads the data from the files in the given path and returns a DataFrame,a
Matrix, or a Raster. There are several assumptions made in this function:
- The files are in the format in_lat.dat, in_alt.dat, and in_<info>.dat
- The longitude for each file is the name of the directory divided by 10
- The in_lat.dat and in_alt.dat files are the same for all directories

- orbital_coordinates: If true, the latitude is fixed so that the north pole is at the top of the raster
"""
function get_data(path::String,Info::Vector{String}=["temp"],skip=13;initial_lat::Integer=1,initial_alt::Integer=1,sink=DataFrame,orbital_coordinates=true)
  cd(path) do
    @assert (length(skip)==1 || length(skip)==length(Info)) "The skip parameter should be a single value or a vector of the same length as the info parameter"

    if length(skip)==1
      skip=fill(skip[1],length(Info))
    end
    Info=copy(Info)
    subdirs=sort(parse.(Int,readdir()))
    # Check if all files are the same
    @assert(GeoUtils.all_files_are_same(".","in_lat.dat"),"in_lat.dat files are not the same")
    @assert(GeoUtils.all_files_are_same(".","in_alt.dat"),"in_alt.dat files are not the same")

    # get lat and alt from first file since they are the same for all files
    data_lat=GeoUtils.fix_latitudes(GeoUtils.convert_to_array(string(subdirs[1])*"/in_lat.dat"),orbital_coordinates);
    #idx_lat=sortperm(data_lat);
    #
    #data_lat=data_lat[idx_lat];

    data_alt=GeoUtils.convert_to_array(string(subdirs[1])*"/in_alt.dat",16);
    T=eltype(data_alt)
    # We assume that the name of the directory is the value of the longitude,
    # and we divide by 10 to get the actual value. According to Bianca the files
    # content is not used, so we do not need to project it.
    #data_lon= @. subdirs/10.0;
    #lengths of the arrays
    lat_length = length(data_lat);
    alt_length = length(data_alt);

    lon_length = length(subdirs);


    data= Array{T,3}(undef,alt_length,lon_length,lat_length);

    all_data=Dict{String,Array{T,3}}();
    for info in Info
      all_data[info]=similar(data);
    end
    spatial_dimension_lon=similar(data);

    fold=similar(data);

    [fold[:,i,:].=i for i in 1:lon_length]

    @inbounds for (i,directory) in enumerate(subdirs)
      for (sk,info) in zip(skip,Info)

        if info !="vmr_prof"

          all_data[info][:,i,:]=GeoUtils.convert_to_array(string(directory)*"/in_"*info*".dat",sk);
        else
          newfiles=parse_vrm_profile(string(directory),"in_vmr_prof.dat");
          for newfile in newfiles
            if !haskey(all_data,newfile[1:end-4])
              push!(Info,newfile[1:end-4]);
              all_data[newfile[1:end-4]]=similar(data);
            end
            all_data[newfile[1:end-4]][:,i,:]=GeoUtils.convert_to_array(string(directory)*"/"*newfile,0);
          end
        end
      end
      tmp=GeoUtils.convert_to_array(string(directory)*"/in_lon.dat");
      [spatial_dimension_lon[j,i,:]=tmp for j in 1:alt_length]
    end

    spatial_dimension_lat=similar(data);
    spatial_dimension_alt=similar(spatial_dimension_lat);




    [spatial_dimension_lat[j,i,:]=data_lat for i in 1:lon_length,j in 1:alt_length]
    [spatial_dimension_alt[:,j,i]=data_alt for i in 1:lat_length,j in 1:lon_length]

    spatial_dimension_lat2,spatial_dimension_lon2= GeocentricLatLonAlt.(spatial_dimension_lat[:],spatial_dimension_lon[:],spatial_dimension_alt[:]) |>
    x-> convert.(LatLonAlt,x) |>
    x-> (latitude.(x),longitude.(x));

    # FILTER INITIAL LAT AND ALT
    spatial_dimension_lat2= reshape(spatial_dimension_lat2,alt_length,lon_length,lat_length) |>
    x-> x[initial_alt:end,:,initial_lat:end] |> x-> x[:];
    spatial_dimension_lon2= reshape(spatial_dimension_lon2,alt_length,lon_length,lat_length) |>
    x-> x[initial_alt:end,:,initial_lat:end] |> x-> x[:];
    spatial_dimension_alt= reshape(spatial_dimension_alt,alt_length,lon_length,lat_length) |>
    x-> x[initial_alt:end,:,initial_lat:end] |> x-> x[:];
    fold= reshape(fold,alt_length,lon_length,lat_length) |>
    x-> x[initial_alt:end,:,initial_lat:end] |> x-> x[:];


    for info in Info
      all_data[info]= all_data[info] |>
      x-> x[initial_alt:end,:,initial_lat:end] ;
    end


    if sink==Matrix
      len=length(all_data);
      data=Array{T,2}(undef,lat_length*lon_length*alt_length,len+3);
      data[:,1]=ustrip.(spatial_dimension_lat2);
      data[:,2]=ustrip.(spatial_dimension_lon2);
      data[:,3]=spatial_dimension_alt;
      data[:,4]=fold[:];
      for (i,info) in enumerate(Info)
        data[:,i+3]=all_data[info][:];
      end
      return data
    end

    if sink==DataFrame


      return DataFrame(
        "lat"=>spatial_dimension_lat2,
        "lon"=>spatial_dimension_lon2,
        "alt"=>spatial_dimension_alt[:].*km,
        "lon_folder"=>fold[:],
        [info=>all_data[info][:] for info in Info]...
      )
    end

    error("Unknown sink type")
  end
end
get_data(path::String,info::String,skip::Integer=13;kwargs...)=get_data(path,[info],skip;kwargs...)

"""
unique_altitudes(data::DataFrame)
unique_altitudes(data::Matrix)

This function returns the unique altitudes in the data.
"""
function unique_altitudes(data::DataFrame)
  return unique(@. altitude(data[!,:coord]))
end
function unique_altitudes(data::Matrix)
  return unique(@. altitude(data[:,1]))
end

"""
unique_longitudes(data::DataFrame)
unique_longitudes(data::Matrix)

This function returns the unique longitudes in the data.
"""
function unique_longitudes(data::DataFrame)
  return unique(@. longitude(data[!,:coord]))
end
function unique_longitudes(data::Matrix)
  return unique(@. longitude(data[:,1]))
end

"""
unique_latitudes(data::DataFrame)
unique_latitudes(data::Matrix)

This function returns the unique latitudes in the data.
"""
function unique_latitudes(data::DataFrame)
  return unique(@. latitude(data[!,:coord]))
end
function unique_latitudes(data::Matrix)
  return unique(@. latitude(data[:,1]))
end
