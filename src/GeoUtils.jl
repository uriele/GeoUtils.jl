module GeoUtils
  using Unitful
  using Unitful: ustrip,unit,°,m,km
  using CoordRefSystems
  using DataFrames
  include("Utils.jl")

  export Lat, Alt, Lon, get_data
  export RealNumber,isRealNumber,isNotRealNumber
  export latitude,longitude,altitude
  """
    get_data(path::String,info::String="temp",skip=13;sink=DataFrame,orbital_coordinates=true)

  This function reads the data from the files in the given path and returns a DataFrame,a
  Matrix, or a Raster. There are several assumptions made in this function:
  - The files are in the format in_lat.dat, in_alt.dat, and in_<info>.dat
  - The longitude for each file is the name of the directory divided by 10
  - The in_lat.dat and in_alt.dat files are the same for all directories

  - orbital_coordinates: If true, the latitude is fixed so that the north pole is at the top of the raster
  """
  function get_data(path::String,info::String="temp",skip=13;sink=DataFrame,orbital_coordinates=true)
    cd(path) do
      subdirs=sort(parse.(Int,readdir()))
      # Check if all files are the same
      @assert(GeoUtils.all_files_are_same(".","in_lat.dat"),"in_lat.dat files are not the same")
      @assert(GeoUtils.all_files_are_same(".","in_alt.dat"),"in_alt.dat files are not the same")

      # get lat and alt from first file since they are the same for all files
      data_lat=GeoUtils.fix_latitudes(GeoUtils.convert_to_array(string(subdirs[1])*"/in_lat.dat"),orbital_coordinates);
      idx_lat=sortperm(data_lat);
      data_lat=data_lat[idx_lat];

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

      data= Array{T,3}(undef,lat_length,lon_length,alt_length);

      spatial_dimension_lon=similar(data);

      @inbounds for (i,directory) in enumerate(subdirs)
        data[:,i,:]=GeoUtils.convert_to_array(string(directory)*"/in_"*info*".dat",skip);
        tmp_lon=@view spatial_dimension_lon[:,i,:];

        tmp_lon=repeat(GeoUtils.convert_to_array(string(directory)*"/in_lon.dat"),alt_length);
      end

      data=data[idx_lat,:,:];
      #if (sink != Raster)
        spatial_dimension_lat=similar(data);
        spatial_dimension_alt=similar(spatial_dimension_lat);
        [spatial_dimension_lat[:,i,j]=data_lat for i in 1:lon_length,j in 1:alt_length]
        [spatial_dimension_alt[i,j,:]=data_alt for i in 1:lat_length,j in 1:lon_length]
        @info "Data dimensions: $(size(data[:]))"
        data=hcat(GeocentricLatLonAlt.(spatial_dimension_lat[:],spatial_dimension_lon[:],spatial_dimension_alt[:]),data[:]);

        @info size(data)
      #end

      if sink==Matrix
        return data
      end

      if sink==DataFrame
        return DataFrame(data,[:coord,Symbol(info)])
      end

      #=
      if (sink == Raster)
        _lat = Lat(data_lat)#.*°);
        _alt = Alt(data_alt)#.*km);
        _lon = Lon(data_lon)#.*°);
        return Raster(data,(_lat,_lon,_alt))
      end
      =#
      error("Unknown sink type")
    end
  end

  """

"""


end
