module GeoUtils
  using Unitful
  using CoordRefSystems
  using DataFrames
  using Rasters
  import DimensionalData: @dim, XDim, YDim, ZDim
  include("Utils.jl")
  # This might be an extension of Rasters.jl
  @dim Lat YDim "Latitude"
  @dim Alt ZDim "Altitude"
  @dim Lon XDim "Longitude"


  """
    get_data(path::String,info::String="temp",skip=13;sink=DataFrame)

  This function reads the data from the files in the given path and returns a DataFrame,a
  Matrix, or a Raster. There are several assumptions made in this function:
  - The files are in the format in_lat.dat, in_alt.dat, and in_<info>.dat
  - The longitude for each file is the name of the directory divided by 10
  - The in_lat.dat and in_alt.dat files are the same for all directories
  """
  function get_data(path::String,info::String="temp",skip=13;sink=DataFrame)
    cd(path) do
      subdirs=readdir()
      # Check if all files are the same
      @assert(GeoUtils.all_files_are_same(".","in_lat.dat"),"in_lat.dat files are not the same")
      @assert(GeoUtils.all_files_are_same(".","in_alt.dat"),"in_alt.dat files are not the same")

      # get lat and alt from first file since they are the same for all files
      data_lat=GeoUtils.convert_to_array(subdirs[1]*"/in_lat.dat")
      data_alt=GeoUtils.convert_to_array(subdirs[1]*"/in_alt.dat",16)
      T=eltype(data_alt)
      # We assume that the name of the directory is the value of the longitude,
      # and we divide by 10 to get the actual value. According to Bianca the files
      # content is not used, so we do not need to project it.
      data_lon= @. parse(T,subdirs)/10;
      #lengths of the arrays
      lat_length = length(data_lat);
      alt_length = length(data_alt);
      lon_length = length(data_lon);

      data= Array{T,3}(undef,lat_length,lon_length.alt_length);

      for (i,directory) in enumerate(subdirs)
        data[:,i,:]=GeoUtils.convert_to_array(directory*"/in_"*info*".dat",skip);
      end

      if (sink != Raster)
        spatial_dimensions= [ [_lat,_lon,_alt] for _lat in data_lat, _lon in data_lon, _alt in data_alt ] |> x-> collect(Iterators.flatten(x))
        data=hcat(spatial_dimensions,data[:]);
      end

      if sink==Matrix
        return data
      end

      if sink==DataFrame
        return DataFrame(data,[:lat,:lon,:alt,Symbol(info)])
      end

      if (sink == Raster)
        _lat = Lat(data_lat);
        _alt = data_alt;
        _lon = data_lon;
        data=Raster(data,(_lat,_lon,_alt))
      end

      error("Unknown sink type")
    end
  end




end
