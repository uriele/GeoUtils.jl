
"""
  files_are_same(file1::String, file2::String)::Bool

Check if two files are the same by comparing their contents and sizes.
"""
function files_are_same(file1::String, file2::String)::Bool
  # Check if both files exist
  if !isfile(file1) || !isfile(file2)
      return false
  end

  # Compare file sizes first for a quick check
  if stat(file1).size != stat(file2).size
      return false
  end

  # Read and compare file contents
  open(file1, "r") do f1
      open(file2, "r") do f2
          return read(f1) == read(f2)
      end
  end
end

"""
  all_files_are_same(directory::String, file::String="in_lat.dat")::Bool

Check if all files with the same name in subdirectories are the same by comparing their contents and sizes.
"""
function all_files_are_same(directory::String,file::String="in_lat.dat")::Bool
  cd(directory) do
    flag = true
    x=readdir()

    for i in 2:length(x)
      file1=x[i]*"/"*file
      file2=x[i-1]*"/"*file
      files_are_same(file1,file2) || return false
    end
    return true
  end
end



function convert_to_array(file::String,skip=15)
  open(file, "r") do f
    data=readlines(f)[skip+1:end] |>
    x-> join(x," ") |> x-> String.(split(x)) |> x-> parse.(Float64,x)
  end
end

abstract type RealNumber end
struct isRealNumber end
struct isNotRealNumber end


RealNumber(x::T) where T = isNotRealNumber()
RealNumber(::Type{<:Real}) = isRealNumber()
RealNumber(x::T) where T<:Real = RealNumber(typeof(x))
RealNumber(x::U) where U<:Unitful.Quantity = RealNumber(Unitful.ustrip(x))

isRealNumber(x) = isa(x,Number) && !isa(x,Complex)
function fix_latitudes(lat::AbstractArray{T},orbital_coordinates=true) where T<:Number
  orbital_coordinates == false && return lat;
  # if the north pole is zero, we need to fix the latitudes
  # so that the north pole is at the top of the raster
  # also if the latitude is smaller than 180 => 90-lat
  # if the latitude is larger than 180 => lat-270
  GeoUtils._fix_latitudes(RealNumber(first(lat)),lat)
end

function fix_latitudes(lat::T,orbital_coordinates=true) where T<:Number
  orbital_coordinates == false && return lat;
  # if the north pole is zero, we need to fix the latitudes
  # so that the north pole is at the top of the raster
  # also if the latitude is smaller than 180 => 90-lat
  # if the latitude is larger than 180 => lat-270
  GeoUtils._fix_latitudes(RealNumber(lat),lat)
end

@inline _fix_latitudes(::isNotRealNumber,lat)= throw(ArgumentError("The latitude should be a real number"))

@inline function _fix_latitudes(::isRealNumber,lat::AbstractArray{T}) where T<:Real
  return [GeoUtils._conversion_latitude(_lat) for _lat in @. mod(lat,360)]
end

@inline function _fix_latitudes(::isRealNumber,lat::T) where T<:Real
  return GeoUtils._conversion_latitude(lat )
end

@inline function _fix_latitudes(::isRealNumber,lat::AbstractArray{U}) where U<:Unitful.Quantity
  _unit=Unitful.unit(first(lat));
  return GeoUtils._fix_latitudes(isRealNumber(),ustrip(lat)).*_unit
end

@inline function _fix_latitudes(::isRealNumber,lat::T) where T<:Unitful.Quantity
  _unit=Unitful.unit(lat);
  return GeoUtils._fix_latitudes(isRealNumber(),ustrip(lat)).*_unit
end

@inline _conversion_latitude(x::T) where T<:Real= x <= 180 ? (90 - x) : (x -270)

latitude(x::CoordRefSystems.Geographic)= x.lat;
longitude(x::CoordRefSystems.Geographic)=x.lon;
altitude(x::GeocentricLatLonAlt)=x.alt;
altitude(x::LatLonAlt)=x.alt;


function parse_vrm_profile(directory::String,file::String)
  _files=[];
  cd(directory) do
    open(file, "r") do f
      text=readlines(f)
      regex_dividers=regex=r"^\s{2,4}[0-9]{1,3}";
      regex_title=r"[0-9]{1,3}\s([a-z,A-Z,0-9\+]*)\s?";
      dividers=findall([!isnothing(match.(regex,line)) for line in text])
      @debug dividers

      titles=[string(title.captures[1])  for title in match.(regex_title,text[dividers])]

      push!(dividers,length(text)+1)

      for (i,t) in enumerate(titles)
        push!(_files,"vrm_"*string(t)*".dat");
        if !isfile("vrm_"*string(t)*".dat")
          open("vrm_"*string(t)*".dat","w") do f
            for str in text[dividers[i]+1:dividers[i+1]-1]
              write(f,str)
            end
          end
        end
      end

    end
  end
  return _files

end
