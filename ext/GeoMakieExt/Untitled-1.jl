using DataFrames


cd("./data") do
  readdir() |>
  x -> for xdir in x
    cd(xdir) do
      readdir() |>
      x -> for xfile in x
        println(xfile)
      end
    end
  end
end


using Base.Filesystem: stat

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

all_files_are_same("data","in_lat.dat")
all_files_are_same("data","in_alt.dat")
all_files_are_same("data","in_lon.dat")


function convert_to_array(file::String,skip=15)
  open(file, "r") do f
    data=readlines(f)[skip+1:end] |>
    x-> join(x," ") |> x-> String.(split(x)) |> x-> parse.(Float64,x)
  end
end


function get_data(path::String,info::String="temp",skip=13;sink=DataFrame)
  cd(path) do
    x=readdir()
    lon_length = length(x)
    # get lat and alt from first file since they are the same for all files
    data_lat=convert_to_array(x[1]*"/in_lat.dat")
    data_alt=convert_to_array(x[1]*"/in_alt.dat",16)
    lat_length = length(data_lat)
    alt_length = length(data_alt)

    data_lat_alt = Array{Float64,2}(undef,lat_length*alt_length,2)
    data_lat_alt[:,1]=repeat(data_lat,alt_length)
    data_lat_alt[:,2]=repeat(data_alt,inner=lat_length)

    data=Array{Float64,2}(undef,lat_length*alt_length*lon_length,4)
    dataview=@view data[:,:];

    slice_length=lat_length*alt_length;
    for directory in x
      dataview[1:slice_length,1:2]=data_lat_alt;
      dataview[1:slice_length,3]=repeat(convert_to_array(directory*"/in_lon.dat",15),alt_length)
      dataview[1:slice_length,4]=convert_to_array(directory*"/in_"*info*".dat",skip)
      dataview= @view dataview[slice_length+1:end,:]
    end
    if (sink == DataFrame)
      data=DataFrame(data,[:lat,:alt,:lon,Symbol(info)])
    end

    return data
  end
end

data=get_data("data")

using Plots

scatter(data[:,3])
