
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
