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






using Plots

scatter(data[:,3])
