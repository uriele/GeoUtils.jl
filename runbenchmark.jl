using PkgBenchmark
using GeoUtils
using Dates
config=BenchmarkConfig(
  juliacmd=`julia -O3`,
  env=Dict("JULIA_NUM_THREADS"=>4)
  )
config_baseline=BenchmarkConfig(
  id="baseline",
  juliacmd=`julia -O3`,
  env=Dict("JULIA_NUM_THREADS"=>1)
  )
  #judgment=judge(GeoUtils,config,config_baseline)
  results=benchmarkpkg("GeoUtils",config)
  cd("benchmark_results") do
    export_markdown("$(today() |> string |> x-> replace(x,"-"=>""))benchmark.md",results)
  end
