using BenchmarkTools
using GeoUtils
using Distributed
using Polyester
N=10^5; #number of rays

rcoord=LinRange(2,10,500) |> x->rand([x;-x],4,N)
rays=map(x-> Ray2D(Vec2(x[1],x[2]),Vec2(x[3],x[4])),eachcol(rcoord))
cdist=zeros(size(rays))
rdist=similar(cdist)
sdist=similar(cdist)
C=Vec2(0.0,0.0)
D=Vec2(1.0,1.0)

function static_parallel_segment!(sdist,rays,C,D)
  Threads.@threads :static for i in 1:N
    @inbounds sdist[i]=@inline distance_from_segment(rays[i],C,D)
  end
end
function static_parallel_unit_circle!(cdist,rays)
  Threads.@threads :static for i in 1:N
    @inbounds cdist[i]=@inline distance_from_unit_circle(rays[i])
  end
end
function static_parallel_circle_radius!(rdist,rays)
  Threads.@threads :static for i in 1:N
    @inbounds rdist[i]=@inline distance_from_radius_new(rays[i],pi/4)
  end
end
function batched_segment!(sdist,rays,C,D)
  @batch for i in 1:N
    @inbounds sdist[i]=@inline @fastmath distance_from_segment(rays[i],C,D)
  end
end
function batched_unit_circle!(cdist,rays)
  @batch for i in 1:N
    @inbounds cdist[i]=@inline @fastmath distance_from_unit_circle(rays[i])
  end
end
function batched_circle_radius!(rdist,rays)
  @batch for i in 1:N
    @inbounds rdist[i]=@inline @fastmath distance_from_radius_new(rays[i],pi/4)
  end
end


SUITE= BenchmarkGroup()
SUITE["distance"] = BenchmarkGroup(["unit_circle","circle_radius","segment"])
SUITE["distance_parallel"] = BenchmarkGroup(["unit_circle","circle_radius","segment"])
SUITE["distance_distributed"] = BenchmarkGroup(["unit_circle","circle_radius","segment"])
SUITE["distance_fused"] = BenchmarkGroup(["unit_circle","circle_radius","segment"])
SUITE["distance_batched"] = BenchmarkGroup(["unit_circle","circle_radius","segment"])
SUITE["distance"]["unit_circle"]=@benchmarkable map!(r->distance_from_unit_circle(r),$cdist,$rays)
SUITE["distance"]["circle_radius"]=@benchmarkable map!(r->distance_from_radius(r,pi/4),$rdist,$rays)
SUITE["distance"]["segment"]=@benchmarkable map!(r->distance_from_segment(r,Vec2(0.0,0.0),Vec2(1.0,1.0)),$sdist,$rays)
SUITE["distance_distributed"]["unit_circle"]=@benchmarkable pmap(r->distance_from_unit_circle(r),$rays)
SUITE["distance_distributed"]["circle_radius"]=@benchmarkable pmap(r->distance_from_radius(r,pi/4),$rays)
SUITE["distance_distributed"]["segment"]=@benchmarkable pmap(r->distance_from_segment(r,Vec2(0.0,0.0),Vec2(1.0,1.0)),$rays)
SUITE["distance_parallel"]["unit_circle"]=@benchmarkable @fastmath static_parallel_unit_circle!($cdist,$rays)
SUITE["distance_parallel"]["circle_radius"]=@benchmarkable @fastmath static_parallel_circle_radius!($rdist,$rays)
SUITE["distance_parallel"]["segment"]=@benchmarkable @fastmath static_parallel_segment!($sdist,$rays,C,D)

SUITE["distance_batched"]["unit_circle"]=@benchmarkable @fastmath batched_unit_circle!($cdist,$rays)
SUITE["distance_batched"]["circle_radius"]=@benchmarkable @fastmath batched_circle_radius!($rdist,$rays)
SUITE["distance_batched"]["segment"]=@benchmarkable @fastmath batched_segment!($sdist,$rays,C,D)

SUITE["distance_fused"]["unit_circle"]=@benchmarkable cdist .= @. distance_from_unit_circle($rays)
SUITE["distance_fused"]["circle_radius"]=@benchmarkable rdist .= @. distance_from_radius($rays,pi/4)
SUITE["distance_fused"]["segment"]=@benchmarkable sdist .= @. distance_from_segment($rays,Vec2(0.0,0.0),Vec2(1.0,1.0))
#tune!(SUITE)
#results=run(SUITE;verbose=true)
#results
