@testset "index style" begin
  @test IndexStyle(SemiCircularMatrix{Float64,Matrix{Float64}}) == IndexCartesian()
end

@testset "construction ($T)" for T=(Float32,Float64,Int)
  data=rand(T,10,5)
  m=SemiCircularMatrix(data)
  @test SemiCircularMatrix{T}(data) == m
  @test SemiCircularMatrix{T,Matrix{T}}(data) == m
end

@testset "type stability ($T)" for T=(Float32,Float64,Int)
  data=rand(T,10,5)
  m=SemiCircularMatrix(data)

  @test @inferred(m[1]) isa T
  @test @inferred(axes(m)) isa Tuple{Vararg{AbstractUnitRange}}
  @test @inferred(similar(m)) isa typeof(m)
  @test @inferred(similar(m,T,(10,5))) isa SemiCircularMatrix{T}
  @test @inferred(similar(typeof(m),axes(m))) isa typeof(m)
end

@testset "indexing" begin
  data=reshape(1:20,10,2)
  m=SemiCircularMatrix(data)
  @test m[1]==1
  @test m[10,2]==20
  @test m[1]==m[21]
  @test m[10,3]==m[10,1]
  @test m[:,-1:1][1,1]==11

  @test IndexStyle(m)==IndexStyle(typeof(m))==IndexCartesian()
  @test isa(m,AbstractMatrix)
  @test isa(m,AbstractArray)
  @test !isa(m,AbstractVector)

end
