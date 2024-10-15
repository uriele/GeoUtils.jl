module GeoMakieExt
  using Makie,GeoMakie
  using CoordRefSystems
  using Makie: convert_arguments
  using Unitful
  @inline _cartesian_x(geo::CoordRefSystems.Cartesian)=geo.x;
  @inline _cartesian_y(geo::CoordRefSystems.Cartesian)=geo.y;
  @inline _cartesian_z(geo::CoordRefSystems.Cartesian3D)=geo.z;

  @inline _cartesian_x(vgeo::AbstractArray{<:CoordRefSystems.Cartesian})=[_cartesian_x(geo) for geo in vgeo];
  @inline _cartesian_y(vgeo::AbstractArray{<:CoordRefSystems.Cartesian})=[_cartesian_y(geo) for geo in vgeo];
  @inline _cartesian_z(vgeo::AbstractArray{<:CoordRefSystems.Cartesian3D})=[_cartesian_z(geo) for geo in vgeo];
  function Makie.convert_arguments(::PointBased, geo::CoordRefSystems.Geographic)
    _geo=convert(Cartesian,geo);
    @info _geo
    @info _cartesian_x(_geo)
    @info _cartesian_y(_geo)
    @info _cartesian_z(_geo)
    return convert_arguments(PointBased(),_cartesian_x(_geo), _cartesian_y(_geo), _cartesian_z(_geo));
  end

  function Makie.convert_arguments(::PointBased, vgeo::Vector{T}) where T<:CoordRefSystems.Geographic
    _vgeo=@. convert(CoordRefSystems.Cartesian,vgeo);
    @info _vgeo
    @info _cartesian_x(_vgeo)
    @info _cartesian_y(_vgeo)
    @info _cartesian_z(_vgeo)
    return convert_arguments(PointBased(),_cartesian_x(_vgeo), _cartesian_y(_vgeo), _cartesian_z(_vgeo));
  end

end
