module LinearAlgebraExt
using LinearAlgebra, CoordRefSystems
import LinearAlgebra: dot, cross, norm, normalize, normalize!
import CoordRefSystems: convert

function convert(::Type{Cartesian}, coords::GeocentricLatLon{Datum})
    return convert(LatLon, coords) |>
              x-> convert(Cartesian,x)
    end
end
