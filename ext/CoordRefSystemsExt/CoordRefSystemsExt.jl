module CoordRefSystemsExt
    using CoordRefSystems.
    using Unitful

    export EarthCenteredInertial, ECI

    struct EarthCenteredInertial <: CoordinateSystem
      altitude::T where T<:Unitful.
    end
end
