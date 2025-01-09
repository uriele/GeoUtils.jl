
#### Trait to define bend type
abstract type IntersectionStyle end
struct IsIntersection <: IntersectionStyle end
struct NoIntersection <: IntersectionStyle end

IntersectionStyle(x)= IntersectionStyle(typeof(x))
IntersectionStyle(::Type) = NoIntersection()

abstract type IntersectionLocation end
abstract type LevelIntersection <: IntersectionLocation end
abstract type RadiusIntersection  <: IntersectionLocation end
abstract type RadiusLevelIntersection <: IntersectionLocation end
abstract type LevelRadiusIntersection <: IntersectionLocation end

struct TopIntersection <: LevelIntersection end
struct BottomIntersection <: LevelIntersection end
struct LeftIntersection <: RadiusIntersection end
struct RightIntersection <: RadiusIntersection end
struct LeftTopIntersection <: RadiusLevelIntersection end
struct LeftBottomIntersection <: RadiusLevelIntersection end
struct RightTopIntersection <: RadiusLevelIntersection end
struct RightBottomIntersection <: RadiusLevelIntersection end
struct TopLeftIntersection <: LevelRadiusIntersection end
struct TopRightIntersection <: LevelRadiusIntersection end
struct BottomLeftIntersection <: LevelRadiusIntersection end
struct BottomRightIntersection <: LevelRadiusIntersection end

IntersectionStyle(::Type{<:LevelIntersection}) = IsIntersection()
IntersectionStyle(::Type{<:RadiusIntersection}) = IsIntersection()
IntersectionStyle(::Type{<:RadiusLevelIntersection}) = IsIntersection()
IntersectionStyle(::Type{<:LevelRadiusIntersection}) = IsIntersection()
