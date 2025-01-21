



# when the element is called shows the information
@inline function Base.show(io::IO,::MIME"text/plain", o::SatOrbit{T}) where T
  name=CoordRefSystems.prettyname(o)


  isnormalized= o._normalized ? "Normalized" : ""
  hasunit= o._normalized ? "" : "km"
  angle_type= islimb(o) ? "limb angle" : "nadir angle"
  _fieldnames=["coords",
              angle_type,
              "altitude",
              "orbital_coordinate"]
  ioctx = IOContext(io, :compact => true)

  _fieldvalues=[
    "(w= $(repr(o.w,context=ioctx))$(hasunit), z= $(repr(o.z,context=ioctx))$(hasunit))",
    "$(repr(o.ang,context=ioctx))°",
    "$(repr(o.h,context=ioctx))$(hasunit)",
    "$(repr(o.orbital_coordinate,context=ioctx))°"]



  print(io,"$isnormalized $name{$T}:")
  ioprintfields(ioctx,_fieldvalues,_fieldnames)
end


# when the element is parto of an array just show the type and if it is normalized
@inline function Base.show(io::IO, o::SatOrbit{T}) where T
  name=CoordRefSystems.prettyname(o)
  isnormalized= o._normalized ? "Normalized" : ""
  print(io,"$isnormalized $name{$T}")
end



@inline function ioprintfields(io::IOContext, fvalues::V,fnames::V; compact=false) where V
  if compact
    vals = map(zip(fvalues,fnames)) do (value, field)

      "$field: $value"
    end
    join(io, vals, ", ")
  else
    len = length(fnames)
    for (i, (value,field)) in enumerate(zip(fvalues,fnames))
      div = i == len ? "\n└─ " : "\n├─ "
      print(io, "$div$field: $value")
    end
  end
end
