



# when the element is called shows the information
@inline function Base.show(io::IO,::MIME"text/plain", o::SatOrbit{T}) where T
  name=CoordRefSystems.prettyname(o)


  normalized= o._normalized ? "Normalized" : ""
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



  print(io,"$normalized $name{$T}:")
  ioprintfields(ioctx,_fieldvalues,_fieldnames)
end

# when the element is part of an array just show the type and if it is normalized
@inline function Base.show(io::IO, o::SatOrbit{T}) where T
  name=CoordRefSystems.prettyname(o)
  normalized= o._normalized ? "Normalized" : ""
  print(io,"$normalized $name{$T}")
end

@inline _write_array_show(ioctx::IOContext,arr,hasunits,f=identity)=join(["$(repr(f(x),context=ioctx))$(hasunits)" for x in arr],",")


# when the element is part of an array just show the type and if it is normalized'

@inline function Base.show(io::IO,::MIME"text/plain", r::AR) where AR<:AbstractRay{T} where T
  name=CoordRefSystems.prettyname(AR)
  normalized=""
  hasunits="km"


  normalized=ifelse(r._normalized, "Normalized" , normalized)
  hasunits  =ifelse(r._normalized, "" , hasunits)

  ioctx = IOContext(io, :compact => true)
  _fieldnames=["origin","direction"]

  _fieldvalues=[_write_array_show(ioctx,r.origin,hasunits),  _write_array_show(ioctx,r.direction,"")]
  print(io,"$normalized $name{$T}:")
  ioprintfields(ioctx,_fieldvalues,_fieldnames)
end


function Base.show(io::IO,r::AR) where AR<:AbstractRay{T} where T
  name=CoordRefSystems.prettyname(AR)
  normalized=""
  hasunits="km"


  normalized=ifelse(r._normalized, "Normalized" , normalized)
  hasunits  =ifelse(r._normalized, "" , hasunits)

  ioctx = IOContext(io, :compact => true)
  _fieldnames=["origin","direction"]

  f=x-> round(x,digits=3)
  _fieldvalues=[_write_array_show(ioctx,r.origin,hasunits,f),  _write_array_show(ioctx,r.direction,"",f)]
  print(io,"$normalized $name{$T} (")
  ioprintfields(io,  _fieldvalues,_fieldnames, compact=true)
  print(io, ")")
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
