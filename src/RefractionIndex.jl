abstract type AirModel end
struct Ciddor <: AirModel end
struct Borzsonyi <: AirModel end
abstract type Mathar<: AirModel end
struct Mathar1 <: Mathar end
struct Mathar2 <: Mathar end
struct Mathar3 <: Mathar end
struct Mathar4 <: Mathar end
#struct Mathar5 <: Mathar end
struct Birch <: AirModel end
struct Peck <: AirModel end

"""
  refractive_index([::AirModel=Ciddor];
  temperature::Union{T,Quantity{T,𝚯}}=15°C,
  pressure::Uniont{T,Quantity{T, 𝐌  𝐋 ^-1 𝐓 ^-2}=101325Pa,
  wavelength::Union{T,Quantity{T,𝐋}=0.532μm,
  humidity=0.0,
  CO2ppm=450)::T where T<:IEEEFloat


Calculate the refractive index of air at a given temperature, pressure, wavelength and humidity using an Air Model.
The units for each parameter, if not explicitly provided, are:
- temperature: °C
- pressure: Pa
- wavelength: μm
Humidity is a dimensionless quantity 0≤humidity≤1. The default value is 0.0. CO2 is the concentration of CO2 in ppm. The default value is 450ppm.


Model Material Refraction Index using the RefractiveIndex.info database[^1]
Current availeble models for air are:
- Ciddor 1996 (0.23-1.69μm) [^2]
- Borzsonyl et al. 2008 (0.4-1.0μm) [^3]
- Mathar 2007 (1.3-2.5μm) [^4]
- Mathar 2007 (2.8-4.2μm) [^4]
- Mathar 2007 (4.35-5.2μm) [^4]
- Mathar 2007 (7.5-14.1μm) [^4]
- Birch and Downs 1994 (0.350–0.650 µm) [^5]
- Peck and Reeder 1972 (0.185–1.7 µm) [^6]

[^1]: https://refrsactiveindex.info/
[^2]: https://refractiveindex.info/?shelf=main&book=Air&page=Ciddor
[^3]: https://refractiveindex.info/?shelf=main&book=Air&page=Borzsonyi
[^4]: https://refractiveindex.info/?shelf=main&book=Air&page=Mathar
[^5]: https://refractiveindex.info/?shelf=main&book=Air&page=Birch
[^6]: https://refractiveindex.info/?shelf=main&book=Air&page=Peck
"""
function refractive_index(model::AirModel=Ciddor();
  temperature::Union{T1,Unitful.Temperature{T1}}=15.0°C,
  pressure::Union{T2,Unitful.Pressure{T2}}=101325.0Pa,
  wavelength::Union{T3,Unitful.Length{T3}}=0.532μm,
  humidity::T4=0.0,
  CO2ppm::U=450) where {T1<:IEEEFloat,T2<:IEEEFloat,T3<:IEEEFloat,T4,U}

  # Convert all inputs to the correct units
  temperature =_convertto(temperature,°C)
  pressure =_convertto(pressure,Pa)
  wavelength =_convertto(wavelength,μm)
  (temperature,pressure,wavelength)=promote(temperature,pressure,wavelength)
  # Find the correct units for the output
  T=typeof(temperature)
  humidity=convert(T,humidity)
  CO2ppm=convert(T,CO2ppm)
  _refractive_index(model,temperature,pressure,wavelength,humidity,CO2ppm)
end

@inline _convertto(x,_unit)=dimension(x)==Unitful.NoDims ? x : ustrip.(uconvert(_unit,x))


@inline _refractive_index(model::AirModel,T,P,λ,h,CO2)= throw(ArgumentError("$(model) not implemented yet"))

@inline function _const_mathar(::Mathar,::Type{T}) where T<:IEEEFloat
  throw(ArgumentError("$(model) not implemented yet"))
end

@inline function _const_mathar(::Type{Mathar1},::Type{T}) where T<:IEEEFloat
  # valide 1.3-2.5μm
  cref = T.([ 0.200192e-3,  0.113474e-9,  -0.424595e-14,  0.100957e-16, -0.293315e-20,  0.307228e-24]) # cm^j
  cT   = T.([ 0.588625e-1, -0.385766e-7,   0.888019e-10, -0.567650e-13,  0.166615e-16, -0.174845e-20]) # cm^j · K
  cTT  = T.([-3.01513,      0.406167e-3,  -0.514544e-6,   0.343161e-9,  -0.101189e-12,  0.106749e-16]) # cm^j · K^2
  cH   = T.([-0.103945e-7,  0.136858e-11, -0.171039e-14,  0.112908e-17, -0.329925e-21,  0.344747e-25]) # cm^j · %^-1
  cHH  = T.([ 0.573256e-12, 0.186367e-16, -0.228150e-19,  0.150947e-22, -0.441214e-26,  0.461209e-30]) # cm^j · %^-2
  cp   = T.([ 0.267085e-8,  0.135941e-14,  0.135295e-18,  0.818218e-23, -0.222957e-26,  0.249964e-30]) # cm^j · Pa^-1
  cpp  = T.([ 0.609186e-17, 0.519024e-23, -0.419477e-27,  0.434120e-30, -0.122445e-33,  0.134816e-37]) # cm^j · Pa^-2
  cTH  = T.([ 0.497859e-4, -0.661752e-8,   0.832034e-11, -0.551793e-14,  0.161899e-17, -0.169901e-21]) # cm^j · K · %^-1
  cTp  = T.([ 0.779176e-6,  0.396499e-12,  0.395114e-16,  0.233587e-20, -0.636441e-24,  0.716868e-28]) # cm^j · K · Pa^-1
  cHp  = T.([-0.206567e-15, 0.106141e-20, -0.149982e-23,  0.984046e-27, -0.288266e-30,  0.299105e-34]) # cm^j · %^-1 · Pa^-1
  σ    = T.(1e4/2.25)    # cm^−1
  return hcat(cref,cT,cTT,cH,cHH,cp,cpp,cTH,cTp,cHp),σ
end

@inline function _const_mathar(::Type{Mathar2},::Type{T}) where T<:IEEEFloat
  # valide 2.8-4.2μm
  cref = T.([ 0.200049e-3,  0.145221e-9,   0.250951e-12, -0.745834e-15, -0.161432e-17,  0.352780e-20]) # cm^j
  cT   = T.([ 0.588432e-1, -0.825182e-7,   0.137982e-9,   0.352420e-13, -0.730651e-15, -0.167911e-18]) # cm^j · K
  cTT  = T.([-3.13579,      0.694124e-3,  -0.500604e-6,  -0.116668e-8,   0.209644e-11,  0.591037e-14]) # cm^j · K^2
  cH   = T.([-0.108142e-7,  0.230102e-11, -0.154652e-14, -0.323014e-17,  0.630616e-20,  0.173880e-22]) # cm^j · %^-1
  cHH  = T.([ 0.586812e-12, 0.312198e-16, -0.197792e-19, -0.461945e-22,  0.788398e-25,  0.245580e-27]) # cm^j · %^-2
  cp   = T.([ 0.266900e-8,  0.168162e-14,  0.353075e-17, -0.963455e-20, -0.223079e-22,  0.453166e-25]) # cm^j · Pa^-1
  cpp  = T.([ 0.608860e-17, 0.461560e-22,  0.184282e-24, -0.524471e-27, -0.121299e-29,  0.246512e-32]) # cm^j · Pa^-2
  cTH  = T.([ 0.517962e-4, -0.112149e-7,   0.776507e-11,  0.172569e-13, -0.320582e-16, -0.899435e-19]) # cm^j · K · %^-1
  cTp  = T.([ 0.778638e-6,  0.446396e-12,  0.784600e-15, -0.195151e-17, -0.542083e-20,  0.103530e-22]) # cm^j · K · Pa^-1
  cHp  = T.([-0.217243e-15, 0.104747e-20, -0.523689e-23,  0.817386e-26,  0.309913e-28, -0.363491e-31]) # cm^j · %^-1 · Pa^-1
  σ    = T.(1e4/3.4)    # cm^−1
  return hcat(cref,cT,cTT,cH,cHH,cp,cpp,cTH,cTp,cHp),σ
end

@inline function _const_mathar(::Type{Mathar3},::Type{T}) where T<:IEEEFloat
  # valide 4.35-5.2μm
  cref = T.([ 0.200020e-3,  0.275346e-9,   0.325702e-12, -0.693603e-14,  0.285610e-17,  0.338758e-18]) # cm^j
  cT   = T.([ 0.590035e-1, -0.375764e-6,   0.134585e-9,   0.124316e-11,  0.508510e-13, -0.189245e-15]) # cm^j · K
  cTT  = T.([-4.09830,      0.250037e-2,   0.275187e-6,  -0.653398e-8,  -0.310589e-9,   0.127747e-11]) # cm^j · K^2
  cH   = T.([-0.140463e-7,  0.839350e-11, -0.190929e-14, -0.121399e-16, -0.898863e-18,  0.364662e-20]) # cm^j · %^-1
  cHH  = T.([ 0.543605e-12, 0.112802e-15, -0.229979e-19, -0.191450e-21, -0.120352e-22,  0.500955e-25]) # cm^j · %^-2
  cp   = T.([ 0.266898e-8,  0.273629e-14,  0.463466e-17, -0.916894e-23,  0.136685e-21,  0.413687e-23]) # cm^j · Pa^-1
  cpp  = T.([ 0.610706e-17, 0.116620e-21,  0.244736e-24, -0.497682e-26,  0.742024e-29,  0.224625e-30]) # cm^j · Pa^-2
  cTH  = T.([ 0.674488e-4, -0.406775e-7,   0.289063e-11,  0.819898e-13,  0.468386e-14, -0.191182e-16]) # cm^j · K · %^-1
  cTp  = T.([ 0.778627e-6,  0.593296e-12,  0.145042e-14,  0.489815e-17,  0.327941e-19,  0.128020e-21]) # cm^j · K · Pa^-1
  cHp  = T.([-0.211676e-15, 0.487921e-20, -0.682545e-23,  0.942802e-25, -0.946422e-27, -0.153682e-29]) # cm^j · %^-1 · Pa^-1
  σ    = T.(1e4/4.8)    # cm^−1
  return hcat(cref,cT,cTT,cH,cHH,cp,cpp,cTH,cTp,cHp),σ
end

@inline function _const_mathar(::Type{Mathar4},::Type{T}) where T<:IEEEFloat
  # valide 7.5-14.1μm
  cref = T.([ 0.199885e-3,  0.344739e-9,  -0.273714e-12,  0.393383e-15, -0.569488e-17,  0.164556e-19]) # cm^j
  cT   = T.([ 0.593900e-1, -0.172226e-5,   0.237654e-8,  -0.381812e-11,  0.305050e-14, -0.157464e-16]) # cm^j · K
  cTT  = T.([-6.50355,      0.103830e-1,  -0.139464e-4,   0.220077e-7,  -0.272412e-10,  0.126364e-12]) # cm^j · K^2
  cH   = T.([-0.221938e-7,  0.347377e-10, -0.465991e-13,  0.735848e-16, -0.897119e-19,  0.380817e-21]) # cm^j · %^-1
  cHH  = T.([ 0.393524e-12, 0.464083e-15, -0.621764e-18,  0.981126e-21, -0.121384e-23,  0.515111e-26]) # cm^j · %^-2
  cp   = T.([ 0.266809e-8,  0.695247e-15,  0.159070e-17, -0.303451e-20, -0.661489e-22,  0.178226e-24]) # cm^j · Pa^-1
  cpp  = T.([ 0.610508e-17, 0.227694e-22,  0.786323e-25, -0.174448e-27, -0.359791e-29,  0.978307e-32]) # cm^j · Pa^-2
  cTH  = T.([ 0.106776e-3, -0.168516e-6,   0.226201e-9,  -0.356457e-12,  0.437980e-15, -0.194545e-17]) # cm^j · K · %^-1
  cTp  = T.([ 0.77368e-6,   0.216404e-12,  0.581805e-15, -0.189618e-17, -0.198869e-19,  0.589381e-22]) # cm^j · K · Pa^-1
  cHp  = T.([-0.206365e-15, 0.300234e-19, -0.426519e-22,  0.684306e-25, -0.467320e-29,  0.126117e-30]) # cm^j · %^-1 · Pa^-1
  σ    = T.(1e4/10.1)    # cm^−1
  return hcat(cref,cT,cTT,cH,cHH,cp,cpp,cTH,cTp,cHp),σ
end

#=
@inline function _const_mathar(::Type{Mathar5},::Type{T}) where T<:IEEEFloat
  # valide 16-24μm
  cref = T.([ 0.199436e-3,  0.299123e-8,  -0.214862e-10,  0.143338e-12,  0.122398e-14, -0.114628e-16]) # cm^j
  # something seems to be wrong with cT...
  cT   = T.([ 0.621723e-1, -0.177074e-4,   0.152213e-6,  -0.954584-9,   -0.996706e-11,  0.921476e-13]) # cm^j · K
  cTT  = T.([-23.2409,      0.108557,     -0.102439e-2,   0.634072e-5,   0.762517e-7,  -0.675587e-9 ]) # cm^j · K^2
  cH   = T.([-0.772707e-7,  0.347237e-9,  -0.272675e-11,  0.170858e-13,  0.156889e-15, -0.150004e-17]) # cm^j · %^-1
  cHH  = T.([-0.326604e-12, 0.463606e-14, -0.364272e-16,  0.228756e-18,  0.209502e-20, -0.200547e-22]) # cm^j · %^-2
  cp   = T.([ 0.266827e-8,  0.120788e-14,  0.522646e-17,  0.783027e-19,  0.753235e-21, -0.228819e-24]) # cm^j · Pa^-1
  cpp  = T.([ 0.613675e-17, 0.585494e-22,  0.286055e-24,  0.425193e-26,  0.413455e-28, -0.812941e-32]) # cm^j · Pa^-2
  cTH  = T.([ 0.375974e-3, -0.171849e-5,   0.146704e-7,  -0.917231e-10, -0.955922e-12,  0.880502e-14]) # cm^j · K · %^-1
  cTp  = T.([ 0.778436e-6,  0.461840e-12,  0.306229e-14, -0.623183e-16, -0.161119e-18,  0.800756e-20]) # cm^j · K · Pa^-1
  cHp  = T.([-0.272614e-15, 0.304662e-18, -0.239590e-20,  0.149285e-22,  0.136086e-24, -0.130999e-26]) # cm^j · %^-1 · Pa^-1
  σ    = T.(1e4/[20.0,16,24])    # cm^−1
  return hcat(cref,cT,cTT,cH,cHH,cp,cpp,cTH,cTp,cHp),σ
end
=#

@inline _mathar_range(::M) where {M<:Mathar}= throw(ArgumentError("$(model) not implemented yet"))

@inline _mathar_range(::Type{Mathar1})= 1e4./(1.3,2.5)
@inline _mathar_range(::Type{Mathar2})= 1e4./(2.8,4.2)
@inline _mathar_range(::Type{Mathar3})= 1e4./(4.35,5.2)
@inline _mathar_range(::Type{Mathar4})= 1e4./(7.5,14.1)
#@inline _mathar_range(::Type{Mathar5})= 1e4./(16,24)

@inline function _refractive_index(::M,temperature::T,pressure::T,wavelength::T,humidity::T,CO2ppm::T)::T where {T<:IEEEFloat,M<:Mathar}
  # adjustable parameters
    temperature_kelvin = ustrip(uconvert(K,temperature*°C))   # Temperature: K

    σ = ustrip(uconvert(cm^-1,T(1/wavelength)*μm^-1))
    # model parameters

    MATHAR_WAVELENGTHS=T.(_mathar_range(M))
    MATHAR_PARAMETERS,MATHAR_WAVELENGTH_REF = _const_mathar(M,T)

    @inline _between(x,a,b)=max(min(x,a),b)

    σ = _between(σ,MATHAR_WAVELENGTHS...)

    # Does not change with the model
    MATHAR_TEMPERATURE_KELVIN_REF = T(273.15+17.5) # K
    MATHAR_PRESSURE_REF = T(75000)       # Pa
    MATHAR_HUMIDITY_REF = T(0.10)          #%
    difftemp=T(1/temperature_kelvin-1/MATHAR_TEMPERATURE_KELVIN_REF)
    diffhum=T((humidity-MATHAR_HUMIDITY_REF)*100)
    diffpress=T(pressure-MATHAR_PRESSURE_REF)

    MATHAR_VECTOR=[T(1),
      difftemp,difftemp^2,
      diffhum,diffhum^2,
      diffpress,diffpress^2,
      difftemp*diffhum,
      difftemp*diffpress,
      diffhum*diffpress]
    #########################################
    val=MATHAR_PARAMETERS*MATHAR_VECTOR
    diffwavel=T(σ-MATHAR_WAVELENGTH_REF) |> x-> map(i->x^i,0:5)

    T(1)+reduce(+,val.*diffwavel)
end

@inline function _refractive_index(::Ciddor,temperature::T,pressure::T,wavelength::T,humidity::T,CO2ppm::T)::T where T<:IEEEFloat
  # constants involved in the calculation of group refractive index of dry air
  k0= T(238.0185) # µm^2;
  k1= T(5792105)  # µm^2;
  k2= T(57.362)   # µm^2;
  k3= T(167917)   # µm^2.
  ONE= T(1.0)
  ZERO= T(0.0)
  TENMINUS8= T(1e-8)
  REFRACTIVE_INDEX_CO2_CONTRIBUTION = T(0.534e-6) # unitless
  REFRACTIVE_INDEX_STANDARD_CO2PPM = T(450) # ppm
  WATER_VAPOR_CORRECTING_FACTOR = T(1.022e-8) # unitless
  w0 =  T(295.235)    # μm^-2 in the paper but it should be unitless
  w1 =  T(2.6422)     # μm^2
  w2 = -T(0.032380)   # μm^4
  w3 =  T(0.0004028)  # μm^6
  MOLAR_MASS_STANDARD_AIR = T(28.9635e-3) # g/mol
  MOLAR_MASS_CO2_CONTRIBUTION = T(12.011e-6) # g/mol
  MOLAR_MASS_STANDARD_CO2PPM = T(400) # ppm
  R_IDEAL_GAS_CONSTANT = T(8.314510) # J/(mol K)

  PRESSURE_AIR_REF= T(101325) # Pa
  PRESSURE_WATER_REF= T(1333) # Pa
  TEMPERATURE_KELVIN_AIR_REF= T(288.15) # K
  TEMPERATURE_CELSIUS_AIR_REF= ustrip(uconvert(°C,TEMPERATURE_KELVIN_AIR_REF*K)) # °C
  TEMPERATURE_WATER_KELVIN_REF= T(293.15) # K
  TEMPERATURE_WATER__CELSIUS_REF= ustrip(uconvert(°C,TEMPERATURE_WATER_KELVIN_REF*K)) # °C

  σ² = ONE/wavelength^2 # vacuum wavenumber µm^-2
  temperature_kelvin = ustrip(uconvert(K,temperature*°C)) # K
  # Equation 12 in the paper
  @inline function _compressibility(temperature_celsius::T,temperature_kelvin::T, pressure::T,humidity::T)::T where T<:IEEEFloat
    a0= T(1.58123e-6) # K/Pa
    a1=T(-2.9331e-8) # 1/Pa
    a2=T(1.1043e-10) # 1/K/Pa
    b0=T(5.707e-6) # K/Pa
    b1=T(-2.051e-8) # 1/Pa
    c0=T(1.9898e-4) # K/Pa
    c1=T(-2.376e-6) # 1/Pa
    d=T(1.83e-11) # K^2/Pa^-2
    e=T(-0.765e-8) # K^2/Pa^-2
    fraction = pressure/temperature_kelvin
    ONE-fraction*(
      a0+a1*temperature_celsius+a2*temperature_celsius^2+
      (b0+b1*temperature_celsius)*humidity+
      (c0+c1*temperature_kelvin)*humidity^2
    )+fraction^2*(d+e*humidity^2)
  end

  @inline function _density(temperature_kelvin::T,pressure::T,molar_mass::T,compressibility::T,CO2ppm::T)::T where T<:IEEEFloat
    pressure* molar_mass/(R_IDEAL_GAS_CONSTANT*temperature_kelvin*compressibility)*(ONE-CO2ppm)
  end

  # refractive index of standard air at 15°C and 101325 Pa 0% humidity and 450 ppm CO2
  n_standard_air_minus_one = (k1/(k0-σ²) + k3/(k2-σ²))*TENMINUS8  #standard dry air (n_as-1)
  n_standard_air_co2ppm_minus_one =n_standard_air_minus_one*(ONE+REFRACTIVE_INDEX_CO2_CONTRIBUTION*(CO2ppm-REFRACTIVE_INDEX_STANDARD_CO2PPM)) # standard dry air

  n_water_vapour_minus_one = WATER_VAPOR_CORRECTING_FACTOR*(w0+w1*σ²+w2*σ²^2+w3*σ²^3)*TENMINUS8 # water vapour

  #M_dry_air= MOLAR_MASS_STANDARD_AIR+MOLAR_MASS_CO2_CONTRIBUTION*(CO2ppm-MOLAR_MASS_STANDARD_CO2PPM) # molar mass of dry air, g/mol
  Z_AIR_REF = _compressibility(TEMPERATURE_CELSIUS_AIR_REF,TEMPERATURE_KELVIN_AIR_REF,PRESSURE_AIR_REF,ZERO) # 0% humidity
  #M_water = T(18.01528e-3) # molar mass of water vapor, g/mol
  Z_WATER_REF = _compressibility(TEMPERATURE_WATER__CELSIUS_REF,TEMPERATURE_WATER_KELVIN_REF,PRESSURE_WATER_REF,ONE) # 100% humidity

  Z_EFFECTIVE = _compressibility(temperature,temperature_kelvin,pressure,humidity) # effective compressibility
  #ρ_air_ref = _density(TEMPERATURE_KELVIN_AIR_REF,PRESSURE_AIR_REF,M_dry_air,Z_dry_air,0) # density of dry air, g/m^3
  #ρ_water_ref = _density(TEMPERATURE_WATER_KELVIN_REF,PRESSURE_WATER_REF,M_water,Z_water,1) # density of water vapor, g/m^3

  # Equation 13 in the paper
  # In theory it should be the ratio between the effective air density and the density of dry air at standard conditions. Where the density is computed as
  #
  # ρ_i = pressure * molar_mass/(R * temperature_kelvin * compressibility_i) * (1-CO2ppm)
  #
  #However, since Ma and Mw does not change, R is a constant and the compressibility is a function of the humidity in the air and not of the composition for a given pressure and temperature,
  # the ratio can be simplified to

  # ρ_air/ρ_dry_air = pressure/pressure_ref * temperature_ref/temperature * compressibility_ref/compressibility * (1-humidity)
  # ρ_water/ρ_water_ref = pressure/pressure_ref * temperature_ref/temperature * compressibility_ref/compressibility * humidity

  # Note: no checks are added since the temperature in kelvin is always positive and the reference pressure is always positive. The compressibility is also ~1 so no division by zero is possible

  ratio_air_density= pressure/PRESSURE_AIR_REF * TEMPERATURE_KELVIN_AIR_REF/temperature_kelvin * Z_AIR_REF/Z_EFFECTIVE * (ONE-humidity)
  ratio_water_density= pressure/PRESSURE_WATER_REF * TEMPERATURE_WATER_KELVIN_REF/temperature_kelvin * Z_WATER_REF/Z_EFFECTIVE * humidity

  return ONE+ratio_air_density*n_standard_air_minus_one+ratio_water_density*n_water_vapour_minus_one

end
