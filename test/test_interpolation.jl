testfile="../data_atm/INP_FILES"

# Read atmosphere
if isdir(testfile)
  #read local atmosphere



  atmosphere=read_local_atmosphere(testfile);
  @testset "Test Ray Tracing File Reading" begin
    @testset "Test new atmosphere reader" begin
      atmosphere_old_reader=GeoUtils.read_local_atmosphere_old(testfile);
      @test atmosphere_old_reader == atmosphere
    end

    orbit_old=GeoUtils.read_orbit_old("$testfile/orbit.dat")
    @testset "Test new orbit" begin
      transposed_orbit=read_orbit("$testfile/orbit.dat")

      orbit=StructArray{SatOrbit{Float64},1}(undef,size(transposed_orbit,2),size(transposed_orbit,1))

      for (o1,to1) in zip(StructArrays.components(orbit),StructArrays.components(transposed_orbit))
        @inbounds o1.=transpose(to1)
      end
        @test orbit == orbit_old
    end
    orbit_normalized=deepcopy(orbit_old)
    normalize_datum!(orbit_normalized)

    for (datumtype,orb) in zip(["WGS84Latest","Normalized"],[orbit_old,orbit_normalized])
      @testset "Test $datumtype orbit from rays" begin
        rays=create_bundle_rays(orb);
        prealloc_rays=StructArray{Ray2D{Float64}}(undef,size(orb));
        rays_from_orbit!(prealloc_rays,orb);
        if datumtype=="Normalized"
          @test all(rays._normalized .== true)
          @test all(prealloc_rays._normalized .== true)
        else
          @test all(rays._normalized .== false)
          @test all(prealloc_rays._normalized .== false)
        end

        @test rays == prealloc_rays
      end
    end

  end

  #define interpolation points
  θᵢ=" 0.00      1.00      2.00      3.00      4.00      5.00      6.00
  7.00      8.00      9.00     10.00     11.00     12.00     13.00
  14.00     15.00     16.00     17.00     18.00     19.00     20.00
  21.00     22.00     23.00     24.00     25.00     26.00     27.00
  28.00     29.00     30.00     31.00     32.00     33.00     34.00
  35.00     36.00     37.00     38.00     39.00     40.00     41.00
  42.00     42.50     43.00     43.50     44.00     44.50     45.00
  45.25     45.50     45.75     46.00     46.10     46.20     46.30
  46.40     46.50     46.60     46.70     46.80     46.90     47.00
  47.10     47.20     47.30     47.40     47.50     47.60     47.70
  47.80     47.90     48.00     48.10     48.20     48.30     48.40
  48.50     48.60     48.70     48.80     48.90     49.00     49.10
  49.20     49.30     49.40     49.50     49.60     49.70     49.80
  49.90     50.00     50.10     50.20     50.30     50.40     50.50
  50.60     50.70     50.80     50.90     51.00     51.10     51.20
  51.30     51.40     51.50     51.60     51.70     51.80     51.90
  52.00     52.10     52.20     52.30     52.40     52.50     52.60
  52.70     52.80     52.90     53.00     53.10     53.20     53.30
  53.40     53.50     53.60     53.70     53.80     53.90     54.00
  54.25     54.50     54.75     55.00     55.25     55.50     56.00
  56.50     57.00     57.50     58.00     59.00     60.00     61.00
  62.00     63.00     64.00     65.00     66.00     67.00     68.00
  69.00     70.00     71.00     72.00     73.00     74.00     75.00
  76.00     77.00     78.00     79.00     80.00     81.00     82.00
  83.00     84.00     85.00     86.00     87.00     88.00     89.00
  90.00     91.00     92.00     93.00     94.00     95.00     96.00
  97.00     98.00     99.00    100.00    101.00    102.00    103.00
  104.00    105.00    106.00    107.00    108.00    109.00    110.00
  111.00    112.00    113.00    114.00    115.00    116.00    117.00
  118.00    119.00    120.00    121.00    122.00    123.00    124.00
  125.00    126.00    127.00    128.00    129.00    130.00    131.00
  132.00    133.00    134.00    135.00    136.00    137.00    138.00
  139.00    140.00    141.00    142.00    143.00    144.00    145.00
  146.00    147.00    148.00    149.00    150.00    151.00    152.00
  153.00    154.00    155.00    156.00    157.00    158.00    159.00
  160.00    161.00    162.00    163.00    164.00    165.00    166.00
  167.00    168.00    169.00    170.00    171.00    172.00    173.00
  174.00    175.00    176.00    177.00    178.00    179.00    180.00
  181.00    182.00    183.00    184.00    185.00    186.00    187.00
  188.00    189.00    190.00    191.00    192.00    193.00    194.00
  195.00    196.00    197.00    198.00    199.00    200.00    201.00
  202.00    203.00    204.00    205.00    206.00    207.00    208.00
  209.00    210.00    211.00    212.00    213.00    214.00    215.00
  216.00    217.00    218.00    219.00    220.00    221.00    222.00
  223.00    224.00    225.00    226.00    227.00    228.00    229.00
  230.00    231.00    232.00    233.00    234.00    235.00    236.00
  237.00    238.00    239.00    240.00    241.00    242.00    243.00
  244.00    245.00    246.00    247.00    248.00    249.00    250.00
  251.00    252.00    253.00    254.00    255.00    256.00    257.00
  258.00    259.00    260.00    261.00    262.00    263.00    264.00
  265.00    266.00    267.00    268.00    269.00    270.00    271.00
  272.00    273.00    274.00    275.00    276.00    277.00    278.00
  279.00    280.00    281.00    282.00    283.00    284.00    285.00
  286.00    287.00    288.00    289.00    290.00    291.00    292.00
  293.00    294.00    295.00    296.00    297.00    298.00    299.00
  300.00    301.00    302.00    303.00    304.00    305.00    306.00
  307.00    308.00    309.00    310.00    311.00    312.00    313.00
  314.00    315.00    316.00    317.00    318.00    319.00    320.00
  321.00    322.00    323.00    324.00    325.00    326.00    327.00
  328.00    329.00    330.00    331.00    332.00    333.00    334.00
  335.00    336.00    337.00    338.00    339.00    340.00    341.00
  342.00    343.00    344.00    345.00    346.00    347.00    348.00
  349.00    350.00    351.00    352.00    353.00    354.00    355.00
  356.00    357.00    358.00    359.00" |> x-> split(x,"\n") |> x-> join(x," ") |>
  x-> split(x," ") |> x-> filter(y-> y≠"",x) |>
  x-> [parse(Float64,y) for y in x]; #.+90.0

  # renormalization
  majoraxis_earth=majoraxis(ellipsoid(WGS84Latest)) |> x-> uconvert(km,x) |> ustrip;


  hᵢ="0.000000000000000E+000   3.60000000000000        3.80000000000000
  4.00000000000000        4.20000000000000        4.40000000000000
  4.60000000000000        4.80000000000000        5.00000000000000
  5.20000000000000        5.40000000000000        5.60000000000000
  5.80000000000000        6.00000000000000        6.20000000000000
  6.40000000000000        6.60000000000000        6.80000000000000
  7.00000000000000        7.20000000000000        7.40000000000000
  7.60000000000000        7.80000000000000        8.00000000000000
  8.20000000000000        8.40000000000000        8.60000000000000
  8.80000000000000        9.00000000000000        9.20000000000000
  9.40000000000000        9.60000000000000        9.80000000000000
  10.0000000000000        10.2000000000000        10.4000000000000
  10.6000000000000        10.8000000000000        11.0000000000000
  11.2000000000000        11.4000000000000        11.6000000000000
  11.8000000000000        12.0000000000000        12.2000000000000
  12.4000000000000        12.6000000000000        12.8000000000000
  13.0000000000000        13.2000000000000        13.4000000000000
  13.6000000000000        13.8000000000000        14.0000000000000
  14.2000000000000        14.4000000000000        14.6000000000000
  14.8000000000000        15.0000000000000        15.2000000000000
  15.4000000000000        15.6000000000000        15.8000000000000
  16.0000000000000        16.2000000000000        16.4000000000000
  16.6000000000000        16.8000000000000        17.0000000000000
  17.2000000000000        17.4000000000000        17.6000000000000
  17.8000000000000        18.0000000000000        18.2000000000000
  18.4000000000000        18.6000000000000        18.8000000000000
  19.0000000000000        19.2000000000000        19.4000000000000
  19.6000000000000        19.8000000000000        20.0000000000000
  21.0000000000000        22.0000000000000        23.0000000000000
  24.0000000000000        25.0000000000000        29.0000000000000
  37.0000000000000        44.0000000000000        68.0000000000000
  68.1500000000000    "  |> x-> split(x,"\n") |> x-> join(x," ") |>
  x-> split(x," ") |> x-> filter(y-> y≠"",x) |>
  x-> [parse(Float64,y) for y in x]./majoraxis_earth;
  # Generate discretized atmosphere using the knots

  #model=Carlotti()
  #interpolation=LinearPressure()


  for model in [Carlotti(), NoAtmosphere(), Ciddor(),Mathar4()]
    for interpolation in [LinearPressure(),LogarithmicPressure()]
      @testset "check $model with interpolation $interpolation" begin
        temperature_knots=copy(atmosphere.temperature)
        pressure_knots=copy(atmosphere.pressure)
        h_knots=copy(atmosphere.h)[:,1]
        θ_knots=copy(atmosphere.θ)[1,:]

        pressure_lowallocation,temperature_lowallocation,refractive_index_lowallocation=GeoUtils._discretize_atmosphere_lowallocation(
        view(pressure_knots,:,:),
        view(temperature_knots,:,:),
        view(h_knots,:),view(θ_knots,:),view(hᵢ,:),view(θᵢ,:);model=model,interpolation_pressure=interpolation)

        pressure,temperature,refractive_index=GeoUtils._discretize_atmosphere(view(atmosphere,:,:),
        view(hᵢ,:),view(θᵢ,:);model=model,interpolation_pressure=interpolation)

        @test pressure_lowallocation == pressure
        @test temperature_lowallocation == temperature
        @test refractive_index_lowallocation == refractive_index
      end
    end
  end


end
