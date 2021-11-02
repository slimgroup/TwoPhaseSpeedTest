function TwoPhase(Kx::Array{T,3}, ϕ::Array{T,3}, qinj::Tuple{T1, T1, T1}, qrate::Number, d::Tuple{T2, T2, T2}, time::Number, dt::Number; Ky=nothing, Kz=nothing, o=nothing) where {T, T1, T2}
    "Kx permeability, ϕ porosity, qinj injection coordinate [m], qrate injection rate [Mt/y]"

    if isnothing(Ky)
        Ky = Kx
    end
    if isnothing(Kz)
        Kz = Kx
    end
    if isnothing(o)
        o = zeros(Float32, length(size(Kx)))
    end

    WriteTxtFile(Kx; name="PERMX")
    WriteTxtFile(Ky; name="PERMY")
    WriteTxtFile(Kz; name="PERMZ")
    WriteTxtFile(ϕ; name="PORO")

    n = size(Kx)
    WriteDATAFile(qinj, qrate, n, d, o, time, dt;
    PERMXtxt="PERMX.txt", PERMYtxt="PERMY.txt", PERMZtxt="PERMZ.txt", POROtxt="PORO.txt", name="CO2SIMULATION")
    run(`flow CO2SIMULATION.DATA`)
    sat, p = ReadResults(; name="CO2SIMULATION")

    return sat, p
end

function ReadResults(;name="CO2SIMULATION")
    py"""
    from ecl.grid import EclGrid, EclRegion
    from ecl.eclfile import EclFile, EclRestartFile
    import matplotlib.pyplot as plt
    import numpy as np
    ###################################################################################################
    # Read eclipse files
    def readecl(name):
        # Load grid
        grid = EclGrid(name + ".EGRID")
        shape = (grid.nx, grid.ny, grid.nz)[::-1]
        # Load snapshots
        rst_file = EclRestartFile(grid, name + ".UNRST")
        sat = [rst_file['SGAS'][i].numpy_view().reshape(shape) for i in range(len(rst_file['SGAS']))]
        p = [rst_file['PRESSURE'][i].numpy_view().reshape(shape) for i in range(len(rst_file['PRESSURE']))]
        return sat, p
    """
    return py"readecl"(name)
end

function ReadInitFile(;name="CO2SIMULATION")
    py"""
    from ecl.grid import EclGrid, EclRegion
    from ecl.eclfile import EclFile, EclRestartFile
    import matplotlib.pyplot as plt
    import numpy as np
    ###################################################################################################
    # Read eclipse files
    def readeclinit(name):
        # Load grid
        grid = EclGrid(name + ".EGRID")
        shape = (grid.nx, grid.ny, grid.nz)[::-1]
        init_file = EclFile(name + ".INIT")
        ϕ = init_file['PORO'][0].numpy_view().reshape(shape)
        K = init_file['PERMX'][0].numpy_view().reshape(shape)
        return K, ϕ
    """
    return py"readeclinit"(name)
end

function WriteTxtFile(var::Array{T,3}; name="PERMX") where T
    
    println("Writing $(name).txt")
    f = open(name*".txt","w")

    write(f, name)
    write(f, "\n")
    global ct = 0
    for k = 1:size(var,3)
        for j = 1:size(var,2)
            for i = 1:size(var,1)
                global ct = ct + 1
                write(f, string(var[i,j,k]))
                write(f, " ")
                (ct%5 == 0) && write(f, "\n")
            end
        end
    end

    write(f, " /")
    close(f)

end

function WriteTxtFile(var::Array{T,2}; name="PERMX.txt") where T
    WriteTxtFile(reshape(var,size(var,1),1,size(var,2)); name=name)
end

function WriteDATAFile(qinj::Tuple{T1, T1, T1}, qrate::Number, n::Tuple{Int, Int, Int}, d::Tuple{T3, T3, T3}, o::Tuple{T2, T2, T2}, time::Number, dt::Number;
    PERMXtxt::String = "PERMX.txt", PERMYtxt::String = "PERMY.txt", PERMZtxt::String = "PERMZ.txt", POROtxt::String = "PORO.txt", name="CO2SIMULATION") where {T1, T2, T3}
    
    println("Writing $(name).txt")
    DIMENS = "DIMENS 
$(n[1]) $(n[2]) $(n[3]) /"

    DX = "DX 
    $(prod(n))*$(d[1]) /"

    DY = "DY  
    $(prod(n))*$(d[2]) /"
    
    DZ = "DZ 
    $(prod(n))*$(d[3]) /"

    TOPS = "TOPS
    $(prod(n))*$(o[3]) /"

    PERMX = "INCLUDE
    '$(PERMXtxt)'
/"

    PERMY = "INCLUDE
    '$(PERMYtxt)'
/"

    PERMZ = "INCLUDE
    '$(PERMZtxt)'
/"

    PORO = "INCLUDE
    '$(POROtxt)'
/"

    qgrid = Int.(round.((qinj.-o)./d))

    rate = qrate * 1f6 * 1f3 / 1.98 / 365

    WELL = "WELSPECS
    Injector I $(qgrid[1]) $(qgrid[2]) 0.0e+00 WATER 0 STD SHUT NO 0 SEG 0/
    /
    
    COMPDAT
    Injector $(qgrid[1]) $(qgrid[2]) $(qgrid[3]) $(qgrid[3]) OPEN -1 8.5107288246779274e+02 2.0e-01 -1.0 0 1* Y -1.0/
    /
    
    WCONINJE
    Injector GAS OPEN RATE $(rate) 3* 0 /
    /"

    TSTEP = "TSTEP
    $(Int(round(365.25*time/dt)))*$(dt)
    /
    "

    str = 
    
"RUNSPEC
    
TITLE
    Simulation 
    
$(DIMENS)
    
EQLDIMS
    1 100 50 1 50
/
    
TABDIMS
    1 1 40 20 2 20
/
    
WELLDIMS
    1 1 1 1
/
    
OIL
GAS
CO2STOR
    
METRIC
UNIFOUT
START 
    1 'JAN' 2000
/
------------------------------------------------------
GRID
    
INIT

$(DX)

$(DY)

$(DZ)

$(TOPS)
    
$(PERMX)

$(PERMY)

$(PERMZ)

$(PORO)
    
------------------------------------------------------
------------------------------------------------------
PROPS
    
ROCK
    1.0e+01 1.6e-06
/
    
-- Brooks-Corey with coefficient 2.8 and
-- entry pressure 2.5kPa (Cavanagh. A. 2013)
SGOF
-- Column 1: gas saturation
-- Column 2: gas relative permeability
-- Column 3: oil relative permeability when oil, gas and connate water are present
-- Column 4: oil-gas capillary pressure
    0.0	0.0	1.0 	0.025
    0.1	0.0     0.740	0.026
    0.2	0.009	0.528	0.027
    0.3	0.030	0.359	0.028
    0.4	0.070	0.230	0.030
    0.5	0.136	0.136	0.032
    0.6	0.230	0.070	0.035
    0.7	0.359	0.030	0.038
    0.8	0.528	0.009	0.044
    0.9	0.740	0.000	0.057 /
    
SALINITY
    0.7 / --35-40g/l  -> 35-40g /kg -> 0.63-0.72 mol/g
     
---------------------------------------------------------
REGIONS
---------------------------------------------------------
    
------------------------------------------------------
SOLUTION
    
EQUIL
    8.6300e+02 1.465583e+02  5.050e+03 0.0 1.0e+02 0.0  1 0 0 / 
/
    
RSVD
    $(o[end]) 0.0
    $(n[end]*d[end]+o[end]) 0.0
/
    
RTEMPVD
    800 32.7134 
    1000 41
/
------------------------------------------------------
SUMMARY
    
-----------------------------------------------------
SCHEDULE
RPTSCHED
PRES SGAS RS WELLS
/
RPTRST
BASIC=1
/
    
$(WELL)

$(TSTEP)

END
"


    f = open(name*".DATA","w")

    write(f, str)
    close(f)
end