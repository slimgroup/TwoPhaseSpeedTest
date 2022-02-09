function TwoPhase(Kx::Array{T,3}, ϕ::Array{T,3}, qinj::Tuple{T1, T1, T1}, qrate::Number, d::Tuple{T2, T2, T2}, time::Number, nt::Int; Ky=nothing, Kz=nothing, TOPS=nothing) where {T, T1, T2}
    "Kx permeability, ϕ porosity, qinj injection coordinate [m], qrate injection rate [Mt/y]"

    if isnothing(Ky)
        Ky = Kx
    end
    if isnothing(Kz)
        Kz = Kx
    end

    WriteTxtFile(Kx; name="PERMX")
    WriteTxtFile(Ky; name="PERMY")
    WriteTxtFile(Kz; name="PERMZ")
    WriteTxtFile(ϕ; name="PORO")

    n = size(Kx)
    WriteDATAFile(qinj, qrate, n, d, time, nt;
    TOPS=TOPS, PERMXtxt="PERMX.txt", PERMYtxt="PERMY.txt", PERMZtxt="PERMZ.txt", POROtxt="PORO.txt", name="CO2SIMULATION")
    run(`flow CO2SIMULATION.DATA`)
    sat, p = ReadResults(; name="CO2SIMULATION")

    return sat, p
end

function TwoPhase(Kx::Array{T,2}, ϕ::Array{T,2}, qinj::Tuple{T1, T1}, qrate::Number, d::Tuple{T2, T2, T2}, time::Number, nt::Int; Ky=nothing, Kz=nothing, TOPS=nothing) where {T, T1, T2}
    n = size(Kx)
    if isnothing(Ky)
        Ky = reshape(Kx,n[1],1,n[2])
    else
        Ky = reshape(Ky,n[1],1,n[2])
    end
    if isnothing(Kz)
        Kz = reshape(Kx,n[1],1,n[2])
    else
        Kz = reshape(Kz,n[1],1,n[2])
    end
    sat, p = TwoPhase(reshape(Kx,n[1],1,n[2]), reshape(ϕ,n[1],1,n[2]), (qinj[1],T1(d[2]/2),qinj[2]), qrate, d, time, nt; Ky=Ky, Kz=Kz, TOPS=TOPS)
    sat = [dropdims(sat[i], dims=2) for i = 1:length(sat)]
    p = [dropdims(p[i], dims=2) for i = 1:length(p)]

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
    sat, p = py"readecl"(name)
    sat = permutedims.(sat, [[3,2,1] for i = 1:length(sat)])
    p = permutedims.(p, [[3,2,1] for i = 1:length(sat)])
    return sat, p
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
    K, ϕ = py"readeclinit"(name)
    K = permutedims(K, [3,2,1])
    ϕ = permutedims(ϕ, [3,2,1])
    return K, ϕ
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

function WriteTxtFile(var::Matrix{T}; name="PERMX") where T # assume X*Z
    WriteTxtFile(reshape(var,size(var,1),1,size(var,2)); name=name)
end

function WriteTxtFile(var::Vector{T}; name="PERMX") where T # assume Z
    WriteTxtFile(reshape(var,1,1,size(var,1),1,1); name=name)
end

function WriteTOPS(TOPS::T, n::Tuple{Int, Int, Int}) where T <: Number  # single value for the top
    TOPS1 = "TOPS
    $(n[1]*n[2])*$TOPS /"

    return TOPS1
end

function WriteTOPS(TOPS::Vector{T}, n::Tuple{Int, Int, Int}) where T <: Number    # vector top, assume 2D X*Z
    
    WriteTxtFile(TOPS; name="TOPS")
    TOPS1 = "INCLUDE
    TOPS.txt /"

    return TOPS1
end

function WriteTOPS(TOPS::Matrix{T}, n::Tuple{Int, Int, Int}) where T <: Number  # matrix top, assume X*Y*Z
    
    WriteTxtFile(TOPS; name="TOPS")
    TOPS1 = "INCLUDE
    TOPS.txt /"

    return TOPS1
end

function WriteTOPS(TOPS, n::Tuple{Int, Int, Int})   # matrix top, assume X*Y*Z
    
    isnothing(TOPS) || @warn "Type of TOPS not supported"
    return WriteTOPS(0, n)
end

function getqgrid(qinj::Tuple{T1, T1, T1}, d::Tuple{T3, T3, T3}, TOPS::T) where {T, T1, T3}
    qgridxy = Int.(ceil.(qinj[1:2]./d[1:2]))
    qgridz = Int(ceil((qinj[end]-TOPS)/d[3]))
    return (qgridxy[1], qgridxy[2], qgridz)
end

function getqgrid(qinj::Tuple{T1, T1, T1}, d::Tuple{T3, T3, T3}, TOPS::Vector{T}) where {T, T1, T3}
    qgridxy = Int.(ceil.(qinj[1:2]./d[1:2]))
    qgridz = Int(ceil((qinj[end]-TOPS[qgridxy[1]])/d[3]))
    return (qgridxy[1], qgridxy[2], qgridz)
end

function getqgrid(qinj::Tuple{T1, T1, T1}, d::Tuple{T3, T3, T3}, TOPS::Matrix{T}) where {T, T1, T3}
    qgridxy = Int.(ceil.(qinj[1:2]./d[1:2]))
    qgridz = Int(ceil((qinj[end]-TOPS[qgridxy[1], qgridxy[2]])/d[3]))
    return (qgridxy[1], qgridxy[2], qgridz)
end

function WriteDATAFile(qinj::Tuple{T1, T1, T1}, qrate::Number, n::Tuple{Int, Int, Int}, d::Tuple{T3, T3, T3}, time::Number, nt::Int;
    TOPS=nothing, PERMXtxt::String = "PERMX.txt", PERMYtxt::String = "PERMY.txt", PERMZtxt::String = "PERMZ.txt", POROtxt::String = "PORO.txt", name="CO2SIMULATION") where {T1, T2, T3}
    
    println("Writing $(name).txt")

    DIMENS = "DIMENS 
$(n[1]) $(n[2]) $(n[3]) /"

    DX = "DX 
    $(prod(n))*$(d[1]) /"

    DY = "DY  
    $(prod(n))*$(d[2]) /"
    
    DZ = "DZ 
    $(prod(n))*$(d[3]) /"

    TOPS1 = WriteTOPS(TOPS, n)

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

    qgrid = getqgrid(qinj, d, TOPS)

    rate = qrate * 1f6 * 1f3 / 1.98 / 365.25

    WELL = "WELSPECS
    Injector I $(qgrid[1]) $(qgrid[2]) 0.0e+00 WATER 0 STD SHUT NO 0 SEG 0/
    /
    
    COMPDAT
    Injector $(qgrid[1]) $(qgrid[2]) $(qgrid[3]) $(qgrid[3]) OPEN -1 8.5107288246779274e+02 2.0e-01 -1.0 0 1* Y -1.0/
    Injector $(qgrid[1]) $(qgrid[2]) $(qgrid[3]) $(qgrid[3]) OPEN -1 4.0622380088359796e+02 2.0e-01 -1.0 0 1* Y -1.0/
    /
    
    WCONINJE
    Injector GAS OPEN RATE $(rate) 3* 0 /
    /"

    TSTEP = "TSTEP
    $(nt)*$(time*365.25/nt)
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
CO2STORE
    
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

$(TOPS1)
    
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
    0.0 0.0
    $(n[end]*d[end]) 0.0
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