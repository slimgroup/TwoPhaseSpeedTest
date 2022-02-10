using PyPlot, JLD2, Polynomials, PyCall, GLMakie, Images, MAT
include("../utils/Simulation.jl")
include("../utils/plotting.jl")
include("../utils/grid.jl")
sat, p = ReadResults(;name="SLEIPNER_ORG")

vars = matread("sleipnertopdz.mat")
dz = vec(Float32.(vars["dz"]))
tops = Float32.(vars["tops"])

sat1 = grid_change(sat, dz, minimum(dz))
p1 = grid_change(p, dz, minimum(dz))

extent = (0, 50*64, cumsum(dz)[end]+tops[31, 43], tops[31, 43])

qinj = (50*31, cumsum(dz)[214]+tops[31, 43])

figure();
imshow(sat1[end][:,43,:]',extent=extent, aspect="auto");
xlabel("X [m]")
ylabel("Z [m]")
colorbar()
PyPlot.scatter(qinj[1], qinj[2], label="injection")
legend()
title("CO2 saturation in original 3D sleipner, take a 2D slice")
savefig("3dsleipner.png", bbox_inches="tight", dpi=300)
