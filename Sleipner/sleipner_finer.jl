### compress 3D sleipner to 2D

using PyPlot, JLD2, Polynomials, PyCall, GLMakie, Images, MAT, Statistics
include("../utils/Simulation.jl")
include("../utils/plotting.jl")
include("../utils/grid.jl")
K, ϕ = ReadInitFile(; name="SLEIPNER_ORG")
n = size(K)

vars = matread("sleipnertopdz.mat")

dz = vec(Float32.(vars["dz"]))
tops = Float32.(vars["tops"])

K = grid_change(K[:,1,:], dz, minimum(dz))
ϕ = grid_change(ϕ[:,1,:], dz, minimum(dz))

d = (50f0, size(tops, 2) * 50f0, minimum(dz))
tops = tops[:,43]
qinj = (d[1]*31, cumsum(dz)[214]+tops[31])
qrate = Float32(1.46/1.38)
time = 15
nt = 100
sat, p = TwoPhase(K, ϕ, qinj, qrate, d, time, nt; TOPS = tops)

extent = (0, d[1]*n[1], d[3]*size(K,2)+tops[31], tops[31])

figure();
imshow(sat[end]',extent=extent, aspect="auto");
xlabel("X [m]")
ylabel("Z [m]")
colorbar()
PyPlot.scatter(qinj[1], qinj[2], label="injection")
legend()
title("CO2 saturation, uniform grid")
savefig("2dsleipneruniform.png", bbox_inches="tight", dpi=300)

matwrite("sleipner2dfine.mat", Dict(
	"sat" => sat,
	"p" => p
); compress = true)