using PyPlot, JLD2, Polynomials, PyCall, GLMakie, MAT
include("../utils/Simulation.jl")
include("../utils/plotting.jl")

vars = matread("../Sleipner/sleipnertopdz.mat")

dz = vec(Float32.(vars["dz"]))
tops = Float32.(vars["tops"])

d = (50f0, size(tops, 2) * 50f0, minimum(dz))
n = (size(tops,1), Int(ceil(sum(dz)/minimum(dz))))

K = zeros(Float32, n)
K[:, 1:800] .= 1500f0
K[:, 801:1100] .= 1f-3
K[:,1101:end] .= 2500f0

phi = 0.35f0 * ones(Float32, n)
tops = tops[:,43]
qinj = (d[1]*32, d[3]*(n[2]-100)+tops[32])
qrate = Float32(1.46/1.38)
time = 43
nt = 10

sat, p = TwoPhase(K, phi, qinj, qrate, d, time, nt; TOPS=tops)
