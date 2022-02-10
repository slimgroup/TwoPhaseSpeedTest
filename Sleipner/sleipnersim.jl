### compress 3D sleipner to 2D

using PyPlot, JLD2, Polynomials, PyCall, GLMakie, Images, MAT, Statistics
include("../utils/Simulation.jl")
include("../utils/plotting.jl")

K, ϕ = ReadInitFile(; name="SLEIPNER_ORG")
n = size(K)

vars = matread("sleipnertopdz.mat")

dz = Float32.(vars["dz"])
tops = Float32.(vars["tops"])

K = K[:,1,:]
ϕ = ϕ[:,1,:]
d = (50f0, size(tops, 2) * 50f0, vec(dz))
tops = tops[:,43]
qinj = (d[1]*31, cumsum(d[3])[214]+tops[31])
qrate = Float32(1.46/1.38)
time = 15
nt = 100
sat, p = TwoPhase(K, ϕ, qinj, qrate, d, time, nt; TOPS = tops)
