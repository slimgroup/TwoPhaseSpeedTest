using PyPlot, JLD2, Polynomials, PyCall, GLMakie, Images, MAT
include("../utils/Simulation.jl")
include("../utils/plotting.jl")

K, ϕ = ReadInitFile(; name="SLEIPNER_ORG")
n = size(K)

vars = matread("sleipnertopdz.mat")

dz = Float32.(vars["dz"])
tops = Float32.(vars["tops"])

d = (50f0, 50f0, minimum(dz))
qinj = (d[1]*31/2, d[2]*43/2,d[3]*214+tops[31,43])
qrate = 1
time = 15
nt = 100

sat, p = TwoPhase(K, ϕ, qinj, qrate, d, time, nt; TOPS=tops)
