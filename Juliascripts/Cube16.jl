using PyPlot, JLD2, Polynomials, PyCall

include("../utils/Simulation.jl")
JLD2.@load "16Cube.jld2"

n = (16, 16, 16)
d = (400f0, 400f0, 256f0)
qinj = (d[1]*n[1]/2, d[2]*n[2]/2, d[3]*n[3]+(341-256)*6f0)
qrate = 7
time = 100
dt = 1

sat, p = TwoPhase(K, phi, qinj, qrate, d, time, dt; o=(0f0,0f0,(341-256)*6f0))
