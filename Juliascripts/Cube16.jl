using PyPlot, JLD2, Polynomials, PyCall, GLMakie
run(`alias flow=/Users/francisyin/Desktop/opm/opm-simulators/build/bin/flow`)
include("../utils/Simulation.jl")
JLD2.@load "16Cube.jld2"
include("../utils/plotting.jl")

n = (16, 16, 16)
d = (400f0, 400f0, 96f0)
qinj = (d[1]*n[1]/2, d[2]*n[2]/2, d[3]*n[3]+(341-256)*6f0)
qrate = 7
time = 30
nt = 1000

sat, p = TwoPhase(K, phi, qinj, qrate, d, time, nt; o=(0f0,0f0,(341-256)*6f0))

plot_3D(sat[end])