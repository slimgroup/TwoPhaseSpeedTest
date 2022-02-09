using PyPlot, JLD2, Polynomials, PyCall, GLMakie, Images
include("../utils/Simulation.jl")
include("../utils/plotting.jl")

n = (256, 256)
d = (0.5f0, 1000f0, 0.5f0)
K = zeros(Float32, n)
K[:,1:183] .= 15f0
K[:,184:190] .= 1f-3
K[:,191:end] .= 2500f0

phi = 0.35f0 * ones(Float32, n)

qinj = (d[1]*n[1]/2f0, d[3]*n[2]+(341-256)*6f0-25f0)
qrate = 7
time = 100
nt = 10

sat, p = TwoPhase(K, phi, qinj, qrate, d, time, nt; TOPS= (341-256)*6f0)
