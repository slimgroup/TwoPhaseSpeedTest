using PyPlot, JLD2, Polynomials, PyCall, MAT
include("../utils/Simulation.jl")
include("../utils/plotting.jl")

vars = matread("../Sleipner/sleipnertopdz.mat")

dz = vec(Float32.(vars["dz"]))
tops = Float32.(vars["tops"])

d = (50f0, 50f0, 5f0)
n = (64, 64, 100)

K = zeros(Float32, n)
K[:, :, 1:50] .= 15f0
K[:, :, 51:60] .= 0.00140825f0
K[:, :, 61:end] .= 2500f0

tops = tops[:, 1:64]

phi = 0.35f0 * ones(Float32, n)
qinj = (d[1]*32, d[2]*32, d[3]*(n[3]-20)+tops[32,32])
qrate = Float32(1.46/1.38)
time = 40
nt = 10

sat, p = TwoPhase(K, phi, qinj, qrate, d, time, nt; TOPS=tops)

matwrite("threelayer3D.mat", Dict(
	"sat" => sat,
	"p" => p
); compress = true)