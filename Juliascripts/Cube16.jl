using PyPlot, JLD2, Polynomials, PyCall, GLMakie, Images
include("../utils/Simulation.jl")
include("../utils/plotting.jl")

n = (64, 64, 64)
K = zeros(Float32, n)
K[:,:,1:26] .= 15f0
K[:,:,27:38] .= 1f-2
K[:,:,39:end] .= 2500f0
phi = 0.36f0*ones(Float32, n)

d = (50f0, 50f0, 0.5f0)
qinj = (d[1]*n[1]/2, d[2]*n[2]/2, d[3]*n[3]+(341-256)*6f0-1f0)
qrate = 7
time = 2
nt = 10

TOPS = 10f0*randn(Float32, n[1], n[2]) .+ (341-256)*6f0

sat, p = TwoPhase(K, phi, qinj, qrate, d, time, nt; TOPS=TOPS)

figure();imshow(K[:,32,:]')
figure();imshow(sat[end][:,32,:]')


#plot_3D(reverse(sat[end],dims=3))
#plot_3D(reverse(K,dims=3))
#plot_3D(reverse(phi,dims=3))
