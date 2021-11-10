using PyPlot, JLD2, Polynomials, PyCall, Images
include("../utils/Simulation.jl")
include("../utils/plotting.jl")

n = (128, 128, 128)
JLD2.@load "256Kphi.jld2"

K = K[1:2:end,1:2:end,1:2:end]
phi = phi[1:2:end,1:2:end,1:2:end]

#K = Float32.(imfilter(K, Kernel.gaussian((2,2,2))
d = (50f0, 50f0, 12f0)
qinj = (d[1]*n[1]/2, d[2]*n[2]/2, d[3]*n[3]+(341-256)*6f0-300f0)
qrate = 7
time = 80
nt = 10

sat, p = TwoPhase(K, phi, qinj, qrate, d, time, nt; o=(0f0,0f0,(341-256)*6f0))

JLD2.@save "satpCompass256.jld2" sat p

#figure();imshow(K[:,32,:]')
#figure();imshow(sat[end][:,32,:]')

#plot_3D(reverse(sat[end],dims=3))
#plot_3D(reverse(K,dims=3))
#plot_3D(reverse(phi,dims=3))
