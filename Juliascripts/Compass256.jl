using PyPlot, JLD2, Polynomials, PyCall, GLMakie, Images
include("../utils/Simulation.jl")
include("../utils/plotting.jl")

n = (256, 256, 256)
JLD2.@load "256K.jld2"

#K = Float32.(imfilter(K, Kernel.gaussian((2,2,2))))

phi = 0.36f0*ones(Float32, n)

for i = 1:n[1]
    for j = 1:n[2]
        for k = 1:n[3]
            p = Polynomial([-0.0314^2*K[i,j,k],2*0.0314^2*K[i,j,k],-0.0314^2*K[i,j,k],1.527^2])
            phi[i,j,k] = minimum(real(roots(p)[findall(real(roots(p)).== roots(p))]))
        end
    end
end

d = (25f0, 25f0, 6f0)
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
