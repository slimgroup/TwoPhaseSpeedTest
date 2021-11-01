using PyPlot, JLD2, Polynomials, PyCall

include("../utils/Simulation.jl")
JLD2.@load "../TakeCompassSubsection/256K.jld2"

K = K[1:16:end,1:16:end,1:16:end]

phi = zeros(Float32, size(K))
for i = 1:16
    for j = 1:16
        for k = 1:16
            p = Polynomial([-0.0314^2*K[i,j,k],2*0.0314^2*K[i,j,k],-0.0314^2*K[i,j,k],1.527^2])
            phi[i,j,k] = minimum(real(roots(p)[findall(real(roots(p)).== roots(p))]))
        end
    end
end

n = (16, 16, 16)
d = (400f0, 400f0, 256f0)
qinj = (d[1]*n[1]/2, d[2]*n[2]/2, d[3]*n[3]/2)
qrate = 7
time = 2
dt = 1

TwoPhase(K, phi, qinj, qrate, d, time, dt; o=(0f0,0f0,(341-256)*6f0))
