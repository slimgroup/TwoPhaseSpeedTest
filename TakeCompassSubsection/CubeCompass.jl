using JLD2
using PyPlot
using SegyIO
using JUDI

run(`wget -r ftp://slim.gatech.edu/data//synth/Compass/final_velocity_model_ieee_6m.sgy`) # this might take a while

# get velocity
block = segy_read("slim.gatech.edu/data/synth/Compass/final_velocity_model_ieee_6m.sgy")

# original compass model is in 25m*25m*6m
n1 = (1911,2730,341)
d = (25f0,25f0,6f0)

n = (637,910,341)

sx = get_header(block, "SourceX")
sy = get_header(block, "SourceY")

v_nogas = zeros(Float32,n)

for i = 1:n[1]
    x = d[1].*(i-1)
    inds = findall(sx.==x)
    slice = block.data[:,inds[sortperm(sy[inds])]]

    v_nogas[i,:,:] = transpose(slice[:,1:Int(end/3)])
end

v = v_nogas[1:256,1:256,end-255:end]
idx_ucfmt = find_water_bottom((v.-3500f0).*(v.>3500f0))

K = zeros(Float32, size(v))

for i = 1:256
    for j = 1:256
        K[:,:,1:idx_ucfmt[i,j]-1] = 1.03*1f-3*v[:,:,1:idx_ucfmt[i,j]-1].+15f0
        K[:,:,idx_ucfmt[i,j]:idx_ucfmt[i,j]+8] .= 1f-4
        K[:,:,idx_ucfmt[i,j]+9:end] = 1.03*1f-3*v[:,:,idx_ucfmt[i,j]+9:end].+200f0
    end
end

phi = zeros(Float32,n)
for i = 1:n[1]
    for j = 1:n[2]
        p = Polynomial([-0.0314^2*Kh[i,j],2*0.0314^2*Kh[i,j],-0.0314^2*Kh[i,j],1.527^2])
        phi[i,j] = minimum(real(roots(p)[findall(real(roots(p)).== roots(p))]))
    end
    for j = idx_ucfmt[i]:idx_ucfmt[i]+3
        phi[i,idx_ucfmt[i]:idx_ucfmt[i]+3] = Float32.(range(0.056,stop=0.1,length=4))
    end
end

JLD2.@save "../model/CompassCube.jld2" v

