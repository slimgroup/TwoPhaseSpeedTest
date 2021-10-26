using JLD2
using PyPlot
using SegyIO

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
JLD2.@save "256cube_v.jld2" v