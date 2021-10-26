using JLD2, JUDI

JLD2.@load "256cube_v.jld2" v
idx_ucfmt = find_water_bottom((v.-3500f0).*(v.>3500f0))

K = zeros(Float32, size(v))

for i = 1:256
    println("i=",i)
    for j = 1:256
        K[:,:,1:idx_ucfmt[i,j]-1] = 1.03*1f-3*v[:,:,1:idx_ucfmt[i,j]-1].+15f0
        K[:,:,idx_ucfmt[i,j]:idx_ucfmt[i,j]+8] .= 1f-4
        K[:,:,idx_ucfmt[i,j]+9:end] = 1.03*1f-3*v[:,:,idx_ucfmt[i,j]+9:end].+200f0
    end
end

JLD2.@save "256K.jld2" K