function grid_change(x::Array{T,3}, dz::Vector{T1}, d::T1) where {T, T1}
    # this happens on the 3rd dimension
    @assert size(x,3) == length(dz)
    stretch = Int.(round.(dz./d))
    y = zeros(T, size(x,1), size(x,2), sum(stretch))
    ct = 1
    for i = 1:size(x,3)
        copyto!(view(y, :, :,ct:ct+stretch[i]-1), repeat(view(x, :, :, i), 1, 1, stretch[i])) 
        ct = ct + stretch[i]
    end
    return y
end

function grid_change(x::Matrix{T}, dz::Vector{T1}, d::T1) where {T, T1}
    return grid_change(reshape(x, size(x,1), 1, size(x,2)), dz, d)[:,1,:]
end

function grid_change(x::Vector{T}, dz::Vector{T1}, d::T1) where {T<:Number, T1}
    return grid_change(reshape(x, 1, 1, size(x,1)), dz, d)[1,1,:]
end

const Param{T} = Union{Vector{T}, Matrix{T}, Array{T,3}}

function grid_change(x::Vector{T}, dz::Vector{T1}, d::T1) where {T<:Param, T1}
    return [grid_change(x[i], dz, d) for i = 1:length(x)]
end
