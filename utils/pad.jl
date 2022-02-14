function pad(x::Array{T, N}, nb::Tuple{Vararg{Int, N}}) where {T, N}
    @assert length(nb) == N
    n = size(x)
    n1 = (n[1:N-1] .+ 2 .*nb[1:N-1])...,n[N] + nb[N]
    A = joKron([joExtend(n[i], :border; pad_upper=nb[end], pad_lower=nb[i], DDT=T, RDT=T) for i = 1:N-1]... , joExtend(n[N], :border; pad_upper=nb[N], DDT=T, RDT=T))
    return reshape(A*vec(x), n1)
end