using Random
using SparseArrays
using KrylovKit

include("pow.jl")
include("make_figs.jl")
Ms = [2^i for i=1:10]
t_KKs, t_pows = zeros(size(Ms)), zeros(size(Ms))
n_KKs, n_pows = zeros(size(Ms)), zeros(size(Ms))

for dist in [rand, randn]
    for (id, M) in enumerate(Ms)
        A = sparse(dist(M,M))
        A = A + A'
        t_pows[id] = @elapsed λ_pow, b_pow, info = pow(A; outputlevel=1)
        n_pows[id] = info["numiter"]
        t_KKs[id] = @elapsed λ_KK, b_KK, info = eigsolve(A, dist(size(A,2)), 1, :LM)
        n_KKs[id] = info.numiter
    end

    make_scatter_line(Ms, t_pows, t_KKs, "power method","Lanczos (KrylovKit.jl)", "matrix dimension", "time [s]", "log", "log", "$(dist)_pow_vs_lanczos_time.png")

    make_scatter_line(Ms, n_pows, n_KKs, "power method","Lanczos (KrylovKit.jl)", "matrix dimension", "number of iterations", "log", "log", "$(dist)_pow_vs_lanczos_iter.png")
end