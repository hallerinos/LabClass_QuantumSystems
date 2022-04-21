using Random
using SparseArrays
using KrylovKit
using PyPlot
pygui(true)

include("pow.jl")
let
    Ms = [2^i for i=1:10]
    t_KKs, t_pows = zeros(size(Ms)), zeros(size(Ms))
    for (id, M) in enumerate(Ms)
        A = sparse(randn(M,M))
        A = A + A'
        t_pows[id] = @elapsed λ_pow, b_pow = pow(A; outputlevel=0)
        t_KKs[id] = @elapsed λ_KK, b_KK = eigsolve(A, randn(size(A,2)), 1, :LM)
    end
    plt.scatter(Ms, t_pows)
    plt.plot(Ms, t_pows, label="power method")
    plt.scatter(Ms, t_KKs)
    plt.plot(Ms, t_KKs, label="Lanczos (KrylovKit.jl)")
    plt.xlabel("matrix dimension")
    plt.ylabel("time [s]")
    plt.legend()
    plt.yscale("log")
end