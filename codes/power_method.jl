using Random
using SparseArrays
using KrylovKit

include("pow.jl")
let
    M = 100
    A = sparse(randn(M,M))
    A = A + A'
    @time λ_pow, b_pow = pow(A)
    @time λ_KK, b_KK = eigsolve(A, randn(size(A,2)), 1, :LM)
    print(λ_pow - λ_KK[1])
end