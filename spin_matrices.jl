using LinearAlgebra, Random, SparseArrays

N = 12
Sx = Vector{SparseMatrixCSC{Float64, Int64}}(undef, N)
A,B = [rand(10,10) for i=1:2]
kron(A,B)