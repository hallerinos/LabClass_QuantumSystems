using SparseArrays, LinearAlgebra

⊗(A,B) = kron(A, B)

function generate_lattice_operators(ô, N::Int64)
    d = size(ô["𝕀"],1)

    Ô = Dict()

    ρ = ô["n̂"]
    N̂ = spzeros(Float64, d^N, d^N)
    for n=1:N, k in keys(ô)
        idL = sparse(I, d^(n-1), d^(n-1))
        idR = sparse(I, d^(N-n), d^(N-n))
        mat = idL ⊗ ô[k] ⊗ idR
        Ô[k, n] = mat

        N̂ += idL ⊗ ρ ⊗ idR
    end
    return Ô, N̂
end