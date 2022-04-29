using SparseArrays

function generate_boson_matrices(;dim=2)

    # compute Id, Sp and Sm
    𝕀 = Matrix(I, dim, dim);
    â = Bidiagonal([0 for i=0:d-1], [sqrt(i) for i=1:d-1], :U)
    â⁺ = â'
    n̂ = â⁺*â

    bosonic_matrices = Dict("â⁺" => sparse(â⁺), "â" => sparse(â), "n̂" => sparse(n̂), "𝕀" => sparse(𝕀))

    # return matrices
    return bosonic_matrices

end