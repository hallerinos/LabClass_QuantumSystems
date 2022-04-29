using SparseArrays

function generate_spin_matrices(;dim=2)

    # initialize spin matrices
    Sx = zeros(ComplexF64,dim,dim);
    Sy = zeros(ComplexF64,dim,dim);
    Sz = zeros(ComplexF64,dim,dim);

    spin = 0.5*(dim-1)

    # construct Sx, Sy and Sz
    for idx = 1 : dim

        for idy = 1 : dim

            entryXY = 0.5 * sqrt((spin + 1) * (idx + idy - 1) - idx * idy);

            if (idx + 1) == idy
                Sx[idx,idy] += entryXY;
                Sy[idx,idy] -= 1im * entryXY;
            end

            if idx == (idy + 1)
                Sx[idx,idy] += entryXY;
                Sy[idx,idy] += 1im * entryXY;
            end

            if idx == idy
                Sz[idx,idy] += spin + 1 - idx;
            end

        end

    end

    # compute Id, Sp and Sm
    Id = one(Sz);
    Sp = Sx + 1im * Sy;
    Sm = Sx - 1im * Sy;
    ρ = Sp*Sm

    spin_matrices = Dict("Sx" => sparse(Sx), "Sy" => sparse(Sy), "Sz" => sparse(Sz), "S+" => sparse(Sp), "S-" => sparse(Sm), "ρ" => sparse(ρ), "Id" => sparse(Id))

    # return spin matrices
    return spin_matrices

end