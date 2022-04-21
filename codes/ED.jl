using LinearAlgebra
using SparseArrays
using KrylovKit
using PyPlot
pygui(true)

σ0 = sparse([1 0; 0  1])
σx = sparse([0 1; 1  0])
σz = sparse([1 0; 0 -1])

function generate_ops(σx, σz, N)
    d = size(σ0,1)
    Sx = Vector{SparseMatrixCSC{Float64, Int64}}(undef, N)
    Sz = Vector{SparseMatrixCSC{Float64, Int64}}(undef, N)
    for n=1:N
        Dn = d^(n-1)
        DN = d^(N-n)
        Sx[n] = 0.5*kron(sparse(I,Dn,Dn), σx, sparse(I,DN,DN))
        Sz[n] = 0.5*kron(sparse(I,Dn,Dn), σz, sparse(I,DN,DN))
    end
    return Sx, Sz
end

function generate_H(Sx, Sz; J=4, h=0, pbc=false)
    H = spzeros(size(Sx[1])) 
    for n=1:length(Sx)-1
        H += J*Sz[n]*Sz[n+1] + h*Sx[n]
    end
    H += h*Sx[end]
    if pbc
        H += J*Sz[1]*Sz[end]
    end
    return H
end

let 
    Ns = [6 12 18]
    pbc = true
    tol = eps()
    nev = 1
    kdmin = 10
    Js = [4]
    hs = 10 .^ (-2:0.2:2)

    
    Mxs, Mzs, SxSxs, SzSzs = [zeros(length(Ns), length(Js), length(hs), nev) for i=1:4]
    for (idJ, J) in enumerate(Js), (idh, h) in enumerate(hs), (idN, N) in enumerate(Ns)
        @info N
        Sx, Sz = generate_ops(σx, σz, N)
        Ĥ = generate_H(Sx, Sz; J=J, h=-h, pbc=pbc)
        E, Ψ, info = eigsolve(Ĥ, nev, :SR, eltype(Ĥ), issymmetric=true, krylovdim=max(nev,kdmin), tol=tol)
        ord = sortperm(E)
        E, Ψ = E[ord][1:nev], Ψ[ord][1:nev]

        M̂x, M̂z = sum(Sx)/N, sum(Sz)/N
        Mxs[idN, idJ, idh, :] .= [ψ'*M̂x*ψ for ψ in Ψ]
        Mzs[idN, idJ, idh, :] .= [ψ'*M̂z*ψ for ψ in Ψ]
        SxSx = Sx[N÷2]*Sx[N÷2+1]
        SzSz = Sz[N÷2]*Sz[N÷2+1]
        SxSxs[idN, idJ, idh, :] .= [ψ'*SxSx*ψ for ψ in Ψ]
        SzSzs[idN, idJ, idh, :] .= [ψ'*SzSz*ψ for ψ in Ψ]
    end

    plt.figure()
    [plt.plot(hs, 2abs.(Mxs[N,1,:,n]), label="2|Mx|,N=$(Ns[N])", marker="o") for N in 1:length(Ns) for n in 1:nev]
    # [plt.plot(hs, 2abs.(Mzs[N,1,:,n]), label="2|Mz|,N=$(Ns[N])", marker=".") for N in 1:length(Ns) for n in 1:nev]
    [plt.plot(hs, 4abs.(SzSzs[N,1,:,n]), label="4|SzSz|, N=$(Ns[N])", marker="h") for N in 1:length(Ns) for n in 1:nev]
    [plt.plot(hs, 4abs.(SxSxs[N,1,:,n]), label="4|SxSx|, N=$(Ns[N])", marker="s") for N in 1:length(Ns) for n in 1:nev]
    plt.axvline(2, color="black", linestyle="-.")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("h")
    plt.ylabel("")
    plt.text(0.01, 0.1, "ordered phase")
    plt.text(5, 0.0075, "disordered phase")
    plt.legend()
    plt.savefig("ising_phases.png")
end