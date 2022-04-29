using KrylovKit, DataFrames, PyPlot
include("spinN.jl")
include("generate_lattice_ops.jl")

L = 5
t = 1.
Us = LinRange(0, 10, 11)
Us = [10]
μs = LinRange(0, 10, 11)
tol = eps()
nev = 1
kdmin = 10

Ô, N̂ = generate_spin_ops(L; dim=6)

obs = DataFrame(U=Float64[], μ=Float64[], n=Int[], ev=Int[], obs=String[], val=ComplexF64[])
for U in Us, μ in μs
    @info "U: $U, μ:$μ"
    Ĥ = convert(SparseMatrixCSC{Float64, Int64}, zero(Ô["Id",1]))
    for n=1:L
        Ĥ += -t*Ô["S+",n]*Ô["S-",mod(n,L)+1]
        Ĥ += -t*Ô["S+",mod(n,L)+1]*Ô["S-",n]
        Ĥ += -μ*Ô["ρ",n] + 0.5*U*Ô["ρ",n]*(Ô["ρ",n]-Ô["Id",n])
    end

    E, Ψ, info = eigsolve(Ĥ, nev, :SR, eltype(Ĥ), issymmetric=true, krylovdim=max(nev,kdmin), tol=tol)
    ord = sortperm(E)
    E, Ψ = E[ord][1:nev], Ψ[ord][1:nev]

    ops = ["S+", "S-", "ρ"]
    ops = ["ρ"]
    for n=1:L
        for o in ops
            [push!(obs, (U, μ, n, id, o, ψ'*Ô[o, n]*ψ)) for (id,ψ) in enumerate(Ψ)]
        end
    end
end