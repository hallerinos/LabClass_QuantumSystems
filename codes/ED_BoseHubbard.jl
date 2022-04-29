using KrylovKit, DataFrames, PyPlot, CSV
include("bosonic_ops.jl")
include("generate_lattice_ops.jl")

L = 6
d = 4
ts = LinRange(0, 0.5, 21)
# ts = [1]
pbc = false
Us = [1]
μs = LinRange(0, 2, 21)
# μs = [0]
tol = eps()
nev = 4
kdmin = 10

ô = generate_boson_matrices(dim=d)
Ô, N̂ = generate_lattice_operators(ô, L)

obs = DataFrame(u=Float64[], mu=Float64[], t=Float64[], i=Int[], ev=Int[], obs=String[], val_re=Float64[], val_im=Float64[])
for U in Us, μ in μs, t in ts
    Ĥ = convert(SparseMatrixCSC{Float64, Int64}, zero(Ô["𝕀",1]))
    for n=1:L-1
        Ĥ += -t*Ô["â⁺",n]*Ô["â",n+1]
        Ĥ += -t*Ô["â⁺",n+1]*Ô["â",n]
        Ĥ += -μ*Ô["n̂",n] + 0.5*U*Ô["n̂",n]*(Ô["n̂",n]-Ô["𝕀",n])
    end
    Ĥ += -μ*Ô["n̂",L] + 0.5*U*Ô["n̂",L]*(Ô["n̂",L]-Ô["𝕀",L])
    if pbc
        Ĥ += -t*Ô["â⁺",L]*Ô["â",1]
        Ĥ += -t*Ô["â⁺",1]*Ô["â",L]
    end

    E, Ψ, info = eigsolve(Ĥ, nev, :SR, eltype(Ĥ), issymmetric=true, krylovdim=max(nev,kdmin), tol=tol)
    ord = sortperm(E)
    E, Ψ = E[ord][1:nev], Ψ[ord][1:nev]

    @info "U: $U, μ:$μ, t:$t, E:$E"

    ops = ["â⁺", "â", "n̂"]
    # ops = ["ρ"]
    for i=1:L
        for o in ops
            vals = [ψ'*Ô[o, i]*ψ for (id,ψ) in enumerate(Ψ)]
            [push!(obs, (U, μ, t, i, id, o, real(v), imag(v))) for (id,v) in enumerate(vals)]
        end
    end
end

CSV.write("observables.csv", obs)