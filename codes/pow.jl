using LinearAlgebra, Random, Printf
# solve eigenvalue problem via power method
function pow(A; b0=sparse(rand(size(A,2))), nmax=10000, outputlevel=0)
    ϵ = 1.0
    λ = 0.0
    bₖ = b0
    numiter = 0
    for k=0:nmax-1
        numiter += 1
        if outputlevel > 0 
            @printf "Iteration | Eigenvalue | Error : %5i | %.6g | %.6g\n" k λ ϵ
        end
        bₖ₊₁ = A*bₖ
        λ = bₖ'*bₖ₊₁
        bₖ₊₁ /= norm(bₖ₊₁)
        ϵ = norm(bₖ₊₁ - sign(λ)*bₖ)
        bₖ = bₖ₊₁
        if ϵ < eps()
            info = Dict("numiter" => numiter)
            return λ, bₖ, info
        end
    end
    info = Dict("numiter" => numiter)
    return λ, bₖ, info
end