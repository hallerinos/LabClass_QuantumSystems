using LinearAlgebra, Random, Printf
# solve eigenvalue problem via power method
function pow(A; b0=rand(size(A,2)), nmax=10000, outputlevel=0)
    ϵ = 1.0
    λ = 0.0
    bₖ = b0
    for k=0:nmax-1
        if outputlevel > 0 
            @printf "Iteration | Eigenvalue | Error : %5i | %.6g | %.6g\n" k λ ϵ
        end
        bₖ₊₁ = A*bₖ
        λ = bₖ'*bₖ₊₁
        bₖ₊₁ /= norm(bₖ₊₁)
        ϵ = norm(bₖ₊₁ - sign(λ)*bₖ)
        bₖ = bₖ₊₁
        if ϵ < eps()
            return λ, bₖ
        end
    end
    return λ, bₖ
end