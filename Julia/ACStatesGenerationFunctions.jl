# Coherent spin state in the Dicke basis
function Coherentdk(N::Int64, θ, φ)
    sinθ = sin(0.5*θ)
    cosθ_expϕ = cos(0.5*θ)*exp(-im*φ)
    if N < 67
        return [(1.0 + 0.0*im) * sqrt(binomial(N,k)) * (sinθ)^k * (cosθ_expϕ)^(N-k) for k in 0:N]
    else
        return [(1.0 + 0.0*im) * sqrt(binomial(BigInt(N),BigInt(k))) * (sinθ)^k * (cosθ_expϕ)^(N-k) for k in 0:N]
    end
end

# Angular momentum operators in the Dicke basis
function colspinjm(k::String,N::Int64)::Array{Complex{Float64},2}
    j = N/2.0;
    if N>=1
        if k=="0"
            res = Matrix{Complex{Float64}}(I, N+1, N+1);
        elseif k=="2"
            res = j * (j + 1.0) .* Matrix{Complex{Float64}}(I, N+1, N+1);
        elseif k=="x"
            res = 0.5 * (colspinjm("p",N) + colspinjm("m",N));
        elseif k=="y"
            res = -0.5 * im * (colspinjm("p",N) - colspinjm("m",N));
        elseif k=="z"
            res = diagm(0 => [1.0*m+0.0*im for m=-j:j]);
        elseif k=="m"
            res = diagm(1 => [(j * (j + 1.0) - m * (m - 1.0))^0.5 + 0.0*im for m=-j+1:j]);
        elseif k=="p"
            res = diagm(-1 => [(j * (j + 1.0) - m * (m + 1.0))^0.5 + 0.0*im for m=-j:j-1]);
        else
            throw(ArgumentError("Error OQLiege: undefined collective spin operator passed to colspinjm"))
            exit(code=1)
        end
        return res
    else
        throw(ArgumentError("Error OQLiege: N<1 passed to colspinjm"))
        exit(code=1)
    end
end

# Compute the Schmidt decomposition of a symmetric multiqubit system (i.e. single spin system)
function schmidt_coefficients_sym(ψ, t::Int)
    N = length(ψ) - 1

    ψ_Mat = [sum(
        k -> begin
            l = i - 1
            m = j - 1
            if l + m == k && 0 ≤ l ≤ t && 0 ≤ m ≤ N-t
                norm = sqrt(binomial(t, l) * binomial(N-t, m) / binomial(N, k))
                return ψ[k + 1] * norm
            else
                return 0.0 + 0.0im
            end
        end,
        0:N
    ) for i in 1:(t+1), j in 1:(N-t+1)]

    return svdvals(ψ_Mat)
end

# Compute the negativity of a symmetric multiqubit state for the bipartition t|N-t
function negativitypuresymqubit(ψ, t::Int)
    α = schmidt_coefficients_sym(ψ, t)
    return sum(i -> sum(j -> α[i] * α[j], 1:i-1), 2:t+1)
end

# Compute the Bures anticoherence measure of a spin state
function acBur(ψ,t::Int)
    neg = negativitypuresymqubit(ψ,t)
    return 1-√((abs(√(t+1)-√(1+2neg)))/(√(t+1)-1))
end

# Compute the exponential e^{i*J*θ} from the eigendecomposition (λ,V) of the J operator
function expJy(θ,λ,V)
    return V*diagm([exp(-im*θ*λ[i]) for i in axes(V,1)])*V'
end

# Compute the exponential e^{i*Jz^2*η} from the eigenvalues of the diagonal Jz^2 operator
function expJz2(η,Jz2)
    diagm([exp(-im*η*Jz2[i,i]) for i in axes(Jz2,1)])
end

# Apply a number of cycles of the protocol based on the parameters ps=(θ,η) and return the final anticoherence measure and the final state
function acStates_SqueezingRotation(t,λy,Vy,Jz2,ψ,ps;acmeasure=acBur)
    θs = vcat(0.0,ps[1:floor(Int,length(ps)/2)])
    χs = ps[ceil(Int,length(ps)/2):end]
    for i in 1:ceil(Int,length(ps)/2)
        ψ = expJz2(χs[i],Jz2)*expJy(θs[i],λy,Vy)*ψ
    end
    return 1-acmeasure(ψ,t),ψ
end