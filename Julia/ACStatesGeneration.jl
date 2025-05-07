using Optim
using DelimitedFiles
using LinearAlgebra

include("ACStatesGenerationFunctions.jl")

project_dir = splitdir(@__DIR__)[1]
data_dir = project_dir*"/Data/"

#######################################################################
#                 Example of protocol optimisation
#######################################################################

N = 6
j = N/2
Jx = colspinjm("x",N)
Jy = colspinjm("y",N)
Jz = colspinjm("z",N)
Jz2 = Jz^2
λy, Vy = eigen(Jy)
ψ = exp(-im*Jx*pi/2)*Coherentdk(N,0.0,0.0)

nCycles = 3
t = 3

χ0 = 0.3
θ0 = pi/2
ps0 = vcat(repeat([θ0],nCycles-1),χ0*rand(nCycles))
res = Optim.optimize(ps->acStates_SqueezingRotation(t,λy,Vy,Jz2,ψ,ps)[1],ps0,NelderMead(),Optim.Options(g_tol=1e-12,iterations=10^5))

res.minimum