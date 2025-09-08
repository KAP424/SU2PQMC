push!(LOAD_PATH,"E:/桌面/JuliaDQMC/code/SU2PQMC/")
using DelimitedFiles

using KAPDQMC
using LinearAlgebra
using Random
rng=MersenneTwister(1)

t=1;   Lattice="HoneyComb120"    
U=8;     Δt=0.1;     Θ=0.3;
BatchSize=10;
  

L=6
site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")
path="E:/桌面/JuliaDQMC/code/SU2PQMC/TEST/"

s=Initial_s(model,rng)

# # Half
indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
println(indexA)

# # HalfHalf
indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))
println(indexB)

λ=0.5
Nλ=2
Sweeps=1
ω=π/2
s=SC_DOP(path,model,ω,indexA,indexB,Sweeps,λ,Nλ,s,true)


path="kap/a"

path=path*"/"