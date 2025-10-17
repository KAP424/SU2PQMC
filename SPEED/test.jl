push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/")
using DelimitedFiles

using KAPDQMC
using LinearAlgebra
using Random
rng=MersenneTwister(1)

t=1;   Lattice="HoneyComb120"    
U=3.8;     Δt=0.1;     Θ=3.0;
BatchSize=10;
Sweeps=1
L=15
# 1025.714 seconds

site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")

# ------------------------------------------------------------------------
# TEST for SCEE
path="C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/SPEED/"

# # Half
indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))

# # HalfHalf
indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))

ss=[Initial_s(model,rng)[:,:],Initial_s(model,rng)[:,:]]
λ=0.5
Nλ=2

ctrl_SCEEicr(path,model,indexA,indexB,Sweeps,λ,Nλ,ss,true)


