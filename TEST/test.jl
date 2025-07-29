push!(LOAD_PATH,"E:/桌面/JuliaDQMC/code/SU2PQMC/")
using DelimitedFiles

using KAPDQMC
using LinearAlgebra
using Random
rng=MersenneTwister(1)

t=1;   Lattice="HoneyComb"    
U=8;     Δt=0.1;     Θ=0.3;
BatchSize=10;
  

L=6
site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")

# ------------------------------------------------------------------------
# TEST for nn2idx
# for i in 1:18
#     print(i)
#     println(nn2idx(Lattice,site,i))
# end

# K=K_Matrix(Lattice,site)
# print(K)
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# TEST for Green function
s=Initial_s(model,rng)

τ=model.Nt
# τ=1
G=Gτ(model,s,div(model.Nt,2))
Gt,G0,Gt0,G0t=G4(model,s,τ,div(model.Nt,2))
Gt0_,G0t_=G12FF(model,s,τ,div(model.Nt,2))
println(norm(Gt0-Gt0_),',',norm(G0t-G0t_))
Gt_=Gτ(model,s,τ)
G0_=Gτ(model,s,div(model.Nt,2))
println(norm(Gt-Gt_),',',norm(G0-G0_))
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# TEST for SCEE
path="E:/桌面/JuliaDQMC/code/SU2PQMC/TEST/"

# # Half
indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
println(indexA)

# # HalfHalf
indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))
println(indexB)

ss=[s[:,:],s[:,:]]
λ=0.5
Nλ=2
Sweeps=1
ctrl_SCEEicr(path,model,indexA,indexB,Sweeps,λ,Nλ,ss,true)


