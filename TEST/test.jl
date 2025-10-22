push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/")
using DelimitedFiles

using KAPDQMC
using LinearAlgebra
using Random
rng=MersenneTwister(1)

t=1;   Lattice="HoneyComb120"    
U=8;     Δt=0.1;     Θ=3.8;
BatchSize=5;
  

L=6
site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")

# println(model.nodes)

s=Initial_s(model,rng)


# ------------------------------------------------------------------------
# BMs=zeros(ComplexF64,length(model.nodes)-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns
# BMinvs=zeros(ComplexF64,length(model.nodes)-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns
# println("--------Test BMs and BMinvs--------")

# for idx in axes(BMs,1)
#     BMs[idx,:,:]=BM_F(model,s,idx)
#     BMinvs[idx,:,:]=BMinv_F(model,s,idx)
#     println(norm(BMs[idx,:,:]*BMinvs[idx,:,:]-I))
# end

# BLMs=zeros(ComplexF64,length(model.nodes),div(model.Ns,2),model.Ns)
# BRMs=zeros(ComplexF64,length(model.nodes),model.Ns,div(model.Ns,2))
# BLMs[end,:,:]=model.Pt'[:,:]
# BRMs[1,:,:]=model.Pt[:,:]
# for i in axes(BMs,1)
#     BLMs[end-i,:,:]=Matrix(qr( (BLMs[end-i+1,:,:]*BMs[end-i+1,:,:])' ).Q)'
#     BRMs[i+1,:,:]=Matrix(qr( BMs[i,:,:]*BRMs[i,:,:] ).Q)
# end

# println("--------Test BLMs and BRMs--------")
# for i in axes(BLMs,1)
#     G1=I(model.Ns)-BRMs[i,:,:] * inv( BLMs[i,:,:] * BRMs[i,:,:] ) * BLMs[i,:,:]
#     G2=Gτ_old(model,s,model.nodes[i])
#     println(norm(G1-G2))
# end

# println("--------Test G4--------")
# for i in eachindex(model.nodes)
#     lt=model.nodes[i]
#     println("lt=",lt)
#     Gt,G0,Gt0,G0t=G4_old(model,s,lt,div(model.Nt,2))
#     Gt_,G0_,Gt0_,G0t_=G4(model.nodes,lt,BLMs,BRMs,BMs,BMinvs)

#     println(norm(Gt-Gt_),',',norm(G0-G0_),',',norm(Gt0-Gt0_),',',norm(G0t-G0t_))

# end

# ------------------------------------------------------------------------


# for lt in 1:model.Nt
#     println(norm(Gτ(model,s,lt)-Gτ_old(model,s,lt)))
# end
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


# ------------------------------------------------------------------------
# TEST for SCEE
path="C:/Users/admin/Desktop/"

# # Half
indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
# println(indexA)

# # HalfHalf
indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))
# println(indexB)

ss=[s[:,:],s[:,:]]
λ=0.5
Nλ=2
Sweeps=1
# ctrl_SCEEicr(path,model,indexA,indexB,Sweeps,λ,Nλ,ss,false)

# ss=ctrl_EEicr(path,model,indexA,10,0.0,1,ss,true)

phy_update(path,model,2,2,s)