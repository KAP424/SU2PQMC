push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/")
using DelimitedFiles
using LinearAlgebra.LAPACK
using KAPDQMC
using LinearAlgebra
using Random

function GroverMatrix(G1,G2)
    println("Calculating Grover Matrix...")
    n = size(G1, 1)
    GM=Matrix{Float64}(undef,n,n)
    
    mul!(GM,G1,G2)
    lmul!(2.0, GM)
    axpy!(-1.0, G1, GM)
    axpy!(-1.0, G2, GM)
    for i in diagind(GM)
        GM[i] += 1.0
    end
    return GM   
    # 2*G1*G2 - G1 - G2 + II
end

using LinearAlgebra

function det_highprec(A::AbstractMatrix{T}; target_type::Type{<:AbstractFloat}=Float64, balance::Bool=false) where T
    # 1) 转为浮点并复制，避免破坏原矩阵
    B = float(copy(A))

    # 2) 可选：对称/厄米矩阵做平衡（改善条件数）
    if balance && (issymmetric(B) || ishermitian(B))
        B, scalefacs = balance(B)
        # 平衡后行列式需除以缩放因子的乘积
        s = sum(log.(abs.(scalefacs)))  # log-scale 避免溢出
        return sign(prod(scalefacs)) * exp(s) * det(B)
    end

    # 3) 目标类型为 BigFloat 时，提升到任意精度
    if target_type <: BigFloat
        B = big.(B)
    end

    # 4) 使用 LAPACK 路径的 det（LU 分解）
    return det(B)
end

# 使用示例：
# A = randn(1000,1000)
# d1 = det_highprec(A)                          # Float64，常规
# d2 = det_highprec(A, target_type=BigFloat)    # 任意精度
# d3 = det_highprec(A, balance=true)            # 对称/厄米病态时更稳

rng=MersenneTwister(1)

t=1;   Lattice="SQUARE"    
U=0;     Δt=0.1;     Θ=0.0;
BatchSize=5;
  

L=32
site=[L,L]

# # # Half
indexA=area_index(Lattice,site,([1,1],[div(L,2),L]))
# # println(indexA)

# # # HalfHalf
indexB=area_index(Lattice,site,([1,1],[div(L,2),div(L,2)]))
# # println(indexB)

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")

A=I(10)


println(sum(model.Pt,dims=1))
# μ = 1e-3
K=model.K[:,:]
# for i in 1:model.Ns
#     x, y = i_xy(Lattice, site, i)
#     K[i, i] += μ * (-1)^(x + y)
# end

println(norm(K-K'))

E,V=LAPACK.syevd!('V', 'L',K[:,:])
E1,V1=eigen(Symmetric(K))
println("eigen is wrong!")
println(norm(E-E1)," ",norm(V-V1))

println("ortho check:")
println(norm(V*V'-I(size(K,1))),"   ",norm(V*Diagonal(E)*V'-K))

Δt=0.1
eK = V * Diagonal(exp.(-Δt.*E)) * V'
println("Diag check:",norm(Diagonal(exp.(-Δt.*E))-I(model.Ns)))
println("eK check:",norm(eK-I(model.Ns)))

println("max E=",maximum(E)," min E=",minimum(E))

BL=model.Pt'*eK
BR=eK*model.Pt
# BL=model.Pt'[:,:]
# BR=model.Pt[:,:]

tmp=BL*BR
ipiv = Vector{LAPACK.BlasInt}(undef, size(tmp, 1))
@allocated  LAPACK.getrf!(tmp,ipiv)
LAPACK.getri!(tmp, ipiv)
G=I(model.Ns)-BR*tmp*BL

GMA=GroverMatrix(G[indexA,indexA],G[indexA,indexA])
GM=2*G[indexA,indexA]*G[indexA,indexA]-2*G[indexA,indexA]+I(length(indexA))

println(norm(GMA-GM))


# println(det(GMA)-det_highprec(GMA, balance=true) ,"   :",det_highprec(GMA))

EE= - log(abs2(det(GMA)))

println("EE= ",EE)




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
# # TEST for SCEE
# path="C:/Users/admin/Desktop/"



# ss=[s[:,:],s[:,:]]
# λ=0.5
# Nλ=2
# Sweeps=1
# # ctrl_SCEEicr(path,model,indexA,indexB,Sweeps,λ,Nλ,ss,false)

# # ss=ctrl_EEicr(path,model,indexA,10,0.0,1,ss,true)

# phy_update(path,model,2,2,s)
