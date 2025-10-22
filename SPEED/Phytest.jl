push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/")
using DelimitedFiles
using BenchmarkTools
using MKL
using KAPDQMC
using LinearAlgebra
using Random 
rng=MersenneTwister(1)

t=1;   Lattice="HoneyComb120"    
U=3.8;     Δt=0.1;     Θ=3.0;
BatchSize=10;
Sweeps=1
L=15
# OLD 89.96 seconds
# NEW 81.889 seconds
# now 58.411 seconds
# 50.744 s (294197 allocations: 157.09 GiB)
# 51.852 s (237912 allocations: 105.04 GiB)
# 44.309 s (640499 allocations: 105.05 GiB)
# two with update 74.987 s (831433 allocations: 179.65 GiB)
# two with no update 15.718 s (2335 allocations: 1.26 GiB)
# two with no update 15.277 s (2041 allocations: 1.13 GiB)
# 77.805 s (1237996 allocations: 192.05 GiB)
#   simple *    160.329 s (1332784 allocations: 287.40 GiB)
# 110.278 s (1237996 allocations: 192.05 GiB)
# 73.387 s (1213511 allocations: 179.14 GiB)

site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")

s=Initial_s(model,rng)
path="C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/SPEED/"

println(@btime phy_update($path,$model,1,0,$Initial_s(model,rng)))
# @benchmark phy_update(path,model,1,0,Initial_s(model,rng))



# function BM_F22(model::_Hubbard_Para, s::Array{UInt8, 2}, idx::Int64)::Array{ComplexF64, 2}
#     Ns=model.Ns
#     nodes=model.nodes
#     eK=model.eK
#     η=model.η
#     α=model.α

#     if idx > length(nodes) || idx < 1
#         error("idx out of range")
#     end
#     D = Array{ComplexF64,1}(undef, Ns)  # 预分配 D 数组
#     tmp=Matrix{ComplexF64}(undef, Ns, Ns)
#     BM=Matrix{ComplexF64}(undef,Ns,Ns)
#     @inbounds for i in diagind(BM)
#         BM[i] = one(eltype(BM))
#     end
#     for lt in nodes[idx] + 1:nodes[idx + 1]
#         @inbounds begin
#             for i in 1:Ns
#                 D[i] =  cis( α *η[s[i, lt]])
#             end
#             mul!(tmp,eK, BM)
#             # (BM, tmp) = (tmp, BM)
#             # for i ∈ 1:Ns
#             #     BM[i, :] .*= D[i]
#             # end
#             for i in 1:Ns
#                 BM[i, :] .= D[i] * tmp[i, :]
#             end
#             # mul!(BM,Diagonal(D), tmp)
#             # BM = diagm(expD) * eK * BM
#         end
#     end
#     return BM
# end

# function BM_F2(model::_Hubbard_Para, s::Array{UInt8, 2}, idx::Int64)::Array{ComplexF64, 2}
#     Ns=model.Ns
#     nodes=model.nodes
#     eK=model.eK
#     η=model.η
#     α=model.α

#     if idx > length(nodes) || idx < 1
#         error("idx out of range")
#     end
#     D = Array{ComplexF64,1}(undef, Ns)  # 预分配 D 数组
#     # tmp=Matrix{ComplexF64}(undef, Ns, Ns)
#     BM=Matrix{ComplexF64}(undef,Ns,Ns)
#     @inbounds for i in diagind(BM)
#         BM[i] = one(eltype(BM))
#     end
#     for lt in nodes[idx] + 1:nodes[idx + 1]
#         @inbounds begin
#             for i in 1:Ns
#                 D[i] =  cis( α *η[s[i, lt]])
#             end
#             # mul!(tmp,eK, BM)
#             # mul!(BM,diagm(D), tmp)
#             BM = Diagonal(D) * eK * BM
#         end
#     end

#     return BM
# end

# function BM_F3(model::_Hubbard_Para, s::Array{UInt8, 2}, idx::Int64)::Array{ComplexF64, 2}
#     Ns=model.Ns
#     nodes=model.nodes
#     eK=model.eK
#     η=model.η
#     α=model.α

#     if idx > length(nodes) || idx < 1
#         error("idx out of range")
#     end
#     D = Array{ComplexF64,1}(undef, Ns)  # 预分配 D 数组
#     # tmp=Matrix{ComplexF64}(undef, Ns, Ns)
#     BM=Matrix{ComplexF64}(undef,Ns,Ns)
#     @inbounds for i in diagind(BM)
#         BM[i] = one(eltype(BM))
#     end
#     for lt in nodes[idx] + 1:nodes[idx + 1]
#         @inbounds begin
#             for i in 1:Ns
#                 @fastmath D[i] =  cis( α *η[s[i, lt]])
#             end
#             # mul!(tmp,eK, BM)
#             # mul!(BM,diagm(D), tmp)
#             @fastmath BM = Diagonal(D) * eK * BM
#         end
#     end

#     return BM
# end

# @benchmark BM_F2($model,$s,1)
# @benchmark BM_F3($model,$s,1)

# Threads.nthreads()
# norm(BM_F22(model,s,1)-BM_F2(model,s,1))

# ------------------------------------------------------------------------
# TEST for SCEE

# phy_update(path,model,0,Sweeps,Initial_s(model,rng))

# # # Half
# indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
# 2.314 s (314 allocations: 793.47 MiB)


# # # HalfHalf
# indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))

# ss=[Initial_s(model,rng)[:,:],Initial_s(model,rng)[:,:]]
# λ=0.5
# Nλ=2

# ctrl_SCEEicr(path,model,indexA,indexB,Sweeps,λ,Nλ,ss,true)



# using LinearAlgebra, BenchmarkTools

# n = 1000
# A = rand(n, n)
# B = rand(n, n)
# C = similar(A)

# # 测试 1: 直接赋值
# @btime $A = $B * $C;

# # 测试 2: 原地操作
# @btime mul!($A, $B, $C);



