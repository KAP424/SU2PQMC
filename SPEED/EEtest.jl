push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/")
using DelimitedFiles
using BenchmarkTools
using KAPDQMC
using LinearAlgebra
using Random

function main()
    rng=MersenneTwister(1)

    t=1;   Lattice="HoneyComb120"    
    U=3.8;     Δt=0.1;     Θ=3.0;
    BatchSize=10;
    Sweeps=1
    L=15
    site=[L,L]

    model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")
    println(model.nodes)
    path="C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/SPEED/"

    # # Half
    indexA=area_index(Lattice,site,([1,1],[div(L,3),L]))
    # println(indexA)
    # # HalfHalf
    indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))
    # println(indexB)
    ss=[Initial_s(model,rng)[:,:],Initial_s(model,rng)[:,:]]
    λ=0.5
    Nλ=2

    println(@btime ctrl_SCEEicr($path,$model,$indexA,$indexB,$Sweeps,$λ,$Nλ,$ss,$true) )
end

main()

#------SCEEicr_PQMC-------#
# OLD 1025.714 seconds
# NEW 886.15 seconds
# !!! 115.434 s (6079507 allocations: 1.37 GiB)





# s=Initial_s(model,rng)
# G1=Gτ_old(model,s,1)
# G2=Gτ_old(model,s,2)

# gmInv_A=GroverMatrix(view(G1,indexA,indexA),view(G2,indexA,indexA))
# GM=gmInv_A[:,:]
# ipivA = Vector{LAPACK.BlasInt}(undef, length(indexA))
# LAPACK.getrf!(gmInv_A,ipivA)
# LAPACK.getri!(gmInv_A, ipivA)

# println(norm(GM*gmInv_A-I(length(indexA))))


# ------------------------------------------------------------------------
# TEST for SCEE



