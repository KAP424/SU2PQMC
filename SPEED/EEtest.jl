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
    # # HalfHalf
    indexB=area_index(Lattice,site,([1,1],[div(L,3),div(2*L,3)]))
    # println(indexB)
    ss=[Initial_s(model,rng)[:,:],Initial_s(model,rng)[:,:]]
    λ=0.5
    Nλ=2

    println(@btime ctrl_SCEEicr($path,$model,$indexA,$indexB,$Sweeps,$λ,$Nλ,$ss,$true) )
    # ctrl_SCEEicr(path,model,indexA,indexB,1,λ,Nλ,ss,true)
end

main()

# time = 163.9698069, bytes = 549133774, alloc = 2176314, gctime = 0.1244426)
# time = 153.8645969, bytes = 249095242, alloc = 2053734, gctime = 0.000782)
# time =  14.2534714, bytes = 168298175, alloc = 403, gctime = 0.0007381)
# 188.140 s (2053734 allocations: 237.56 MiB)

#------SCEEicr_PQMC-------#
# U=3.8;     Δt=0.1;     Θ=3.0;   L=15;    Sweeps=1;
# OLD 1025.714 seconds
# NEW 886.15 seconds
# !!! 115.434 s (6079507 allocations: 1.37 GiB)
# 128.842 s (6038300 allocations: 835.83 MiB)
# Final 95.108 s (5242872 allocations: 662.00 MiB)

# L=9
# no update 2.281 s (142399 allocations: 126.97 MiB)
# update 8.909 s (2188951 allocations: 240.55 MiB)
# update 8.779 s (2213106 allocations: 169.78 MiB)
# 8.625 s (2176897 allocations: 167.88 MiB)

using LinearAlgebra,LinearAlgebra.BLAS

a=rand(2,1)+rand(2,1)*im
b=rand(1,2)+rand(1,2)*im


dotu(a,b)

(b*a)[1,1]

dotc(a,b)
dotc(b,a)

