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
BatchSize=1000;
Sweeps=1
L=15

#------------------phy_PQMC-------------------#
# U=3.8;     Δt=0.1;     Θ=3.0;   L=15;       Sweeps=1;
# OLD 89.96 seconds
# NEW 81.889 seconds
# two with no update 15.718 s (2335 allocations: 1.26 GiB)
# two with no update 15.277 s (2041 allocations: 1.13 GiB)
# 77.805 s (1237996 allocations: 192.05 GiB)
#   simple *    160.329 s (1332784 allocations: 287.40 GiB)
# 110.278 s (1237996 allocations: 192.05 GiB)
# 73.387 s (1213511 allocations: 179.14 GiB)
# 扫描不更新 9.385 s (433099 allocations: 550.44 MiB)
# 扫描更新 44.226 s (640088 allocations: 104.76 GiB)
# 11.916 s (519345 allocations: 672.02 MiB)
# 扫描更新！！！   19.314 s (1005949 allocations: 954.42 MiB)
# 扫描不更新！！！ 14.464 s (1566 allocations: 729.15 MiB)
# 18.346 s (1005819 allocations: 898.78 MiB)
# 17.965 s (1162997 allocations: 326.09 MiB)

# Final 17.931 s (1155434 allocations: 256.35 MiB)

#------SCEEicr_PQMC-------#
# U=3.8;     Δt=0.1;     Θ=3.0;   L=15;    Sweeps=1;
# OLD 1025.714 seconds
# NEW 886.15 seconds
# !!! 115.434 s (6079507 allocations: 1.37 GiB)
# 128.842 s (6038300 allocations: 835.83 MiB)
# Final 95.108 s (5242872 allocations: 662.00 MiB)

site=[L,L]

model=Hubbard_Para(t,U,Lattice,site,Δt,Θ,BatchSize,"V")

s=Initial_s(model,rng)
path="C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/SPEED/"
println(@btime phy_update($path,$model,1,0,$Initial_s(model,rng)))
# phy_update(path,model,1,0,Initial_s(model,rng))
