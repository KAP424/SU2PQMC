using BenchmarkTools, LinearAlgebra,Random
using LinearAlgebra.BLAS,LinearAlgebra.LAPACK
push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/")



n = 512
A = randn(n, n)
Ainv = similar(A)
ipiv = Vector{BLAS.BlasInt}(undef, n)
work = Vector{Float64}(undef, n * 64)

# workspace query
A_tmp = copy(A)
LAPACK.getrf!(A_tmp)
# lwork = LAPACK.getri_workspace_query('N', A_tmp)
# resize!(work, Int(real(lwork)))

# repeat inverse many times
for k = 1:1000
    copy!(Ainv, A)
    LAPACK.getrf!(Ainv, ipiv)
    LAPACK.getri!(Ainv, ipiv)
end
using  Base.Threads
@Threads.thread_local L=1