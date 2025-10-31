using BenchmarkTools, LinearAlgebra,Random
using LinearAlgebra.BLAS,LinearAlgebra.LAPACK
push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/")

A=rand(10,100)*1im
A_copy=copy(A)
tau = Vector{ComplexF64}(undef, 10)
LAPACK.gerqf!(A,tau)
LAPACK.orgrq!(A,tau,10)
A*A'


tmpNn=A_copy'[:,:]
LAPACK.geqrf!(tmpNn, tau)
LAPACK.orgqr!(tmpNn, tau, 10)

norm(A-tmpNn')

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