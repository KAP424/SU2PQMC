# ------------------------------------------------------------------------------------------------------------
using BenchmarkTools, LinearAlgebra,Random
using LinearAlgebra.BLAS
push!(LOAD_PATH,"C:/Users/admin/Desktop/JuliaDQMC/code/SU2PQMC/")

A=rand(1000,1000)
B=rand(1000,1000)
C=rand(1000,1000)

@btime A.=B'
# 605.400 μs (2 allocations: 32 bytes)

function QRRR(A)
    tau=Vector{Float64}(undef, size(A,2))
    LAPACK.geqrf!(A, tau)
    LAPACK.orgqr!(A, tau,  size(A,2))
end
function QRRR1(A)
    Matrix(qr!(A).Q)
end

A=rand(5000,3000)
@btime QRRR($A)
@btime QRRR1($A)

#   1.421 s (11 allocations: 1.49 MiB)
# 1.920 ms (11 allocations: 915.93 KiB)


@btime A .= B-C

BLAS.geam!('N', 'N', 1.0, A, -1.0, B, C)

@btime geam!('N', 'N', 1.0, II, -1.0, tmpNN, G)


A=rand(5,3)
LAPACK.geqrf!(A,tau)
LAPACK.orgqr!(A, tau, 3)

norm(A'*A - I(3))

# 步骤2: 使用 orgqr! 从 geqrf! 的结果中生成完整的Q矩阵
    # orgqr! 同样会原地修改传入的矩阵，并返回生成的Q矩阵
    # 我们传入 A_copy 和 tau，函数会根据分解数据生成Q
A = LinearAlgebra.LAPACK.orgqr!(A, tau)
    
tau

include("../GlobalVars.jl")

using .GlobalVars


a

function Initial()
    global a=1
    global b=[1,2,3]
    global c=rand(3,3)
end

function main()
    Initial()
    
    print(a,"\n",b,"\n",c,"\n")

    a=a+1
    print(a,"\n",b,"\n",c,"\n")
end

main()


# function aa1()
#     A=Matrix{Float64}(undef, 10000,10000)
#     rand!(A)
#     A=5*A
# end

# function aa2()
#     A=Matrix{Float64}(undef, 10000,10000)
#     rand!(A)
#     A.=5*A
# end

# @benchmark aa1()
# @benchmark aa2()

# A=Matrix{Float64}(undef, 10,5)  
# B=rand(10,5)  

# A.= qr(B).Q


# qr(B).Q

# function d1()
#     A = Matrix{Float64}(undef, 10000,10000)
#     rand!(A)
#     D=rand(10000)
#     A*=diagm(D)
# end
# function d2()
#     A = Matrix{Float64}(undef, 10000,10000)
#     rand!(A)
#     D=rand(10000)
#     A*=Diagonal(D)
# end
# @benchmark d1()
# @benchmark d2()






# function Initial_0() 
#     A = Array{Float64}(undef, 10000,10000)
#     B = reshape(collect(1:10000*10000),10000,10000)
#     A.=B
# end

# function Initial_1()
#     A = Array{Float64}(undef, 10000,10000)
#     B = reshape(collect(1:10000*10000),10000,10000)
#     A.=B[:,:]
# end

# function Initial_2() 
#     A = Array{Float64}(undef, 10000,10000)
#     B = reshape(collect(1:10000*10000),10000,10000)
#     A=B[:,:]
# end

# function Initial_3() 
#     A = Array{Float64}(undef, 10000,10000)
#     B = reshape(collect(1:10000*10000),10000,10000)
#     @views A.=B
# end

# @benchmark Initial_0()
# @benchmark Initial_1()
# @benchmark Initial_2()
# @benchmark Initial_3()



# function t1()
#     D = Array{Float64,1}(undef, 10000000)  # 预分配 D 数组
#     expD =  Array{ComplexF64,1}(undef, 10000000) 
#     rand!(D)
#     expD .= exp.(1im * 1.23 * D)
# end

# function t2()
#     D = Array{ComplexF64,1}(undef, 10000000)  # 预分配 D 数组
#     rand!(D)
#     D .= exp.(1im * 1.23 * D)
# end

# function t3()   # fastest
#     D = Array{ComplexF64,1}(undef, 10000000)  # 预分配 D 数组
#     rand!(D)
#     @. D= exp(1im * 1.23 * D)
# end

# @benchmark t1()
# @benchmark t2()
# @benchmark t3()

# using BenchmarkTools, LinearAlgebra ,Random
# function MMP1()
#     a = Matrix{Float64}(undef, 1000, 1000)
#     b = Matrix{Float64}(undef, 1000, 1000)
#     c = Matrix{Float64}(undef, 1000, 1000)
#     @inbounds for i in diagind(c)
#         c[i] = one(eltype(c))
#     end
#     for i in 1:100
#         rand!(a)
#         rand!(b)
#         c=a*b*c
#         n = norm(c)
#         @inbounds c .= c / n
#     end
#     return c
# end

# function MMP2()
#     a = Matrix{Float64}(undef, 1000, 1000)
#     b = Matrix{Float64}(undef, 1000, 1000)
#     tmp = Matrix{Float64}(undef, 1000, 1000)
#     c = Matrix{Float64}(undef, 1000, 1000)
#     @inbounds for i in diagind(c)
#         c[i] = one(eltype(c))
#     end
#     for i in 1:100
#         rand!(a)
#         rand!(b)
#         mul!(tmp,b,c)
#         mul!(c,a,tmp)
#         n = norm(c)
#         @inbounds c .= c / n
#     end
#     return c
# end

# function MMP3() # fastest and easiest
#     a = Matrix{Float64}(undef, 1000, 1000)
#     b = Matrix{Float64}(undef, 1000, 1000)
#     tmp = Matrix{Float64}(undef, 1000, 1000)
#     c = Matrix{Float64}(undef, 1000, 1000)
#     @inbounds for i in diagind(c)
#         c[i] = one(eltype(c))
#     end
#     for i in 1:100
#         rand!(a)
#         rand!(b)
#         mul!(tmp,b,c)
#         mul!(c,a,tmp)
#         c .= c / norm(c)
#     end
#     return c
# end

# @benchmark MMP1() #   2.952 s (911 allocations: 2.26 GiB)
# @benchmark MMP3() 
# @benchmark MMP2() 

# ------------------------------------------------------------------

# using BenchmarkTools, LinearAlgebra,Random


# function mmm1()
#     A=Matrix{Float64}(undef,10000,1)
#     B=Matrix{Float64}(undef,1,10000)
#     C=Matrix{Float64}(undef,10000,10000)
#     rand!(A)
#     rand!(B)
#     C = A * B
# end

# function mmm2()
#     A=Matrix{Float64}(undef,10000,1)
#     B=Matrix{Float64}(undef,1,10000)
#     C=Matrix{Float64}(undef,10000,10000)
#     rand!(A)
#     rand!(B)
#     mul!(C,A,B)
# end

# function mmm3()
#     A=Matrix{Float64}(undef,10000,1)
#     B=Matrix{Float64}(undef,1,10000)
#     C=Matrix{Float64}(undef,10000,10000)
#     rand!(A)
#     rand!(B)
#     @inbounds begin
#         for i in axes(A,1)
#             for j in axes(B,2)
#                 C[i,j] = A[i,1] * B[1,j]
#             end
#         end
#     end 
# end

# @benchmark mmm1()
# @benchmark mmm2()
# @benchmark mmm3()


# function add1()
#     A=rand(10000,10000)
#     B=rand(10000,10000)
#     C=similar(A)
#     C=A+B
# end

# function add2()
#     A=rand(10000,10000)
#     B=rand(10000,10000)
#     C=similar(A)
#     C.=A+B
# end
# function add3()
#     A=rand(10000,10000)
#     B=rand(10000,10000)
#     C=similar(A)
#     @fastmath C.=A+B
# end
# @benchmark add1()
# @benchmark add2()
# @benchmark add3()

