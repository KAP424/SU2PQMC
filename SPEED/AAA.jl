# ------------------------------------------------------------------------------------------------------------

using BenchmarkTools, LinearAlgebra,Random


# 预先为矩阵 A 分配内存
A = Matrix{Float64}(undef, 5, 3) # 创建一个 5x3 的未初始化矩阵
# 接下来可以初始化 A，例如用随机数填充
rand!(A) # 使用随机数覆盖 A 的内容

# 执行原地 QR 分解，F 的类型由编译器自动推断
F = qr!(A) # 现在 F 包含了分解结果，A 的内容已被覆盖为分解因子

# 查看 F 的具体类型
Matrix(F.Q)'
println("Type of F: ", typeof(F))

B=Matrix(qr((tmp*tmpNN)').Q)'

A*A'
B*B'


tmpNn'*tmpNn

C=rand(5,3)

R=qr!(C)

R.R

Matrix(R.Q)'*Matrix(R.Q)
C
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

