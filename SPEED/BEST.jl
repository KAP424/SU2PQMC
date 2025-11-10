using LinearAlgebra,LinearAlgebra.LAPACK, BenchmarkTools

B=rand(500,500)
C=rand(1000,500)
tau = Vector{Float64}(undef, 500)

@btime begin
    A=rand(1000,500)
    for i in 1:100
        mul!(C, A, B)
        copyto!(A, C)
        LAPACK.geqrf!(A, tau)
        LAPACK.orgqr!(A, tau)
    end
end
# (time = 3.3495654, bytes = 25628184, alloc = 900, gctime = 0.0)

@btimed begin
    A=rand(1000,500)
    for i in 1:100
        mul!(C, A, B)
        copyto!(A, C)
        A,tau=LAPACK.geqrf!(A)
        LAPACK.orgqr!(A,tau)
    end
end
# (time = 3.5229168, bytes = 30038155, alloc = 1103, gctime = 0.0)

@btimed begin
    A=rand(1000,500)
    for i in 1:100
        mul!(C, A, B)
        A = qr(C).Q
    end
end
# (time = 1.927754, bytes = 630833776, alloc = 1201, gctime = 0.0390994)


LAPACK.gerqf!
LAPACK.geqrt!

LAPACK.getri!


# 定义一个矩阵 A
A = rand(10,3)
Acopy=copy(A)


tau = Vector{Float64}(undef, size(A,2))
jpvt = Vector{LAPACK.BlasInt}(undef, size(A, 2))

LAPACK.geqp3!(A,jpvt,tau)
R = triu(A)[jpvt,:]
LAPACK.orgqr!(A, tau)
Acopy ≈ A*R

# 验证分解结果
println("列置换向量 jpvt: ", jpvt)
println("Q = ", Q)  # 正交矩阵 Q
println("R = ", R)  # 上三角矩阵 R
println("A[:, jpvt] ≈ Q * R: ", A[:, jpvt] ≈ Q * R)


@allocated A,tau=LAPACK.geqrf!(A)
@allocated LAPACK.geqrf!(A,tau)


@allocated LAPACK.orgqr!(A, tau)

@allocated qr(A).Q


A=qr(A).Q

A=A*C


Q,R,P=qr(A,Val(true))

A[:,P]≈Q * R

Matrix(qr(A).Q)

C=rand(500,500)

mul!(A, B, C)


A'*A

@allocated a.=A[:,1]


@btime C.= A.- Diagonal(II) 
@allocated @views C[1:100,1:100].=A[1:100,1:100].-B[1:100,1:100]

@allocated @views mul!(C, view(A,:,1), view(B,1:1,:))

axpy!(-Δ1/r1, tmpNN, G)

A = rand(1000,1000)+1im*rand(1000,1000)
A= A+ A'
AA = A[:,:]
ipiv = Vector{LAPACK.BlasInt}(undef, size(A, 1))
@allocated  LAPACK.getrf!(A,ipiv)
LAPACK.getri!(A, ipiv)
norm(AA*A-I(1000))


@allocated lmul!(2.0, A)

ns=div(model.Ns, 2)
NN=length(model.nodes)
tau = Vector{ComplexF64}(undef, ns)
LAPACK.geqrf!(tmpNn, tau)
LAPACK.orgqr!(tmpNn, tau, ns)


A=rand(100,100)

qr(A) 

QRPivoted()