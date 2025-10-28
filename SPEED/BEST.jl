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
