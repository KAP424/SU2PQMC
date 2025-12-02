# using SU(2) ±1,±2 HS transformation
mutable struct UpdateBuffer_
    r::Matrix{ComplexF64}                  
    subidx::Vector{Int64}      
end

function UpdateBuffer()
    return UpdateBuffer_(
        Matrix{ComplexF64}(undef, 1, 1),
        [0],
    )    
end

struct Hubbard_Para_
    Lattice::String
    t::Float64
    U::Float64
    site::Vector{Int64}
    Θ::Float64
    Ns::Int64
    Nt::Int64
    K::Array{Float64,2}
    BatchSize::Int64
    WrapTime::Int64
    Δt::Float64
    α::Float64
    γ::Vector{Float64}
    η::Vector{Float64}
    Pt::Array{Float64,2}
    HalfeK::Array{Float64,2}
    eK::Array{Float64,2}
    HalfeKinv::Array{Float64,2}
    eKinv::Array{Float64,2}
    nodes::Vector{Int64}
    samplers_dict::Dict{UInt8, Random.Sampler}
end

function Hubbard_Para(t, U, Lattice::String, site, Δt, Θ, BatchSize, Initial::String)
    Nt = 2 * cld(Θ, Δt)
    WrapTime = div(BatchSize, 2)
    
    α = sqrt(Δt * U / 2)
    γ = [1 + sqrt(6) / 3, 1 + sqrt(6) / 3, 1 - sqrt(6) / 3, 1 - sqrt(6) / 3]
    η = [sqrt(2 * (3 - sqrt(6))), -sqrt(2 * (3 - sqrt(6))), sqrt(2 * (3 + sqrt(6))), -sqrt(2 * (3 + sqrt(6)))]


    K = K_Matrix(Lattice, site)
    Ns = size(K, 1)

    E, V = LAPACK.syevd!('V', 'L',K[:,:])
    HalfeK=V*Diagonal(exp.(-Δt.*E./2))*V'
    eK=V*Diagonal(exp.(-Δt.*E))*V'
    HalfeKinv=V*Diagonal(exp.(Δt.*E./2))*V'
    eKinv=V*Diagonal(exp.(Δt.*E))*V'

    Pt = zeros(Float64, Ns, div(Ns, 2))  # 预分配 Pt
    if Initial == "H0"
        KK = copy(K)
        μ = 1e-5
        if occursin("HoneyComb", Lattice)
            KK .+= μ * diagm(repeat([-1, 1], div(Ns, 2)))
        elseif Lattice == "SQUARE"
            for i in 1:Ns
                x, y = i_xy(Lattice, site, i)
                KK[i, i] += μ * (-1)^(x + y)
            end
        end
        E, V = LAPACK.syevd!('V', 'L',KK)
        Pt .= V[:, 1:div(Ns, 2)]
    elseif Initial=="V" 
        if occursin("HoneyComb", Lattice)
            for i in 1:div(Ns,2)
                Pt[i*2,i]=1
            end
        else
            count=1
            for i in 1:Ns
                x,y=i_xy(Lattice,site,i)
                if (x+y)%2==1
                    Pt[i,count]=1
                    count+=1
                end
            end
        end
    end
    
    if div(Nt, 2) % BatchSize == 0
        nodes = collect(0:BatchSize:Nt)
    else
        nodes = vcat(0, reverse(collect(div(Nt, 2) - BatchSize:-BatchSize:1)), collect(div(Nt, 2):BatchSize:Nt), Nt)
    end

    rng=MersenneTwister(Threads.threadid()+time_ns())
    elements = (1, 2, 3, 4)
    samplers_dict = Dict{UInt8, Random.Sampler}()
    for excluded in elements
        allowed = [i for i in elements if i != excluded]
        samplers_dict[excluded] = Random.Sampler(rng, allowed)
    end

    return Hubbard_Para_(Lattice, t, U, site, Θ, Ns, Nt, K, BatchSize, WrapTime, Δt, α, γ, η, Pt, HalfeK, eK, HalfeKinv, eKinv, nodes,samplers_dict)
end

# Buffers for phy_update workflow
mutable struct PhyBuffer_
	tau::Vector{ComplexF64}
	ipiv::Vector{LAPACK.BlasInt}

	G::Matrix{ComplexF64}
	BM::Matrix{ComplexF64}
	BLs::Array{ComplexF64,3}
	BRs::Array{ComplexF64,3}

	# generic temporaries
    N::Vector{ComplexF64}
	NN::Matrix{ComplexF64}
	Nn::Matrix{ComplexF64}
	nn::Matrix{ComplexF64}
	nN::Matrix{ComplexF64}
	zN::Matrix{ComplexF64}   # 1 x Ns
end

function PhyBuffer(Ns,NN)
    ns = div(Ns,2)
    return PhyBuffer_(
        Vector{ComplexF64}(undef, ns),
        Vector{LAPACK.BlasInt}(undef, ns),

        Matrix{ComplexF64}(undef, Ns, Ns),
        Matrix{ComplexF64}(undef, Ns, Ns),
        Array{ComplexF64}(undef, ns, Ns, NN),
        Array{ComplexF64}(undef, Ns, ns, NN),

        Vector{ComplexF64}(undef, Ns),
        Matrix{ComplexF64}(undef, Ns, Ns),
        Matrix{ComplexF64}(undef, Ns, ns),
        Matrix{ComplexF64}(undef, ns, ns),
        Matrix{ComplexF64}(undef, ns, Ns),
        Matrix{ComplexF64}(undef, 1, Ns),
    )
end
# ---------------------------------------------------------------------------------------
# Buffers for SCEE workflow
mutable struct SCEEBuffer_
    II::Matrix{ComplexF64}                 # Ns x Ns identity matrix (dense)
    N::Vector{ComplexF64}                 # Ns
    N_::Vector{ComplexF64}                # Ns
    zN::Matrix{ComplexF64}                # 1 x Ns
    nn::Matrix{ComplexF64}             
    NN::Matrix{ComplexF64}             
    NN_::Matrix{ComplexF64}            
    Nn::Matrix{ComplexF64}             
    nN::Matrix{ComplexF64}             
    ipiv::Vector{LAPACK.BlasInt}        # length ns
	tau::Vector{ComplexF64}                # length ns
end

mutable struct G4Buffer_
    Gt::Matrix{ComplexF64}
    G0::Matrix{ComplexF64}
    Gt0::Matrix{ComplexF64}
    G0t::Matrix{ComplexF64}
    BLMs::Array{ComplexF64,3}
    BRMs::Array{ComplexF64,3}
    BMs::Array{ComplexF64,3}
    BMinvs::Array{ComplexF64,3}
end

mutable struct AreaBuffer_
    index::Vector{Int64}          # length nA
    detg::ComplexF64
    gmInv::Matrix{ComplexF64}          # nA x nA
    NN::Matrix{ComplexF64}              # nA x nA
    Nz::Matrix{ComplexF64}              # nA x 2
    zN::Matrix{ComplexF64}              # 1 x nA
    a::Matrix{ComplexF64}               # nA x 1
    b::Matrix{ComplexF64}               # 1 x nA
    Tau::Matrix{ComplexF64}             # 1 x 1
	ipiv::Vector{LAPACK.BlasInt}        # length ns
end

function SCEEBuffer(Ns)
    ns=div(Ns, 2)
    return SCEEBuffer_(
        Matrix{ComplexF64}(I, Ns, Ns),
        Vector{ComplexF64}(undef, Ns),
        Vector{ComplexF64}(undef, Ns),
        Matrix{ComplexF64}(undef, 1, Ns),
        Matrix{ComplexF64}(undef, ns, ns),
        Matrix{ComplexF64}(undef, Ns, Ns),
        Matrix{ComplexF64}(undef, Ns, Ns),
        Matrix{ComplexF64}(undef, Ns, ns),
        Matrix{ComplexF64}(undef, ns, Ns),
        Vector{LAPACK.BlasInt}(undef, ns),
        Vector{ComplexF64}(undef, ns),
    )
end

function G4Buffer(Ns,NN)
    ns=div(Ns, 2)
    return G4Buffer_(
        Matrix{ComplexF64}(undef, Ns, Ns),
        Matrix{ComplexF64}(undef, Ns, Ns),
        Matrix{ComplexF64}(undef, Ns, Ns),
        Matrix{ComplexF64}(undef, Ns, Ns),

        Array{ComplexF64,3}(undef, ns, Ns, NN),
        Array{ComplexF64,3}(undef, Ns, ns, NN),
        Array{ComplexF64,3}(undef, Ns, Ns, NN),
        Array{ComplexF64,3}(undef, Ns, Ns, NN),
    )
end

function AreaBuffer(index)
    nA = length(index)
    return AreaBuffer_(
        index,
        0.0,
        Matrix{ComplexF64}(undef, nA, nA),
        Matrix{ComplexF64}(undef, nA, nA),
        Matrix{ComplexF64}(undef, nA, 1),
        Matrix{ComplexF64}(undef, 1, nA),
        Matrix{ComplexF64}(undef, nA, 1),
        Matrix{ComplexF64}(undef, 1, nA),
        Matrix{ComplexF64}(undef, 1, 1),
        Vector{LAPACK.BlasInt}(undef, nA),
    )
end


