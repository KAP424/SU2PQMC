# using SU(2) ±1,±2 HS transformation

mutable struct _Hubbard_Para
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
end

function Hubbard_Para(t,U,Lattice::String,site,Δt,Θ,BatchSize,Initial::String)
    Nt::Int64=2*cld(Θ,Δt)
    WrapTime::Int64=div(BatchSize,2)
    
    α::Float64=sqrt(Δt*U/2)
    γ::Vector{Float64}=[1+sqrt(6)/3,1+sqrt(6)/3,1-sqrt(6)/3,1-sqrt(6)/3]
    η::Vector{Float64}=[sqrt(2*(3-sqrt(6))),-sqrt(2*(3-sqrt(6))),sqrt(2*(3+sqrt(6))),-sqrt(2*(3+sqrt(6)))]
    
    K::Array{Float64,2}=K_Matrix(Lattice,site)
    Ns::Int64=size(K)[1]

    E,V=eigen(t*K)
    
    HalfeK=V*diagm(exp.(-Δt.*E./2))*V'
    eK=V*diagm(exp.(-Δt.*E))*V'
    HalfeKinv=V*diagm(exp.(Δt.*E./2))*V'
    eKinv=V*diagm(exp.(Δt.*E))*V'

    if Initial=="H0"
        # 交错化学势，打开gap，去兼并
        KK=K[:,:]
        μ=1e-5
        if occursin("HoneyComb", Lattice)
            KK+=μ*diagm(repeat([-1, 1], div(Ns, 2)))
        elseif Lattice=="SQUARE"
            for i in 1:Ns
                x,y=i_xy(Lattice,site,i)
                KK[i,i]+=μ*(-1)^(x+y)
            end
        end
        E,V=eigen(KK)
        Pt=V[:,1:div(Ns,2)]
        # K[K .!= 0] .+=( rand(size(K)...) * 0.1)[K.!= 0]
        # K=(K+K')./2
    elseif Initial=="V" 
        Pt=zeros(Float64,Ns,Int(Ns/2))
        for i in 1:Int(Ns/2)
            Pt[i*2,i]=1
        end
    end

    return _Hubbard_Para(Lattice,t,U,site,Θ,Ns,Nt,K,BatchSize,WrapTime,Δt,α,γ,η,Pt,HalfeK,eK,HalfeKinv,eKinv)

end


function setμ(model::_Hubbard_Para,μ)
    # fig1:1d
    # km=abs(acos(μ/2))
    # N_particle=Int(round(km/π*model.Ns))
    
    # fig2:2d-circle Fermi surface
    N_particle=Int(round( μ^2/4/π *model.Ns ))
    
    E,V=eigen(model.K)
    model.Pt=V[:,1:N_particle]
end

# function setNN_hopping(t)
    
# end