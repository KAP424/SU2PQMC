# 2d Trotter Decomposition

function Initial_s(model::_Hubbard_Para,rng::MersenneTwister)::Array{UInt8,2}
    sp=Random.Sampler(rng,[1,2,3,4])

    s::Array{UInt8,2}=zeros(model.Ns,model.Nt)

    for i = 1:model.Ns
        for j = 1:model.Nt
            # 从elements中随机选择一个元素来填充当前位置  
            s[i, j] =rand(rng,sp)
        end  
    end  

    return s
end



function Gτ(model::_Hubbard_Para,s::Array{UInt8,2},τ::Int64)::Array{ComplexF64,2}
    """
    equal time Green function
    """
    BL::Array{ComplexF64,2}=model.Pt'[:,:]
    BR::Array{ComplexF64,2}=model.Pt[:,:]

    counter=0
    for i in model.Nt:-1:τ+1
        D=[model.η[x] for x in s[:,i]]
        BL=BL*diagm(exp.(1im*model.α.*D))*model.eK
        counter+=1
        if counter==model.BatchSize
            counter=0
            BL=Matrix(qr(BL').Q)'
        end
    end
    counter=0
    for i in 1:1:τ
        D=[model.η[x] for x in s[:,i]]
        BR=diagm(exp.(1im*model.α.*D))*model.eK*BR
        counter+=1
        if counter==model.BatchSize
            counter=0
            BR=Matrix(qr(BR).Q)
        end
    end

    BL=Matrix(qr(BL').Q)'
    BR=Matrix(qr(BR).Q)

    return I(model.Ns)-BR*inv(BL*BR)*BL
end


function G4(model::_Hubbard_Para,s::Array{UInt8,2},τ1::Int64,τ2::Int64)
    """
    displaced Green function
    return:
        G(τ₁),G(τ₂),G(τ₁,τ₂),G(τ₂,τ₁)
    """
    if τ1>τ2
        BBs=zeros(ComplexF64,cld(τ1-τ2,model.BatchSize),model.Ns,model.Ns)
        BBsInv=zeros(ComplexF64,size(BBs))
        
        UL=zeros(ComplexF64,1+size(BBs)[1],div(model.Ns,2),model.Ns)
        UR=zeros(ComplexF64,size(UL)[1],model.Ns,div(model.Ns,2))
        G=zeros(ComplexF64,size(UL)[1],model.Ns,model.Ns)

        UL[end,:,:]=model.Pt'[:,:]
        UR[1,:,:]=model.Pt[:,:]
    
        counter=0
        for i in 1:τ2
            D=[model.η[x] for x in s[:,i]]
            UR[1,:,:]=diagm(exp.(1im*model.α.*D))*model.eK*UR[1,:,:]
            counter+=1
            if counter==model.BatchSize
                counter=0
                UR[1,:,:]=Matrix(qr(UR[1,:,:]).Q)
            end
        end
        # UR[1,:,:]=UR[1,:,:]
        UR[1,:,:]=Matrix(qr(UR[1,:,:]).Q)
    
        counter=0
        for i in model.Nt:-1:τ1+1
            D=[model.η[x] for x in s[:,i]]
            UL[end,:,:]=UL[end,:,:]*diagm(exp.(1im*model.α.*D))*model.eK
            counter+=1
            if counter==model.BatchSize
                counter=0
                UL[end,:,:]=Matrix(qr(UL[end,:,:]').Q)'
            end
        end
        # UL[end,:,:]=UL[end,:,:]
        UL[end,:,:]=Matrix(qr(UL[end,:,:]').Q)'
    
        for i in 1:size(BBs)[1]-1
            BBs[i,:,:]=I(model.Ns)
            BBsInv[i,:,:]=I(model.Ns)
            for j in 1:model.BatchSize
                D=[model.η[x] for x in s[:,τ2+(i-1)*model.BatchSize+j]]
                BBs[i,:,:]=diagm(exp.(1im*model.α.*D))*model.eK*BBs[i,:,:]
                BBsInv[i,:,:]=BBsInv[i,:,:]*model.eKinv*diagm(exp.(-1im*model.α.*D))
            end
        end
    
        BBs[end,:,:]=I(model.Ns)
        BBsInv[end,:,:]=I(model.Ns)
        for j in τ2+(size(BBs)[1]-1)*model.BatchSize+1:τ1
            D=[model.η[x] for x in s[:,j]]
            BBs[end,:,:]=diagm(exp.(1im*model.α.*D))*model.eK*BBs[end,:,:]
            BBsInv[end,:,:]=BBsInv[end,:,:]*model.eKinv*diagm(exp.(-1im*model.α.*D))
        end
    
        for i in 1:size(BBs)[1]
            UL[end-i,:,:]=Matrix(qr( (UL[end-i+1,:,:]*BBs[end-i+1,:,:])' ).Q)'
            UR[i+1,:,:]=Matrix(qr(BBs[i,:,:]*UR[i,:,:]).Q)
        end
    
        for i in 1:size(G)[1]
            G[i,:,:]=I(model.Ns)-UR[i,:,:]*inv(UL[i,:,:]*UR[i,:,:])*UL[i,:,:]
            #####################################################################
            # if i <size(G)[1]
            #     if norm(Gτ(model,s,τ2+(i-1)*model.BatchSize)-G[i,:,:])>1e-3
            #         error("$i Gt")
            #     end
            # else
            #     if norm(Gτ(model,s,τ1)-G[i,:,:])>1e-3
            #         error("$i Gt")
            #     end
            # end
            #####################################################################
        end


        G12=I(model.Ns)
        G21=-I(model.Ns)
        for i in 1:size(BBs)[1]
            G12=G12*BBs[end-i+1,:,:]*G[end-i,:,:]
            G21=G21*( I(model.Ns)-G[i,:,:] )*BBsInv[i,:,:]
        end
        
        return G[end,:,:],G[1,:,:],G12,G21
    
    elseif τ1<τ2
        G2,G1,G21,G12=G4(model,s,τ2,τ1)
        return G1,G2,G12,G21
    else
        G=Gτ(model,s,τ1)
        return G,G,-(I(model.Ns)-G),G
    
    end

end


function GroverMatrix(G1,G2)
    II=I(size(G1)[1])
    return G1*G2+(II-G1)*(II-G2)
end

function Free_G(Lattice,site,Θ,Initial)
    """
    input:
        Lattice: "HoneyComb" or "SQUARE"
        site: [Int64,Int64]
        Θ: Float64
        Initial: "H0" or "V"
    return Green function with Free Hubbard Hamiltonian from H0 initial state or SDW initial state 
    对于得到H0初态(平衡态结果)
    必须要对H0加一个极其微弱的交错化学势,以去除基态简并,从而得到正确的结果
    """
    K=K_Matrix(Lattice,site)
    Ns=size(K)[1]
    

    Δt=0.1
    Nt=Int(round(Θ/Δt))

    if Initial=="H0"
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
    end

    E,V=eigen(K)
    eK=V*diagm(exp.(-Δt.*E))*V'

    if Initial=="H0"
        Pt=V[:,1:div(Ns,2)]
    elseif Initial=="V" 
        Pt=zeros(Float64,Ns,Int(Ns/2))
        for i in 1:Int(Ns/2)
            Pt[i*2,i]=1
        end
    end

    BL=Pt'[:,:]
    BR=Pt[:,:]

    count=0
    for i in 1:Nt
        BL=BL*eK
        BR=eK*BR
        count+=1
        if count==10
            counter=0
            BL=Matrix(qr(BL').Q)'
            BR=Matrix(qr(BR).Q)
        end
    end

    return I(Ns)-BR*inv(BL*BR)*BL
end

function G12FF(model,s,τ1,τ2)
    """
    Debug G(τ1,τ2) Green function
    Without numerical stablility
    Only for short time debug!!!
    """
    if τ1>τ2
        G=Gτ(model,s,τ2)
        BBs=I(model.Ns)
        BBsInv=I(model.Ns)

        for i in τ2+1:τ1
            D=[model.η[x] for x in s[:,i]]
            BBs=diagm(exp.(1im*model.α.*D))*model.eK*BBs
            BBsInv=BBsInv*model.eKinv*diagm(exp.(-1im*model.α.*D))
        end


        G12=BBs*G
        G21=-( I(model.Ns)-G ) * BBsInv

        return G12,G21
    elseif τ1<τ2
        G21,G12=G12FF(model,s,τ2,τ1)
        return G12,G21
    end
    
end

