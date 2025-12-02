# attractive-U and repulsive-U get the same S_2


function ctrl_SCEEicr(path::String,model::Hubbard_Para_,indexA::Vector{Int64},indexB::Vector{Int64},Sweeps::Int64,λ::Float64,Nλ::Int64,ss::Vector{Matrix{UInt8}},record)
    global LOCK=ReentrantLock()
    Ns=model.Ns
    ns=div(Ns, 2)
    NN=length(model.nodes)
    Θidx=div(NN,2)+1

    UPD = UpdateBuffer()
    SCEE=SCEEBuffer(model.Ns)
    A=AreaBuffer(indexA)
    B=AreaBuffer(indexB)
    G1=G4Buffer(model.Ns,NN)
    G2=G4Buffer(model.Ns,NN)

    name = if model.Lattice=="SQUARE" "□" 
        elseif model.Lattice=="HoneyComb60" "HC" 
        elseif model.Lattice=="HoneyComb120" "HC120" 
        else error("Lattice: $(model.Lattice) is not allowed !") end    
    file="$(path)/SCEEicr$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)N$(Nλ)BS$(model.BatchSize).csv"
    rng=MersenneTwister(Threads.threadid()+time_ns())
    
    atexit() do 
        if record
            lock(LOCK) do
                open(file, "a") do io
                    writedlm(io, O', ',')
                end
            end
        end
        # writedlm("$(path)ss/SS$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)λ$(Int(round(Nλ*λ))).csv", [ss[1] ss[2]],",")
    end
    
    Gt1,Gt01,G0t1,BLMs1,BRMs1,BMs1,BMsinv1 =
        G1.Gt, G1.Gt0, G1.G0t, G1.BLMs, G1.BRMs, G1.BMs, G1.BMinvs
    Gt2,Gt02,G0t2,BLMs2,BRMs2,BMs2,BMsinv2 =
        G2.Gt, G2.Gt0, G2.G0t, G2.BLMs, G2.BRMs, G2.BMs, G2.BMinvs

    tmpN ,tmpN_,tmpNN,tmpNn,tmpnN,tau = SCEE.N,SCEE.N_,SCEE.NN,SCEE.Nn,SCEE.nN,SCEE.tau

    
    tmpO=0.0
    counter=0
    O=zeros(Float64,Sweeps+1)
    O[1]=λ

    for idx in 1:NN-1
        BM_F!(tmpN,tmpNN,view(BMs1,:, : , idx),model,ss[1],idx)
        BM_F!(tmpN,tmpNN,view(BMs2,:,:,idx),model,ss[2],idx)
        BMinv_F!(tmpN,tmpNN,view(BMsinv1,:,:,idx),model,ss[1],idx)
        BMinv_F!(tmpN,tmpNN,view(BMsinv2,:,:,idx),model,ss[2],idx)
    end

    transpose!(view(BLMs1,:,:,NN) , model.Pt)
    copyto!(view(BRMs1,:,:,1) , model.Pt)
    
    transpose!(view(BLMs2,:,:,NN) , model.Pt)
    copyto!(view(BRMs2,:,:,1) , model.Pt)

    for i in 1:NN-1
        mul!(tmpnN,view(BLMs1,:,:,NN-i+1),view(BMs1,:,:,NN-i))
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau, ns)
        copyto!(view(BLMs1,:,:,NN-i) , tmpnN)
        
        mul!(tmpNn, view(BMs1,:,:,i), view(BRMs1,:,:,i))
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau, ns)
        copyto!(view(BRMs1,:,:,i+1) , tmpNn)
        # ---------------------------------------------------------------
        mul!(tmpnN,view(BLMs2,:,:,NN-i+1),view(BMs2,:,:,NN-i))
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau, ns)
        copyto!(view(BLMs2,:,:,NN-i) , tmpnN)

        mul!(tmpNn, view(BMs2,:,:,i), view(BRMs2,:,:,i))
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau, ns)
        copyto!(view(BRMs2,:,:,i+1) , tmpNn)

    end

    idx=1
    get_ABGM!(G1,G2,A,B,SCEE,model.nodes,idx,"Forward")
    for loop in 1:Sweeps
        # println("\n ====== Sweep $loop / $Sweeps ======")
        for lt in 1:model.Nt
            @inbounds @simd for iii in 1:Ns
                @fastmath tmpN[iii] = cis( model.α *model.η[ss[1][iii, lt]] ) 
                @fastmath tmpN_[iii] = cis( model.α *model.η[ss[2][iii, lt]] ) 
            end

            WrapKV!(tmpNN,model.eK,model.eKinv,tmpN,Gt01,"Forward", "L")
            WrapKV!(tmpNN,model.eK,model.eKinv,tmpN_,Gt02,"Forward", "L")
            WrapKV!(tmpNN,model.eK,model.eKinv,tmpN,Gt1,"Forward", "B")
            WrapKV!(tmpNN,model.eK,model.eKinv,tmpN_,Gt2,"Forward", "B")
            WrapKV!(tmpNN,model.eK,model.eKinv,tmpN,G0t1,"Forward", "R")
            WrapKV!(tmpNN,model.eK,model.eKinv,tmpN_,G0t2,"Forward", "R")

            #####################################################################
            Gt1_,G01_,Gt01_,G0t1_=G4_old(model,ss[1],lt,div(model.Nt,2))
            Gt2_,G02_,Gt02_,G0t2_=G4_old(model,ss[2],lt,div(model.Nt,2))
                
            if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
                println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
                error("WrapTime=$lt ")
            end
            #####################################################################

            UpdateSCEELayer!(rng,view(ss[1],:,lt),view(ss[2],:,lt),G1,G2,A,B,model,UPD,SCEE,λ)

            ##------------------------------------------------------------------------
            tmpO+=(A.detg/B.detg)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------

            if  any(model.nodes .== lt) 
                idx+=1
                BM_F!(tmpN,tmpNN,view(BMs1,:,:,idx-1),model,ss[1],idx-1)
                BMinv_F!(tmpN,tmpNN,view(BMsinv1,:,:,idx-1),model,ss[1],idx-1)
                BM_F!(tmpN,tmpNN,view(BMs2,:,:,idx-1),model,ss[2],idx-1)
                BMinv_F!(tmpN,tmpNN,view(BMsinv2,:,:,idx-1),model,ss[2],idx-1)
                for i in idx:max(Θidx,idx)
                    # println("update BR i=",i)
                    mul!(tmpNn, view(BMs1,:,:,i-1), view(BRMs1,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    copyto!(view(BRMs1,:,:,i) , tmpNn)
                    # ---------------------------------------------------------------
                    mul!(tmpNn, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    copyto!(view(BRMs2,:,:,i) , tmpNn)
                end

                for i in idx-1:-1:min(Θidx,idx)
                    # println("update BL i=",i)
                    mul!(tmpnN,view(BLMs1,:,:,i+1),view(BMs1,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau, ns)
                    copyto!(view(BLMs1,:,:,i) , tmpnN)
                    # ---------------------------------------------------------------
                    mul!(tmpnN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau, ns)
                    copyto!(view(BLMs2,:,:,i) , tmpnN)
                end
                get_ABGM!(G1,G2,A,B,SCEE,model.nodes,idx,"Forward")
            end

        end

        # println("\n ----------------reverse update ----------------")

        for lt in model.Nt:-1:1
            
            #####################################################################
            Gt1_,G01_,Gt01_,G0t1_=G4_old(model,ss[1],lt,div(model.Nt,2))
            Gt2_,G02_,Gt02_,G0t2_=G4_old(model,ss[2],lt,div(model.Nt,2))
            if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
                println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
                error("WrapTime=$lt ")
            end
            #####################################################################
            
            UpdateSCEELayer!(rng,view(ss[1],:,lt),view(ss[2],:,lt),G1,G2,A,B,model,UPD,SCEE,λ)

            ##------------------------------------------------------------------------
            tmpO+=(A.detg/B.detg)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------

            if  any(model.nodes.== (lt-1)) 
                idx-=1
                BM_F!(tmpN,tmpNN,view(BMs1,:,:,idx),model,ss[1],idx)
                BM_F!(tmpN,tmpNN,view(BMs2,:,:,idx),model,ss[2],idx)
                BMinv_F!(tmpN,tmpNN,view(BMsinv1,:,:,idx),model,ss[1],idx)
                BMinv_F!(tmpN,tmpNN,view(BMsinv2,:,:,idx),model,ss[2],idx)
                for i in idx:-1:min(Θidx,idx)
                    # println("update BL i=",i)
                    mul!(tmpnN,view(BLMs1,:,:,i+1),view(BMs1,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau, ns)
                    copyto!(view(BLMs1,:,:,i) , tmpnN)

                    mul!(tmpnN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    LAPACK.gerqf!(tmpnN,tau)
                    LAPACK.orgrq!(tmpnN, tau, ns)
                    copyto!(view(BLMs2,:,:,i) , tmpnN)
                end
                for i in idx+1:max(Θidx,idx)
                    # println("update BR i=",i)
                    mul!(tmpNn, view(BMs1,:,:,i-1), view(BRMs1,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    copyto!(view(BRMs1,:,:,i) , tmpNn)

                    mul!(tmpNn, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    copyto!(view(BRMs2,:,:,i) , tmpNn)
                end
                get_ABGM!(G1,G2,A,B,SCEE,model.nodes,idx,"Backward")
            else
                @inbounds @simd for iii in 1:Ns
                    @fastmath tmpN[iii] = cis(- model.α *model.η[ss[1][iii, lt]] ) 
                    @fastmath tmpN_[iii] = cis(- model.α *model.η[ss[2][iii, lt]] ) 
                end
    
                WrapKV!(tmpNN,model.eK,model.eKinv,tmpN,Gt01,"Backward", "L")
                WrapKV!(tmpNN,model.eK,model.eKinv,tmpN_,Gt02,"Backward", "L")
                WrapKV!(tmpNN,model.eK,model.eKinv,tmpN,Gt1,"Backward", "B")
                WrapKV!(tmpNN,model.eK,model.eKinv,tmpN_,Gt2,"Backward", "B")
                WrapKV!(tmpNN,model.eK,model.eKinv,tmpN,G0t1,"Backward", "R")
                WrapKV!(tmpNN,model.eK,model.eKinv,tmpN_,G0t2,"Backward", "R")
            end
            
        end

        O[loop+1]=tmpO/counter
        tmpO=0.0
        counter=0
    end

    return ss
end


"""
    Return det(Tau)
    Update s1 and overwrite a , b , Tau.
        a = G0t1(:,subidx) ⋅ r
        b = Gt01[subidx,:] ⋅ (2G02-I) ⋅ gmInv
    with r ≡ inv(r) ⋅ Δ
    Warning : G02 here !!!  Gt01,G0t1
    ------------------------------------------------------------------------------
"""
function get_abTau1!(A::AreaBuffer_,UPD::UpdateBuffer_,G0,Gt0,G0t)
    copyto!(A.NN, view(G0,A.index,A.index))
    lmul!(2.0, A.NN)
    for i in diagind(A.NN)
        A.NN[i] -= 1
    end
    mul!(A.zN,view(Gt0,UPD.subidx,A.index),A.NN)
    mul!(A.b, A.zN, A.gmInv)
    mul!(A.a,view(G0t,A.index,UPD.subidx),UPD.r)
    mul!(A.Tau, A.b, A.a)
    for i in diagind(A.Tau)
        A.Tau[i] += 1
    end
    ans=abs2(det(A.Tau))
    return ans
end


"""
    Return det(Tau)
    Update s2 and overwrite a , b , Tau.
        a = (2G01-I) ⋅ G0t2(:,subidx)
        b = r ⋅ Gt02[subidx,:] ⋅ gmInv
    with r ≡ inv(r) ⋅ Δ
    Warning : G01 here !!!  Gt02,G0t2
    ------------------------------------------------------------------------------
"""
function get_abTau2!(A::AreaBuffer_,UPD::UpdateBuffer_,G0,Gt0,G0t)
    copyto!(A.NN , view(G0,A.index,A.index))
    lmul!(2.0, A.NN)
    for i in diagind(A.NN)
        A.NN[i] -= 1
    end
    mul!(A.a,A.NN,view(G0t,A.index,UPD.subidx))
    mul!(A.zN,UPD.r,view(Gt0,UPD.subidx,A.index))
    mul!(A.b,A.zN,A.gmInv)
    mul!(A.Tau, A.b, A.a)
    for i in diagind(A.Tau)
        A.Tau[i] += 1
    end
    ans=abs2(det(A.Tau))
    return ans
end

"""
    Update G4 by overwriting Gt,G0,Gt0,G0t with r and subidx.
        G0 = G0 + G0t(:,subidx) ⋅ r ⋅ Gt0[subidx,:]
        Gt0 = Gt0 + r ⋅ Gt[subidx,:] ⋅ G0t(:,subidx)

        G0t = G0t - G0t(:,subidx) ⋅ r ⋅ (I-Gt[subidx,:])
        Gt = Gt - r ⋅ Gt[subidx,:] ⋅ r ⋅ (I-Gt[subidx,:])
    with r ≡ inv(r) ⋅ Δ
    ------------------------------------------------------------------------------
"""
function G4update!(SCEE::SCEEBuffer_,UPD::UpdateBuffer_,G::G4Buffer_)
    mul!(SCEE.zN,UPD.r,view(G.Gt0,UPD.subidx,:))   # useful for GΘ,GτΘ
    mul!(SCEE.NN, view(G.G0t,:,UPD.subidx),SCEE.zN)
    axpy!(1.0, SCEE.NN, G.G0)
    mul!(SCEE.NN, view(G.Gt,:,UPD.subidx),SCEE.zN)
    axpy!(1.0, SCEE.NN, G.Gt0)

    mul!(SCEE.zN,UPD.r,view(G.Gt,UPD.subidx,:))
    lmul!(-1.0,SCEE.zN)
    axpy!(1.0,UPD.r,view(SCEE.zN,:,UPD.subidx))   # useful for GΘτ,Gτ
    mul!(SCEE.NN, view(G.G0t,:,UPD.subidx),SCEE.zN)
    axpy!(-1.0, SCEE.NN, G.G0t)
    mul!(SCEE.NN, view(G.Gt,:,UPD.subidx),SCEE.zN)
    axpy!(-1.0, SCEE.NN, G.Gt)
end

"""
    Update gmInv with a,b,Tau.
        gmInv = gmInv - gmInv ⋅ a ⋅ inv(Tau) ⋅ b
    Universal for s1 and s2.
    ------------------------------------------------------------------------------
"""
function GMupdate!(A::AreaBuffer_)
    mul!(A.Nz, A.gmInv,A.a )
    # inv22!(A.Tau)
    A.Tau .= inv(A.Tau)
    mul!(A.zN,A.Tau,A.b)
    mul!(A.NN, A.Nz, A.zN)
    axpy!(-1.0, A.NN, A.gmInv)
end

function get_ABGM!(G1::G4Buffer_,G2::G4Buffer_,A::AreaBuffer_,B::AreaBuffer_,SCEE::SCEEBuffer_,nodes,idx,direction::String="Backward")
    G4!(SCEE,G1,nodes,idx,direction)
    G4!(SCEE,G2,nodes,idx,direction)
    GroverMatrix!(A.gmInv,view(G1.G0,A.index,A.index),view(G2.G0,A.index,A.index))
    A.detg=abs2(det(A.gmInv))
    LAPACK.getrf!(A.gmInv, A.ipiv)
    LAPACK.getri!(A.gmInv, A.ipiv)

    GroverMatrix!(B.gmInv,view(G1.G0,B.index,B.index),view(G2.G0,B.index,B.index))
    B.detg=abs2(det(B.gmInv))
    LAPACK.getrf!(B.gmInv, B.ipiv)
    LAPACK.getri!(B.gmInv, B.ipiv)
end 

function UpdateSCEELayer!(rng,s1,s2,G1::G4Buffer_,G2::G4Buffer_,A::AreaBuffer_,B::AreaBuffer_,model::Hubbard_Para_,UPD::UpdateBuffer_,SCEE::SCEEBuffer_,λ)
    for i in axes(s1,1)
        UPD.subidx=[i]

        # update s1
        begin
            sx = rand(rng,  model.samplers_dict[s1[i]])
            p=get_r!(UPD,model.α * (model.η[sx]- model.η[s1[i]]),G1.Gt)
            p*=model.γ[sx]/model.γ[s1[i]]

            detTau_A=get_abTau1!(A,UPD,G2.G0,G1.Gt0,G1.G0t)
            detTau_B=get_abTau1!(B,UPD,G2.G0,G1.Gt0,G1.G0t)

            @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
            if rand(rng)<p
                A.detg*=detTau_A
                B.detg*=detTau_B

                GMupdate!(A)
                GMupdate!(B)
                G4update!(SCEE,UPD,G1)
                s1[i]=sx
            end
        end

        # update ss[2]
        begin
            sx = rand(rng,  model.samplers_dict[s2[i]])
            p=get_r!(UPD,model.α * (model.η[sx]- model.η[s2[i]]),G2.Gt)
            p*=model.γ[sx]/model.γ[s2[i]]

            detTau_A=get_abTau2!(A,UPD,G1.G0,G2.Gt0,G2.G0t)
            detTau_B=get_abTau2!(B,UPD,G1.G0,G2.Gt0,G2.G0t)

            @fastmath p*= (detTau_A)^λ * (detTau_B)^(1-λ)
            if rand(rng)<p
                A.detg*=detTau_A
                B.detg*=detTau_B

                GMupdate!(A)
                GMupdate!(B)
                G4update!(SCEE,UPD,G2)
                s2[i]=sx
            end
        end
    end
end