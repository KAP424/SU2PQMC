# attractive-U and repulsive-U get the same S_2


function ctrl_SCEEicr(path::String,model::_Hubbard_Para,indexA::Vector{Int64},indexB::Vector{Int64},Sweeps::Int64,λ::Float64,Nλ::Int64,ss::Vector{Matrix{UInt8}},record)
    Ns=model.Ns
    ns=div(model.Ns, 2)
    NN=length(model.nodes)
    tau = Vector{ComplexF64}(undef, ns)
    global ipiv = Vector{LAPACK.BlasInt}(undef, ns)
    ipivA = Vector{LAPACK.BlasInt}(undef, length(indexA))
    ipivB = Vector{LAPACK.BlasInt}(undef, length(indexB))


    name = if model.Lattice=="SQUARE" "□" 
        elseif model.Lattice=="HoneyComb60" "HC" 
        elseif model.Lattice=="HoneyComb120" "HC120" 
        else error("Lattice: $(model.Lattice) is not allowed !") end    
    file="$(path)SCEEicr$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)N$(Nλ)BS$(model.BatchSize).csv"
    
    atexit() do
        if record
            open(file, "a") do io
                lock(io)
                writedlm(io, O', ',')
                unlock(io)
            end
        end
        writedlm("$(path)ss/SS$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)λ$(Int(round(Nλ*λ))).csv", [ss[1] ss[2]],",")
    end
    
    Gt1= Matrix{ComplexF64}(undef ,model.Ns, model.Ns)
    Gt2= Matrix{ComplexF64}(undef ,model.Ns, model.Ns)
    G01= Matrix{ComplexF64}(undef ,model.Ns, model.Ns)
    G02= Matrix{ComplexF64}(undef ,model.Ns, model.Ns)
    Gt01= Matrix{ComplexF64}(undef ,model.Ns, model.Ns)
    Gt02= Matrix{ComplexF64}(undef ,model.Ns, model.Ns)
    G0t1= Matrix{ComplexF64}(undef ,model.Ns, model.Ns)
    G0t2= Matrix{ComplexF64}(undef ,model.Ns, model.Ns)
    gmInv_A=Matrix{ComplexF64}(undef ,length(indexA),length(indexA))
    gmInv_B=Matrix{ComplexF64}(undef ,length(indexB),length(indexB))
    detg_A=detg_B=0 

    b_A= Matrix{ComplexF64}(undef ,1,length(indexA))
    tmp1A= Matrix{ComplexF64}(undef ,1,length(indexA))
    tmpA1= Matrix{ComplexF64}(undef ,length(indexA),1)
    a_A= Matrix{ComplexF64}(undef ,length(indexA),1)

    b_B= Matrix{ComplexF64}(undef ,1,length(indexB))
    tmp1B= Matrix{ComplexF64}(undef ,1,length(indexB))
    tmpB1= Matrix{ComplexF64}(undef ,length(indexB),1)
    a_B= Matrix{ComplexF64}(undef ,length(indexB),1)

    # 预分配临时数组
    global tmpN = Vector{ComplexF64}(undef, Ns)
    tmpN2 = Vector{ComplexF64}(undef, Ns)
    global tmpNN = Matrix{ComplexF64}(undef, Ns, Ns)
    global tmpNN2 = Matrix{ComplexF64}(undef, Ns, Ns)
    global tmpNn = Matrix{ComplexF64}(undef, Ns, ns)
    global tmpnn = Matrix{ComplexF64}(undef, ns, ns)
    tmpnN = Matrix{ComplexF64}(undef, ns, Ns)
    tmp1N = Matrix{ComplexF64}(undef ,1, Ns)
    tmpAA = Matrix{ComplexF64}(undef ,length(indexA),length(indexA))
    tmpBB = Matrix{ComplexF64}(undef ,length(indexB),length(indexB))

    rng=MersenneTwister(Threads.threadid()+time_ns())
    elements=(1, 2, 3, 4)
    samplers_dict = Dict{UInt8, Random.Sampler}()
    for excluded in elements
        allowed = [i for i in elements if i != excluded]
        samplers_dict[excluded] = Random.Sampler(rng, allowed)
    end

    tmpO=0
    counter=0
    O=zeros(Float64,Sweeps+1)
    O[1]=λ

    global II=Diagonal(ones(ComplexF64,model.Ns))

    BMs1=Array{ComplexF64}(undef,model.Ns,model.Ns,NN-1)  # Number_of_BM*Ns*Ns
    BMs2=Array{ComplexF64}(undef,model.Ns,model.Ns,NN-1)  # Number_of_BM*Ns*Ns
    BMsinv1=Array{ComplexF64}(undef,model.Ns,model.Ns,NN-1)  # Number_of_BM*Ns*Ns
    BMsinv2=Array{ComplexF64}(undef,model.Ns,model.Ns,NN-1)  # Number_of_BM*Ns*Ns

    for idx in 1:NN-1
        BM_F!(view(BMs1,:, : , idx),model,ss[1],idx)
        BM_F!(view(BMs2,:,:,idx),model,ss[2],idx)
        BMinv_F!(view(BMsinv1,:,:,idx),model,ss[1],idx)
        BMinv_F!(view(BMsinv2,:,:,idx),model,ss[2],idx)
    end

    BLMs1=Array{ComplexF64}(undef,ns,model.Ns,NN)
    BRMs1=Array{ComplexF64}(undef,model.Ns,ns,NN)
    view(BLMs1,:,:,NN) .= model.Pt'
    view(BRMs1,:,:,1) .= model.Pt
    
    BLMs2=Array{ComplexF64}(undef,ns,model.Ns,NN)
    BRMs2=Array{ComplexF64}(undef,model.Ns,ns,NN)
    view(BLMs2,:,:,NN) .= model.Pt'
    view(BRMs2,:,:,1) .= model.Pt

    for i in 1:NN-1
        mul!(tmpnN,view(BLMs1,:,:,NN-i+1),view(BMs1,:,:,NN-i))
        tmpNn.=tmpnN'
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau, ns)
        view(BLMs1,:,:,NN-i) .= tmpNn'
        
        mul!(tmpNn, view(BMs1,:,:,i), view(BRMs1,:,:,i))
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau, ns)
        view(BRMs1,:,:,i+1) .= tmpNn
        # ---------------------------------------------------------------
        mul!(tmpnN,view(BLMs2,:,:,NN-i+1),view(BMs2,:,:,NN-i))
        tmpNn.=tmpnN'
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau, ns)
        view(BLMs2,:,:,NN-i) .= tmpNn'

        mul!(tmpNn, view(BMs2,:,:,i), view(BRMs2,:,:,i))
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau, ns)
        view(BRMs2,:,:,i+1) .= tmpNn

    end
    Θidx=div(NN,2)+1


    for loop in 1:Sweeps
        println("\n ====== Sweep $loop / $Sweeps ======")

        for lt in 1:model.Nt
            if  any(model.nodes.==(lt-1)) 
                # println("\n Wrap Time: $lt")
                idx= (lt==1) ? 2 : findfirst(model.nodes .== (lt-1))
                BM_F!(view(BMs1,:,:,idx-1),model,ss[1],idx-1)
                BMinv_F!(view(BMsinv1,:,:,idx-1),model,ss[1],idx-1)
                BM_F!(view(BMs2,:,:,idx-1),model,ss[2],idx-1)
                BMinv_F!(view(BMsinv2,:,:,idx-1),model,ss[2],idx-1)
                for i in idx:max(Θidx,idx)
                    # println("update BR i=",i)
                    mul!(tmpNn, view(BMs1,:,:,i-1), view(BRMs1,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BRMs1,:,:,i) .= tmpNn
                    # ---------------------------------------------------------------
                    mul!(tmpNn, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BRMs2,:,:,i) .= tmpNn
                end

                for i in idx:-1:min(Θidx,idx)-1
                    # println("update BL i=",i)
                    mul!(tmpnN,view(BLMs1,:,:,i+1),view(BMs1,:,:,i))
                    tmpNn.=tmpnN'
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BLMs1,:,:,i) .= tmpNn'
                    # ---------------------------------------------------------------
                    mul!(tmpnN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    tmpNn.=tmpnN'
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BLMs2,:,:,i) .= tmpNn'
                end
                idx=findfirst(model.nodes .== (lt-1))
                G4!(Gt1,G01,Gt01,G0t1,model.nodes,idx,BLMs1,BRMs1,BMs1,BMsinv1)
                G4!(Gt2,G02,Gt02,G0t2,model.nodes,idx,BLMs2,BRMs2,BMs2,BMsinv2)
                GroverMatrix!(gmInv_A,view(G01,indexA,indexA),view(G02,indexA,indexA))
                detg_A=abs2(det(gmInv_A))
                LAPACK.getrf!(gmInv_A,ipivA)
                LAPACK.getri!(gmInv_A, ipivA)
                GroverMatrix!(gmInv_B,view(G01,indexB,indexB),view(G02,indexB,indexB))
                detg_B=abs2(det(gmInv_B))
                LAPACK.getrf!(gmInv_B,ipivB)
                LAPACK.getri!(gmInv_B, ipivB)
                # #####################################################################
                # GM_A_=GroverMatrix_old(view(G01,indexA,indexA),view(G02,indexA,indexA))
                # GM_B_=GroverMatrix_old(view(G01,indexB,indexB),view(G02,indexB,indexB))
                # GM_A_=inv(GM_A_)
                # GM_B_=inv(GM_B_)
                # if norm(gmInv_A-GM_A_)+norm(gmInv_B-GM_B_)>1e-8
                #     println(norm(gmInv_A-GM_A_)," ",norm(gmInv_B-GM_B_))
                #     error("WrapTime=$lt GM_A wrong!")
                # end
                # #####################################################################
            end
            @inbounds @simd for iii in 1:Ns
                @fastmath tmpN[iii] = cis( model.α *model.η[ss[1][iii, lt]] ) 
                @fastmath tmpN2[iii] = cis( model.α *model.η[ss[2][iii, lt]] ) 
            end
            mul!(tmpNN, model.eK, Gt1)
            mul!(Gt1,tmpNN,model.eKinv)
            mul!(tmpNN,Diagonal(tmpN),Gt1)
            conj!(tmpN)
            mul!(Gt1,tmpNN,Diagonal(tmpN))

            mul!(tmpNN, model.eK, Gt2)
            mul!(Gt2,tmpNN,model.eKinv)
            mul!(tmpNN,Diagonal(tmpN2),Gt2)
            conj!(tmpN2)
            mul!(Gt2,tmpNN,Diagonal(tmpN2))

            if lt==div(model.Nt,2)+1
                mul!(tmpNN, model.eKinv, Diagonal(tmpN))
                tmpNN2 .= G01
                for i in diagind(tmpNN2)
                    tmpNN2[i] -= 1
                end
                mul!(G0t1, tmpNN2 ,tmpNN)
                conj!(tmpN)
                mul!(tmpNN, model.eK, G01)
                mul!(Gt01,Diagonal(tmpN),tmpNN)
                
                mul!(tmpNN, model.eKinv, Diagonal(tmpN2))
                tmpNN2 .= G02
                for i in diagind(tmpNN2)
                    tmpNN2[i] -= 1
                end
                mul!(G0t2, tmpNN2 ,tmpNN)
                conj!(tmpN2)
                mul!(tmpNN, model.eK, G02)
                mul!(Gt02,Diagonal(tmpN2),tmpNN)

                # Gt01=diagm(exp.(1im*model.α.*D1))*model.eK*G01
                # Gt02=diagm(exp.(1im*model.α.*D2))*model.eK*G02
                # G0t1=-(II-G01)*model.eKinv*diagm(exp.(-1im*model.α.*D1))
                # G0t2=-(II-G02)*model.eKinv*diagm(exp.(-1im*model.α.*D2))
            else
                mul!(tmpNN, G0t1, model.eKinv)
                mul!(G0t1, tmpNN , Diagonal(tmpN))
                conj!(tmpN)
                mul!(tmpNN, model.eK, Gt01)
                mul!(Gt01,Diagonal(tmpN),tmpNN)

                mul!(tmpNN, G0t2, model.eKinv)
                mul!(G0t2, tmpNN , Diagonal(tmpN2))
                conj!(tmpN2)
                mul!(tmpNN, model.eK, Gt02)
                mul!(Gt02,Diagonal(tmpN2),tmpNN)
                # G0t1=G0t1*model.eKinv*diagm(exp.(-1im*model.α.*D1))
                # G0t2=G0t2*model.eKinv*diagm(exp.(-1im*model.α.*D2))
                # Gt01=diagm(exp.(1im*model.α.*D1))*model.eK*Gt01
                # Gt02=diagm(exp.(1im*model.α.*D2))*model.eK*Gt02
            end

            #####################################################################
            # Gt1_,G01_,Gt01_,G0t1_=G4_old(model,ss[1],lt,div(model.Nt,2))
            # Gt2_,G02_,Gt02_,G0t2_=G4_old(model,ss[2],lt,div(model.Nt,2))
                
            # if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
            #     println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
            #     error("WrapTime=$lt ")
            # end
            #####################################################################


            for x in 1:model.Ns
                # @views tmpAA .= 2* G02[indexA,indexA] .- IA
                tmpAA .= view(G02,indexA,indexA)
                lmul!(2.0, tmpAA)
                for i in diagind(tmpAA)
                    tmpAA[i] -= 1
                end
                b_A .= view(Gt01,x:x,indexA)
                mul!(tmp1A, b_A, tmpAA)
                mul!(b_A, tmp1A, gmInv_A)
                a_A .= view(G0t1,indexA,x)
                Tau_A=(b_A*a_A)[1,1]

                # @views tmpBB .= 2*G02[indexB,indexB] .- IB
                tmpBB .= view(G02,indexB,indexB)
                lmul!(2.0, tmpBB)
                for i in diagind(tmpBB)
                    tmpBB[i] -= 1
                end
                b_B .= view(Gt01,x:x,indexB)
                mul!(tmp1B, b_B, tmpBB)
                mul!(b_B, tmp1B, gmInv_B)
                a_B .= view(G0t1,indexB,x)
                Tau_B=(b_B*a_B)[1,1]

                # b_A=transpose(Gt01[x,view(A,:)]) *(2*G02[view(A,:),view(A,:)]-IA)*gmInv_A
                # a_A=G0t1[view(A,:),x]
                # Tau_A=b_A*a_A
                # b_B=transpose(Gt01[x,view(B,:)]) *(2*G02[view(B,:),view(B,:)]-IB)*gmInv_B
                # a_B=G0t1[view(B,:),x]
                # Tau_B=b_B*a_B
                
                sx1=rand(rng,samplers_dict[ss[1][x,lt]])
                
                @fastmath Δ1=cis(model.α*(model.η[sx1]-model.η[ss[1][x,lt]]))-1
                @fastmath r1=1+Δ1*(1-Gt1[x,x])

                @fastmath p=model.γ[sx1]/model.γ[ss[1][x,lt]]*abs2(r1+Δ1*Tau_A)^λ*abs2(r1+Δ1*Tau_B)^(1-λ)
                
                if rand(rng)<p
                    rho_A=Δ1/(r1+Tau_A*Δ1)
                    detg_A*=abs2(1+Δ1/r1*Tau_A)
                    rho_B=Δ1/(r1+Tau_B*Δ1)
                    detg_B*=abs2(1+Δ1/r1*Tau_B)

                    mul!(tmpA1, gmInv_A,a_A )
                    mul!(tmpAA, tmpA1, b_A)
                    axpy!(-rho_A, tmpAA, gmInv_A)
                    mul!(tmpB1, gmInv_B,a_B )
                    mul!(tmpBB, tmpB1, b_B)
                    axpy!(-rho_B, tmpBB, gmInv_B)
                    # gmInv_A-=rho_A* ( gmInv_A*a_A .* b_A)
                    # gmInv_B-=rho_B* ( gmInv_B*a_B .* b_B)

                    mul!(tmpNN, view(G0t1,:,x),view(Gt01,x:x,:))
                    axpy!(Δ1/r1, tmpNN, G01)
                    mul!(tmpNN, view(Gt1,:,x),view(Gt01,x:x,:))
                    axpy!(Δ1/r1, tmpNN, Gt01)
                    tmp1N[1, :] .= -view(Gt1,x, :)
                    tmp1N[1, x] += 1
                    mul!(tmpNN, view(G0t1,:,x),tmp1N)
                    axpy!(-Δ1/r1, tmpNN, G0t1)
                    mul!(tmpNN, view(Gt1,:,x),tmp1N)
                    axpy!(-Δ1/r1, tmpNN, Gt1)

                    ss[1][x,lt]=sx1
                    # G01+=Δ1/r1* (G0t1[:,x] .* transpose(Gt01[x,:]))
                    # Gt01+=Δ1/r1* (Gt1[:,x] .* transpose(Gt01[x,:]))
                    # G0t1-=Δ1/r1* (G0t1[:,x] .* transpose( (II-Gt1)[x,:] ) )
                    # Gt1-=Δ1/r1* (Gt1[:,x] .* transpose( (II-Gt1)[x,:]) )         
                    #####################################################################
                    # print('-')
                    # Gt1_,G01_,Gt01_,G0t1_=G4_old(model,ss[1],lt,div(model.Nt,2))
                    # GM_A_=GroverMatrix(G01_[indexA,indexA],G02[indexA,indexA])
                    # gmInv_A_=inv(GM_A_)
                    # GM_B_=GroverMatrix(G01_[indexB,indexB],G02[indexB,indexB])
                    # gmInv_B_=inv(GM_B_)
                    # detg_A_=abs2(det(GM_A_))
                    # detg_B_=abs2(det(GM_B_))
                    # if norm(Gt1-Gt1_)+norm(G01-G01_)+norm(Gt01-Gt01_)+norm(G0t1-G0t1_)+
                    #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                    #     println('\n',norm(Gt1-Gt1_),'\n',norm(G01-G01_),'\n',norm(Gt01-Gt01_),'\n',norm(G0t1-G0t1_))
                    #     println(norm(gmInv_A_-gmInv_A),' ',norm(gmInv_B-gmInv_B_),"\n",abs(detg_A-detg_A_),' ',abs(detg_B-detg_B_))
                    #     error("$lt  $x:,,,asdasdasd")
                    # end
                    #####################################################################
                end


                tmp1A .= view(Gt02,x:x,indexA)
                mul!(b_A, tmp1A, gmInv_A)
                # @views tmpAA .= 2*G01[indexA,indexA] .- IA
                tmpAA .= view(G01,indexA,indexA)
                lmul!(2.0, tmpAA)
                for i in diagind(tmpAA)
                    tmpAA[i] -= 1
                end
                mul!(a_A, tmpAA, view(G0t2,indexA,x))
                Tau_A=(b_A*a_A)[1,1]


                tmp1B .= view(Gt02,x:x,indexB)
                mul!(b_B, tmp1B, gmInv_B)
                # @views tmpBB .= 2*G01[indexB,indexB] .- IB
                tmpBB .= view(G01,indexB,indexB)
                lmul!(2.0, tmpBB)
                for i in diagind(tmpBB)
                    tmpBB[i] -= 1
                end
                mul!(a_B, tmpBB, view(G0t2,indexB,x))
                Tau_B=(b_B*a_B)[1,1]

                # b_A=transpose(Gt02[x,view(A,:)]) *gmInv_A
                # a_A=(2*G01[view(A,:),view(A,:)]-IA)*G0t2[view(A,:),x]
                # b_B=transpose(Gt02[x,view(B,:)]) *gmInv_B
                # a_B=(2*G01[view(B,:),view(B,:)]-IB)*G0t2[view(B,:),x]
                # Tau_B=b_B*a_B

                sx2=rand(rng,samplers_dict[ss[2][x,lt]])

                @fastmath Δ2=cis(model.α*(model.η[sx2]-model.η[ss[2][x,lt]]))-1
                @fastmath r2=(1+Δ2*(1-Gt2[x,x]))
                @fastmath p=model.γ[sx2]/model.γ[ss[2][x,lt]]*abs2(r2+Δ2*Tau_A)^λ*abs2(r2+Δ2*Tau_B)^(1-λ)

                if rand(rng)<p
                    rho_A=Δ2/(r2+Tau_A*Δ2)
                    detg_A*=abs2(1+Δ2/r2*Tau_A)
                    rho_B=Δ2/(r2+Tau_B*Δ2)
                    detg_B*=abs2(1+Δ2/r2*Tau_B)

                    mul!(tmpA1, gmInv_A,a_A )
                    mul!(tmpAA, tmpA1, b_A)
                    axpy!(-rho_A, tmpAA, gmInv_A)
                    mul!(tmpB1, gmInv_B,a_B )
                    mul!(tmpBB, tmpB1, b_B)
                    axpy!(-rho_B, tmpBB, gmInv_B)
                    # gmInv_A-=rho_A* ( gmInv_A*a_A .* b_A)
                    # gmInv_B-=rho_B* ( gmInv_B*a_B .* b_B)

                    mul!(tmpNN, view(G0t2,:,x),view(Gt02,x:x,:))
                    axpy!(Δ2/r2, tmpNN, G02)
                    mul!(tmpNN, view(Gt2,:,x),view(Gt02,x:x,:))
                    axpy!(Δ2/r2, tmpNN, Gt02)
                    tmp1N[1, :] .= -view(Gt2,x, :)
                    tmp1N[1, x] += 1
                    mul!(tmpNN, view(G0t2,:,x),tmp1N)
                    axpy!(-Δ2/r2, tmpNN, G0t2)
                    mul!(tmpNN, view(Gt2,:,x),tmp1N)
                    axpy!(-Δ2/r2, tmpNN, Gt2)
                    ss[2][x,lt]=sx2
                    # G02+=Δ2/r2* (G0t2[:,x] .* transpose( Gt02[x,:]))
                    # Gt02+=Δ2/r2* (Gt2[:,x] .* transpose( Gt02[x,:]))
                    # G0t2-=Δ2/r2* (G0t2[:,x] .* transpose( (II-Gt2)[x,:]))
                    # Gt2-=Δ2/r2* (Gt2[:,x] .* transpose( (II-Gt2)[x,:])   )      
                    # #####################################################################
                    # print('*')
                    # Gt2_,G02_,Gt02_,G0t2_=G4_old(model,ss[2],lt,div(model.Nt,2))
                    # GM_A_=GroverMatrix(G01[indexA,indexA],G02_[indexA,indexA])
                    # gmInv_A_=inv(GM_A_)
                    # GM_B_=GroverMatrix(G01[indexB,indexB],G02_[indexB,indexB])
                    # gmInv_B_=inv(GM_B_)
                    # detg_A_=abs2(det(GM_A_))
                    # detg_B_=abs2(det(GM_B_))

                    # if norm(Gt2-Gt2_)+norm(G02-G02_)+norm(Gt02-Gt02_)+norm(G0t2-G0t2_)+
                    #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                    #     println('\n',norm(Gt2-Gt2_),'\n',norm(G02-G02_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t2-G0t2_))
                    #     error("$lt  $x:,,,asdasdasd")
                    # end
                    #####################################################################
                end

            end

            ##------------------------------------------------------------------------

            tmpO+=(detg_A/detg_B)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
        end

        G4!(Gt1,G01,Gt01,G0t1,model.nodes,NN,BLMs1,BRMs1,BMs1,BMsinv1)
        G4!(Gt2,G02,Gt02,G0t2,model.nodes,NN,BLMs2,BRMs2,BMs2,BMsinv2)
        GroverMatrix!(gmInv_A,view(G01,indexA,indexA),view(G02,indexA,indexA))
        detg_A=abs2(det(gmInv_A))
        LAPACK.getrf!(gmInv_A,ipivA)
        LAPACK.getri!(gmInv_A, ipivA)
        GroverMatrix!(gmInv_B,view(G01,indexB,indexB),view(G02,indexB,indexB))
        detg_B=abs2(det(gmInv_B))
        LAPACK.getrf!(gmInv_B,ipivB)
        LAPACK.getri!(gmInv_B, ipivB)

        # println("\n ----------------reverse update ----------------")

        for lt in model.Nt:-1:1
            if  any(model.nodes.==lt) 
                idx= (lt==model.nodes[end]) ? NN : findfirst(model.nodes .== lt)+1
                # println("\n Wrap Time: $lt")
                BM_F!(view(BMs1,:,:,idx-1),model,ss[1],idx-1)
                BM_F!(view(BMs2,:,:,idx-1),model,ss[2],idx-1)
                BMinv_F!(view(BMsinv1,:,:,idx-1),model,ss[1],idx-1)
                BMinv_F!(view(BMsinv2,:,:,idx-1),model,ss[2],idx-1)

                for i in idx-1:-1:min(Θidx,idx)-1
                    # println("update BL i=",i)
                    mul!(tmpnN,view(BLMs1,:,:,i+1),view(BMs1,:,:,i))
                    tmpNn.=tmpnN'
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BLMs1,:,:,i) .= tmpNn'

                    mul!(tmpnN,view(BLMs2,:,:,i+1),view(BMs2,:,:,i))
                    tmpNn.=tmpnN'
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BLMs2,:,:,i) .= tmpNn'
                    # view(BLMs1,:,:i) .=BLMs1[i+1,:,:]*BMs1[i,:,:]
                    # BLMs1[i,:,:]=Matrix(qr(BLMs1[i,:,:]').Q)'
                    # BLMs2[i,:,:]=BLMs2[i+1,:,:]*BMs2[i,:,:]
                    # BLMs2[i,:,:]=Matrix(qr(BLMs2[i,:,:]').Q)'
                end
                for i in idx:max(Θidx,idx)
                    # println("update BR i=",i)
                    mul!(tmpNn, view(BMs1,:,:,i-1), view(BRMs1,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BRMs1,:,:,i) .= tmpNn

                    mul!(tmpNn, view(BMs2,:,:,i-1), view(BRMs2,:,:,i-1))
                    LAPACK.geqrf!(tmpNn,tau)
                    LAPACK.orgqr!(tmpNn, tau, ns)
                    view(BRMs2,:,:,i) .= tmpNn
                    # BRMs1[i,:,:]=BMs1[i-1,:,:]*BRMs1[i-1,:,:]
                    # BRMs2[i,:,:]=BMs2[i-1,:,:]*BRMs2[i-1,:,:]
                    # BRMs1[i,:,:]=Matrix(qr(BRMs1[i,:,:]).Q)
                    # BRMs2[i,:,:]=Matrix(qr(BRMs2[i,:,:]).Q)
                end
                idx=findfirst(model.nodes .== lt)
                G4!(Gt1,G01,Gt01,G0t1,model.nodes,idx,BLMs1,BRMs1,BMs1,BMsinv1)
                G4!(Gt2,G02,Gt02,G0t2,model.nodes,idx,BLMs2,BRMs2,BMs2,BMsinv2)
                GroverMatrix!(gmInv_A,view(G01,indexA,indexA),view(G02,indexA,indexA))
                detg_A=abs2(det(gmInv_A))
                LAPACK.getrf!(gmInv_A,ipivA)
                LAPACK.getri!(gmInv_A, ipivA)
                GroverMatrix!(gmInv_B,view(G01,indexB,indexB),view(G02,indexB,indexB))
                detg_B=abs2(det(gmInv_B))
                LAPACK.getrf!(gmInv_B,ipivB)
                LAPACK.getri!(gmInv_B, ipivB)

                #####################################################################
                # println("--------Test BMs and BMinvs--------")
                # BMs=zeros(ComplexF64,NN-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns
                # BMinvs=zeros(ComplexF64,NN-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns

                # for idxx in axes(BMs,1)
                #     BMs[:,:,idxx]=BM_F(model,ss[1],idxx)
                #     BMinvs[:,:,idxx]=BMinv_F(model,ss[1],idxx)
                #     println(norm(BMs[:,:,idxx]-BMs1[:,:,idxx]),",",norm(BMinvs[:,:,idxx]-BMsinv1[:,:,idxx]))
                # end
                # println("--------Test BLMs and BRMs --------")
                # BLMs=zeros(ComplexF64,NN,ns,model.Ns)
                # BRMs=zeros(ComplexF64,NN,model.Ns,ns)
                # BLMs[end,:,:]=model.Pt'[:,:]
                # BRMs[1,:,:]=model.Pt[:,:]
                # for i in axes(BMs,1)
                #     BLMs[end-i,:,:]=Matrix(qr( (BLMs[end-i+1,:,:]*BMs[end-i+1,:,:])' ).Q)'
                #     BRMs[i+1,:,:]=Matrix(qr( BMs[i,:,:]*BRMs[i,:,:] ).Q)
                #     println(norm(BLMs1[end-i,:,:]-BLMs[end-i,:,:]),",",norm(BRMs1[i+1,:,:]-BRMs[i+1,:,:]))
                # end
                #####################################################################

            else
                @inbounds @simd for iii in 1:Ns
                    @fastmath tmpN[iii] = cis( model.α *model.η[ss[1][iii, lt+1]] ) 
                    @fastmath tmpN2[iii] = cis( model.α *model.η[ss[2][iii, lt+1]] ) 
                end
                mul!(tmpNN,Gt1,Diagonal(tmpN))
                conj!(tmpN)
                mul!(Gt1,Diagonal(tmpN),tmpNN)
                mul!(tmpNN,model.eKinv,Gt1)
                mul!(Gt1,tmpNN,model.eK)

                mul!(tmpNN,Gt2,Diagonal(tmpN2))
                conj!(tmpN2)
                mul!(Gt2,Diagonal(tmpN2),tmpNN)
                mul!(tmpNN,model.eKinv,Gt2)
                mul!(Gt2,tmpNN,model.eK)

                # D1=[model.η[x] for x in ss[1][:,lt+1]]
                # D2=[model.η[x] for x in ss[2][:,lt+1]]
                # Gt1=model.eKinv*diagm(exp.(-1im*model.α.*D1)) *Gt1* diagm(exp.(1im*model.α.*D1))*model.eK 
                # Gt2=model.eKinv*diagm(exp.(-1im*model.α.*D2)) *Gt2* diagm(exp.(1im*model.α.*D2))*model.eK
                
                # if lt==div(model.Nt,2)+1
                #     Gt01=model.eKinv*diagm(exp.(-1im*model.α.*D1))*G01
                #     Gt02=model.eKinv*diagm(exp.(-1im*model.α.*D2))*G02
                #     G0t1=-(II-G01)*diagm(exp.(1im*model.α.*D1))*model.eK
                #     G0t2=-(II-G02)*diagm(exp.(1im*model.α.*D2))*model.eK
                # else
                mul!(tmpNN,Diagonal(tmpN),Gt01)
                conj!(tmpN)
                mul!(Gt01,model.eKinv,tmpNN)
                mul!(tmpNN,G0t1,Diagonal(tmpN))
                mul!(G0t1,tmpNN,model.eK)

                mul!(tmpNN,Diagonal(tmpN2),Gt02)
                conj!(tmpN2)
                mul!(Gt02,model.eKinv,tmpNN)
                mul!(tmpNN,G0t2,Diagonal(tmpN2))
                mul!(G0t2,tmpNN,model.eK)

                # G0t1=G0t1*diagm(exp.(1im*model.α.*D1))*model.eK
                # G0t2=G0t2*diagm(exp.(1im*model.α.*D2))*model.eK
                # Gt01=model.eKinv*diagm(exp.(-1im*model.α.*D1))*Gt01
                # Gt02=model.eKinv*diagm(exp.(-1im*model.α.*D2))*Gt02
            end
            #####################################################################
            # Gt1_,G01_,Gt01_,G0t1_=G4_old(model,ss[1],lt,div(model.Nt,2))
            # Gt2_,G02_,Gt02_,G0t2_=G4_old(model,ss[2],lt,div(model.Nt,2))
            # if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
            #     println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
            #     error("WrapTime=$lt ")
            # end
            #####################################################################
            
            for x in 1:model.Ns
                # @views tmpAA .= 2* G02[indexA,indexA] .- IA
                tmpAA .= view(G02,indexA,indexA)
                lmul!(2.0, tmpAA)
                for i in diagind(tmpAA)
                    tmpAA[i] -= 1
                end
                b_A .= view(Gt01,x:x,indexA)
                mul!(tmp1A, b_A, tmpAA)
                mul!(b_A, tmp1A, gmInv_A)
                a_A .= view(G0t1,indexA,x)
                Tau_A=(b_A*a_A)[1,1]

                # @views tmpBB .= 2*G02[indexB,indexB] .- IB
                tmpBB .= view(G02,indexB,indexB)
                lmul!(2.0, tmpBB)
                for i in diagind(tmpBB)
                    tmpBB[i] -= 1
                end
                b_B .= view(Gt01,x:x,indexB)
                mul!(tmp1B, b_B, tmpBB)
                mul!(b_B, tmp1B, gmInv_B)
                a_B .= view(G0t1,indexB,x)
                Tau_B=(b_B*a_B)[1,1]

                # b_A=transpose(Gt01[x,view(A,:)]) *(2*G02[view(A,:),view(A,:)]-IA)*gmInv_A
                # a_A=G0t1[view(A,:),x]
                # Tau_A=b_A*a_A
                # b_B=transpose(Gt01[x,view(B,:)]) *(2*G02[view(B,:),view(B,:)]-IB)*gmInv_B
                # a_B=G0t1[view(B,:),x]
                # Tau_B=b_B*a_B
                
                sx1=rand(rng,samplers_dict[ss[1][x,lt]])
                
                @fastmath Δ1=cis(model.α*(model.η[sx1]-model.η[ss[1][x,lt]]))-1
                @fastmath r1=1+Δ1*(1-Gt1[x,x])

                @fastmath p=model.γ[sx1]/model.γ[ss[1][x,lt]]*abs2(r1+Δ1*Tau_A)^λ*abs2(r1+Δ1*Tau_B)^(1-λ)
                
                if rand(rng)<p
                    rho_A=Δ1/(r1+Tau_A*Δ1)
                    detg_A*=abs2(1+Δ1/r1*Tau_A)
                    rho_B=Δ1/(r1+Tau_B*Δ1)
                    detg_B*=abs2(1+Δ1/r1*Tau_B)

                    mul!(tmpA1, gmInv_A,a_A )
                    mul!(tmpAA, tmpA1, b_A)
                    axpy!(-rho_A, tmpAA, gmInv_A)
                    mul!(tmpB1, gmInv_B,a_B )
                    mul!(tmpBB, tmpB1, b_B)
                    axpy!(-rho_B, tmpBB, gmInv_B)
                    # gmInv_A-=rho_A* ( gmInv_A*a_A .* b_A)
                    # gmInv_B-=rho_B* ( gmInv_B*a_B .* b_B)

                    mul!(tmpNN, view(G0t1,:,x),view(Gt01,x:x,:))
                    axpy!(Δ1/r1, tmpNN, G01)
                    mul!(tmpNN, view(Gt1,:,x),view(Gt01,x:x,:))
                    axpy!(Δ1/r1, tmpNN, Gt01)
                    tmp1N[1, :] .= -view(Gt1,x, :)
                    tmp1N[1, x] += 1
                    mul!(tmpNN, view(G0t1,:,x),tmp1N)
                    axpy!(-Δ1/r1, tmpNN, G0t1)
                    mul!(tmpNN, view(Gt1,:,x),tmp1N)
                    axpy!(-Δ1/r1, tmpNN, Gt1)

                    ss[1][x,lt]=sx1
                    # G01+=Δ1/r1* (G0t1[:,x] .* transpose(Gt01[x,:]))
                    # Gt01+=Δ1/r1* (Gt1[:,x] .* transpose(Gt01[x,:]))
                    # G0t1-=Δ1/r1* (G0t1[:,x] .* transpose( (II-Gt1)[x,:] ) )
                    # Gt1-=Δ1/r1* (Gt1[:,x] .* transpose( (II-Gt1)[x,:]) )         
                    #####################################################################
                    # print('-')
                    # Gt1_,G01_,Gt01_,G0t1_=G4_old(model,ss[1],lt,div(model.Nt,2))
                    # GM_A_=GroverMatrix(G01_[indexA,indexA],G02[indexA,indexA])
                    # gmInv_A_=inv(GM_A_)
                    # GM_B_=GroverMatrix(G01_[indexB,indexB],G02[indexB,indexB])
                    # gmInv_B_=inv(GM_B_)
                    # detg_A_=abs2(det(GM_A_))
                    # detg_B_=abs2(det(GM_B_))
                    # if norm(Gt1-Gt1_)+norm(G01-G01_)+norm(Gt01-Gt01_)+norm(G0t1-G0t1_)+
                    #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                    #     println('\n',norm(Gt1-Gt1_),'\n',norm(G01-G01_),'\n',norm(Gt01-Gt01_),'\n',norm(G0t1-G0t1_))
                    #     println(norm(gmInv_A_-gmInv_A),' ',norm(gmInv_B-gmInv_B_),"\n",abs(detg_A-detg_A_),' ',abs(detg_B-detg_B_))
                    #     error("$lt  $x:,,,asdasdasd")
                    # end
                    #####################################################################
                end


                tmp1A .= view(Gt02,x:x,indexA)
                mul!(b_A, tmp1A, gmInv_A)
                # @views tmpAA .= 2*G01[indexA,indexA] .- IA
                tmpAA .= view(G01,indexA,indexA)
                lmul!(2.0, tmpAA)
                for i in diagind(tmpAA)
                    tmpAA[i] -= 1
                end
                mul!(a_A, tmpAA, view(G0t2,indexA,x))
                Tau_A=(b_A*a_A)[1,1]


                tmp1B .= view(Gt02,x:x,indexB)
                mul!(b_B, tmp1B, gmInv_B)
                # @views tmpBB .= 2*G01[indexB,indexB] .- IB
                tmpBB .= view(G01,indexB,indexB)
                lmul!(2.0, tmpBB)
                for i in diagind(tmpBB)
                    tmpBB[i] -= 1
                end
                mul!(a_B, tmpBB, view(G0t2,indexB,x))
                Tau_B=(b_B*a_B)[1,1]

                # b_A=transpose(Gt02[x,view(A,:)]) *gmInv_A
                # a_A=(2*G01[view(A,:),view(A,:)]-IA)*G0t2[view(A,:),x]
                # b_B=transpose(Gt02[x,view(B,:)]) *gmInv_B
                # a_B=(2*G01[view(B,:),view(B,:)]-IB)*G0t2[view(B,:),x]
                # Tau_B=b_B*a_B

                sx2=rand(rng,samplers_dict[ss[2][x,lt]])

                @fastmath Δ2=cis(model.α*(model.η[sx2]-model.η[ss[2][x,lt]]))-1
                @fastmath r2=(1+Δ2*(1-Gt2[x,x]))
                @fastmath p=model.γ[sx2]/model.γ[ss[2][x,lt]]*abs2(r2+Δ2*Tau_A)^λ*abs2(r2+Δ2*Tau_B)^(1-λ)

                if rand(rng)<p
                    rho_A=Δ2/(r2+Tau_A*Δ2)
                    detg_A*=abs2(1+Δ2/r2*Tau_A)
                    rho_B=Δ2/(r2+Tau_B*Δ2)
                    detg_B*=abs2(1+Δ2/r2*Tau_B)

                    mul!(tmpA1, gmInv_A,a_A )
                    mul!(tmpAA, tmpA1, b_A)
                    axpy!(-rho_A, tmpAA, gmInv_A)
                    mul!(tmpB1, gmInv_B,a_B )
                    mul!(tmpBB, tmpB1, b_B)
                    axpy!(-rho_B, tmpBB, gmInv_B)
                    # gmInv_A-=rho_A* ( gmInv_A*a_A .* b_A)
                    # gmInv_B-=rho_B* ( gmInv_B*a_B .* b_B)

                    mul!(tmpNN, view(G0t2,:,x),view(Gt02,x:x,:))
                    axpy!(Δ2/r2, tmpNN, G02)
                    mul!(tmpNN, view(Gt2,:,x),view(Gt02,x:x,:))
                    axpy!(Δ2/r2, tmpNN, Gt02)
                    tmp1N[1, :] .= -view(Gt2,x, :)
                    tmp1N[1, x] += 1
                    mul!(tmpNN, view(G0t2,:,x),tmp1N)
                    axpy!(-Δ2/r2, tmpNN, G0t2)
                    mul!(tmpNN, view(Gt2,:,x),tmp1N)
                    axpy!(-Δ2/r2, tmpNN, Gt2)
                    ss[2][x,lt]=sx2
                    # G02+=Δ2/r2* (G0t2[:,x] .* transpose( Gt02[x,:]))
                    # Gt02+=Δ2/r2* (Gt2[:,x] .* transpose( Gt02[x,:]))
                    # G0t2-=Δ2/r2* (G0t2[:,x] .* transpose( (II-Gt2)[x,:]))
                    # Gt2-=Δ2/r2* (Gt2[:,x] .* transpose( (II-Gt2)[x,:])   )      
                    # #####################################################################
                    # print('*')
                    # Gt2_,G02_,Gt02_,G0t2_=G4_old(model,ss[2],lt,div(model.Nt,2))
                    # GM_A_=GroverMatrix(G01[indexA,indexA],G02_[indexA,indexA])
                    # gmInv_A_=inv(GM_A_)
                    # GM_B_=GroverMatrix(G01[indexB,indexB],G02_[indexB,indexB])
                    # gmInv_B_=inv(GM_B_)
                    # detg_A_=abs2(det(GM_A_))
                    # detg_B_=abs2(det(GM_B_))

                    # if norm(Gt2-Gt2_)+norm(G02-G02_)+norm(Gt02-Gt02_)+norm(G0t2-G0t2_)+
                    #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                    #     println('\n',norm(Gt2-Gt2_),'\n',norm(G02-G02_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t2-G0t2_))
                    #     error("$lt  $x:,,,asdasdasd")
                    # end
                    #####################################################################
                end

            end

            ##------------------------------------------------------------------------
            tmpO+=(detg_A/detg_B)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
        end

        O[loop+1]=tmpO/counter
        tmpO=counter=0
    end
    return ss
end
