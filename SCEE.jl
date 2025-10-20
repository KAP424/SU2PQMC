# attractive-U and repulsive-U get the same S_2


function ctrl_SCEEicr(path::String,model::_Hubbard_Para,indexA::Vector{Int64},indexB::Vector{Int64},Sweeps::Int64,λ::Float64,Nλ::Int64,ss::Vector{Matrix{UInt8}},record)
    if model.Lattice=="SQUARE"
        name="□"
    elseif occursin("HoneyComb", model.Lattice)
        name="HC"
    end
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
    
    # Gt1=zeros(ComplexF64,model.Ns,model.Ns)
    # Gt2=zeros(ComplexF64,model.Ns,model.Ns)
    # G01=zeros(ComplexF64,model.Ns,model.Ns)
    # G02=zeros(ComplexF64,model.Ns,model.Ns)
    # Gt01=zeros(ComplexF64,model.Ns,model.Ns)
    # Gt02=zeros(ComplexF64,model.Ns,model.Ns)
    # G0t1=zeros(ComplexF64,model.Ns,model.Ns)
    # G0t2=zeros(ComplexF64,model.Ns,model.Ns)
    # gmInv_A=zeros(ComplexF64,length(indexA),length(indexA))
    # gmInv_B=zeros(ComplexF64,length(indexB),length(indexB))
    # detg_A=detg_B=0

    # rng=MersenneTwister(Threads.threadid()+round(Int,time()*1000))
    rng=MersenneTwister(1)
    elements=(1, 2, 3, 4)

    tmpO=0
    counter=0
    O=zeros(Sweeps+1)
    O[1]=λ

    II=I(model.Ns)
    IA=I(length(indexA))
    IB=I(length(indexB))

    BMs1=zeros(ComplexF64,length(model.nodes)-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns
    BMs2=zeros(ComplexF64,length(model.nodes)-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns
    BMsinv1=zeros(ComplexF64,length(model.nodes)-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns
    BMsinv2=zeros(ComplexF64,length(model.nodes)-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns

    for idx in axes(BMs1,1)
        BMs1[idx,:,:]=BM_F(model,ss[1],idx)
        BMs2[idx,:,:]=BM_F(model,ss[2],idx)
        BMsinv1[idx,:,:]=BMinv_F(model,ss[1],idx)
        BMsinv2[idx,:,:]=BMinv_F(model,ss[2],idx)
    end

    BLMs1=zeros(ComplexF64,length(model.nodes),div(model.Ns,2),model.Ns)
    BRMs1=zeros(ComplexF64,length(model.nodes),model.Ns,div(model.Ns,2))
    BLMs1[end,:,:]=model.Pt'[:,:]
    BRMs1[1,:,:]=model.Pt[:,:]
    BLMs2=zeros(ComplexF64,length(model.nodes),div(model.Ns,2),model.Ns)
    BRMs2=zeros(ComplexF64,length(model.nodes),model.Ns,div(model.Ns,2))
    BLMs2[end,:,:]=model.Pt'[:,:]
    BRMs2[1,:,:]=model.Pt[:,:]

    for i in axes(BMs1,1)
        BLMs1[end-i,:,:]=Matrix(qr( (BLMs1[end-i+1,:,:]*BMs1[end-i+1,:,:])' ).Q)'
        BRMs1[i+1,:,:]=Matrix(qr( BMs1[i,:,:]*BRMs1[i,:,:] ).Q)
        BLMs2[end-i,:,:]=Matrix(qr( (BLMs2[end-i+1,:,:]*BMs2[end-i+1,:,:])' ).Q)'
        BRMs2[i+1,:,:]=Matrix(qr( BMs2[i,:,:]*BRMs2[i,:,:] ).Q)
    end

    Θidx=div(length(model.nodes),2)+1

    for loop in 1:Sweeps
        Gt1,G01,Gt01,G0t1=G4(model.nodes,1,BLMs1,BRMs1,BMs1,BMsinv1)
        Gt2,G02,Gt02,G0t2=G4(model.nodes,1,BLMs2,BRMs2,BMs2,BMsinv2)
        GM_A=GroverMatrix(G01[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
        gmInv_A=inv(GM_A)
        GM_B=GroverMatrix(G01[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
        gmInv_B=inv(GM_B)
        detg_A=abs2(det(GM_A))
        detg_B=abs2(det(GM_B))
        for lt in 1:model.Nt
            if  any(model.nodes[2:end].==(lt-1)) 
                idx=findfirst(model.nodes .== (lt-1))
                println("\n Wrap Time: $lt")
                BMs1[idx-1,:,:]=BM_F(model,ss[1],idx-1)
                BMsinv1[idx-1,:,:]=BMinv_F(model,ss[1],idx-1)
                BMs2[idx-1,:,:]=BM_F(model,ss[2],idx-1)
                BMsinv2[idx-1,:,:]=BMinv_F(model,ss[2],idx-1)

                for i in idx-1:-1:1
                    # println("update BL i=",i)
                    BLMs1[i,:,:]=BLMs1[i+1,:,:]*BMs1[i,:,:]
                    BLMs2[i,:,:]=BLMs2[i+1,:,:]*BMs2[i,:,:]
                    BLMs1[i,:,:]=Matrix(qr(BLMs1[i,:,:]').Q)'
                    BLMs2[i,:,:]=Matrix(qr(BLMs2[i,:,:]').Q)'
                end
                for i in idx:size(BLMs1,1)
                    # println("update BR i=",i)
                    BRMs1[i,:,:]=BMs1[i-1,:,:]*BRMs1[i-1,:,:]
                    BRMs2[i,:,:]=BMs2[i-1,:,:]*BRMs2[i-1,:,:]
                    BRMs1[i,:,:]=Matrix(qr(BRMs1[i,:,:]).Q)
                    BRMs2[i,:,:]=Matrix(qr(BRMs2[i,:,:]).Q)
                end
                Gt1,G01,Gt01,G0t1=G4(model.nodes,idx,BLMs1,BRMs1,BMs1,BMsinv1)
                Gt2,G02,Gt02,G0t2=G4(model.nodes,idx,BLMs2,BRMs2,BMs2,BMsinv2)
                GM_A=GroverMatrix(G01[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
                gmInv_A=inv(GM_A)
                GM_B=GroverMatrix(G01[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
                gmInv_B=inv(GM_B)
                detg_A=abs2(det(GM_A))
                detg_B=abs2(det(GM_B))

                # ------------------------------------------------------------------------------
                for idxx in axes(BMs1,1)
                    println(norm(BM_F(model,ss[1],idxx)-BMs1[idxx,:,:]),",",norm(BMinv_F(model,ss[1],idxx)-BMsinv1[idxx,:,:]))
                    println(norm(BM_F(model,ss[2],idxx)-BMs2[idxx,:,:]),",",norm(BMinv_F(model,ss[2],idxx)-BMsinv2[idxx,:,:]))
                end
                # for i in axes(BLMs1,1)
                #     G1=I(model.Ns)-BRMs1[i,:,:] * inv( BLMs1[i,:,:] * BRMs1[i,:,:] ) * BLMs1[i,:,:]
                #     G2=I(model.Ns)-BRMs2[i,:,:] * inv( BLMs2[i,:,:] * BRMs2[i,:,:] ) * BLMs2[i,:,:]
                #     G1_=Gτ_old(model,ss[1],model.nodes[i])
                #     G2_=Gτ_old(model,ss[2],model.nodes[i])
                #     println(norm(G1-G1_)," ",norm(G2-G2_))
                # end

                # BMs=zeros(ComplexF64,length(model.nodes)-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns
                # BMinvs=zeros(ComplexF64,length(model.nodes)-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns
                # println("--------Test BMs and BMinvs--------")

                # for idxx in axes(BMs,1)
                #     BMs[idxx,:,:]=BM_F(model,ss[1],idxx)
                #     BMinvs[idxx,:,:]=BMinv_F(model,ss[1],idxx)
                #     println(norm(BMs[idxx,:,:]-BMs1[idxx,:,:]),",",norm(BMinvs[idxx,:,:]-BMsinv1[idxx,:,:]))
                # end
                # # println("--------Test BLMs and BRMs 2--------")
                # # BLMs=zeros(ComplexF64,length(model.nodes),div(model.Ns,2),model.Ns)
                # # BRMs=zeros(ComplexF64,length(model.nodes),model.Ns,div(model.Ns,2))
                # # BLMs[end,:,:]=model.Pt'[:,:]
                # # BRMs[1,:,:]=model.Pt[:,:]
                # # for i in axes(BMs,1)
                # #     BLMs[end-i,:,:]=Matrix(qr( (BLMs[end-i+1,:,:]*BMs[end-i+1,:,:])' ).Q)'
                # #     BRMs[i+1,:,:]=Matrix(qr( BMs[i,:,:]*BRMs[i,:,:] ).Q)
                # #     println(norm(BLMs1[end-i,:,:]-BLMs[end-i,:,:]),",",norm(BRMs1[i+1,:,:]-BRMs[i+1,:,:]))
                # # end

                # ------------------------------------------------------------------------------
            end
            D1=[model.η[x] for x in ss[1][:,lt]]
            D2=[model.η[x] for x in ss[2][:,lt]]
            Gt1=diagm(exp.(1im*model.α.*D1))*model.eK *Gt1* model.eKinv*diagm(exp.(-1im*model.α.*D1))
            Gt2=diagm(exp.(1im*model.α.*D2))*model.eK *Gt2* model.eKinv*diagm(exp.(-1im*model.α.*D2))
            
            if lt==div(model.Nt,2)+1
                Gt01=diagm(exp.(1im*model.α.*D1))*model.eK*G01
                Gt02=diagm(exp.(1im*model.α.*D2))*model.eK*G02
                G0t1=-(II-G01)*model.eKinv*diagm(exp.(-1im*model.α.*D1))
                G0t2=-(II-G02)*model.eKinv*diagm(exp.(-1im*model.α.*D2))
            else
                G0t1=G0t1*model.eKinv*diagm(exp.(-1im*model.α.*D1))
                G0t2=G0t2*model.eKinv*diagm(exp.(-1im*model.α.*D2))
                Gt01=diagm(exp.(1im*model.α.*D1))*model.eK*Gt01
                Gt02=diagm(exp.(1im*model.α.*D2))*model.eK*Gt02
            end

            #####################################################################
            Gt1_,G01_,Gt01_,G0t1_=G4_old(model,ss[1],lt,div(model.Nt,2))
            Gt2_,G02_,Gt02_,G0t2_=G4_old(model,ss[2],lt,div(model.Nt,2))
                
            if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
                println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
                error("WrapTime=$lt ")
            end
            #####################################################################


            for x in 1:model.Ns
                b_A=transpose(Gt01[x,indexA[:]]) *(2*G02[indexA[:],indexA[:]]-IA)*gmInv_A
                a_A=G0t1[indexA[:],x]
                Tau_A=b_A*a_A
                
                b_B=transpose(Gt01[x,indexB[:]]) *(2*G02[indexB[:],indexB[:]]-IB)*gmInv_B
                a_B=G0t1[indexB[:],x]
                Tau_B=b_B*a_B
                
                sp=Random.Sampler(rng,[i for i in elements if i != ss[1][x,lt]])
                sx1=rand(rng,sp)
                
                Δ1=exp(1im*model.α*(model.η[sx1]-model.η[ss[1][x,lt]]))-1
                r1=1+Δ1*(1-Gt1[x,x])

                p=model.γ[sx1]/model.γ[ss[1][x,lt]]*abs2(r1+Δ1*Tau_A)^λ*abs2(r1+Δ1*Tau_B)^(1-λ)
                
                if rand(rng)<p
                    rho_A=Δ1/(r1+Tau_A*Δ1)
                    gmInv_A-=rho_A* ( gmInv_A*a_A .* b_A)
                    detg_A*=abs2(1+Δ1/r1*Tau_A)

                    rho_B=Δ1/(r1+Tau_B*Δ1)
                    gmInv_B-=rho_B* ( gmInv_B*a_B .* b_B)
                    detg_B*=abs2(1+Δ1/r1*Tau_B)

                    G01+=Δ1/r1* (G0t1[:,x] .* transpose(Gt01[x,:]))
                    Gt01+=Δ1/r1* (Gt1[:,x] .* transpose(Gt01[x,:]))
                    G0t1-=Δ1/r1* (G0t1[:,x] .* transpose( (II-Gt1)[x,:] ) )
                    Gt1-=Δ1/r1* (Gt1[:,x] .* transpose( (II-Gt1)[x,:]) )         
                    ss[1][x,lt]=sx1
                    #####################################################################
                    print('-')
                    Gt1_,G01_,Gt01_,G0t1_=G4_old(model,ss[1],lt,div(model.Nt,2))
                    GM_A_=GroverMatrix(G01_[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
                    gmInv_A_=inv(GM_A_)
                    GM_B_=GroverMatrix(G01_[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
                    gmInv_B_=inv(GM_B_)
                    detg_A_=abs2(det(GM_A_))
                    detg_B_=abs2(det(GM_B_))

                    if norm(Gt1-Gt1_)+norm(G01-G01_)+norm(Gt01-Gt01_)+norm(G0t1-G0t1_)+
                       norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                        println('\n',norm(Gt1-Gt1_),'\n',norm(G01-G01_),'\n',norm(Gt01-Gt01_),'\n',norm(G0t1-G0t1_))
                        error("$lt  $x:,,,asdasdasd")
                    end
                    #####################################################################
                end



                b_A=transpose(Gt02[x,indexA[:]]) *gmInv_A
                a_A=(2*G01[indexA[:],indexA[:]]-IA)*G0t2[indexA[:],x]
                Tau_A=b_A*a_A

                b_B=transpose(Gt02[x,indexB[:]]) *gmInv_B
                a_B=(2*G01[indexB[:],indexB[:]]-IB)*G0t2[indexB[:],x]
                Tau_B=b_B*a_B

                sp=Random.Sampler(rng,[i for i in elements if i != ss[2][x,lt]])
                sx2=rand(rng,sp)

                Δ2=(exp(1im*model.α*(model.η[sx2]-model.η[ss[2][x,lt]]))-1)
                r2=(1+Δ2*(1-Gt2[x,x]))

                p=model.γ[sx2]/model.γ[ss[2][x,lt]]*abs2(r2+Δ2*Tau_A)^λ*abs2(r2+Δ2*Tau_B)^(1-λ)

                if rand(rng)<p
                    rho_A=Δ2/(r2+Tau_A*Δ2)
                    gmInv_A-=rho_A* ( gmInv_A*a_A .* b_A)
                    detg_A*=abs2(1+Δ2/r2*Tau_A)

                    rho_B=Δ2/(r2+Tau_B*Δ2)
                    gmInv_B-=rho_B* ( gmInv_B*a_B .* b_B)
                    detg_B*=abs2(1+Δ2/r2*Tau_B)

                    G02+=Δ2/r2* (G0t2[:,x] .* transpose( Gt02[x,:]))
                    Gt02+=Δ2/r2* (Gt2[:,x] .* transpose( Gt02[x,:]))
                    G0t2-=Δ2/r2* (G0t2[:,x] .* transpose( (II-Gt2)[x,:]))
                    Gt2-=Δ2/r2* (Gt2[:,x] .* transpose( (II-Gt2)[x,:])   )      
                    ss[2][x,lt]=sx2
                    # #####################################################################
                    print('*')
                    Gt2_,G02_,Gt02_,G0t2_=G4_old(model,ss[2],lt,div(model.Nt,2))
                    GM_A_=GroverMatrix(G01[indexA[:],indexA[:]],G02_[indexA[:],indexA[:]])
                    gmInv_A_=inv(GM_A_)
                    GM_B_=GroverMatrix(G01[indexB[:],indexB[:]],G02_[indexB[:],indexB[:]])
                    gmInv_B_=inv(GM_B_)
                    detg_A_=abs2(det(GM_A_))
                    detg_B_=abs2(det(GM_B_))

                    if norm(Gt2-Gt2_)+norm(G02-G02_)+norm(Gt02-Gt02_)+norm(G0t2-G0t2_)+
                       norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                        println('\n',norm(Gt2-Gt2_),'\n',norm(G02-G02_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t2-G0t2_))
                        error("$lt  $x:,,,asdasdasd")
                    end
                    #####################################################################
                end

            end

            ##------------------------------------------------------------------------

            tmpO+=(detg_A/detg_B)^(1/Nλ)
            counter+=1
            ##------------------------------------------------------------------------
        end

        Gt1,G01,Gt01,G0t1=G4(model.nodes,length(model.nodes),BLMs1,BRMs1,BMs1,BMsinv1)
        Gt2,G02,Gt02,G0t2=G4(model.nodes,length(model.nodes),BLMs2,BRMs2,BMs2,BMsinv2)
        GM_A=GroverMatrix(G01[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
        gmInv_A=inv(GM_A)
        GM_B=GroverMatrix(G01[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
        gmInv_B=inv(GM_B)
        detg_A=abs2(det(GM_A))
        detg_B=abs2(det(GM_B))

        println("\n ----------------reverse update ----------------")

        for lt in model.Nt:-1:1
            if  any(model.nodes.==lt) 
                idx= (lt==model.nodes[end]) ? length(model.nodes) : findfirst(model.nodes .== lt)+1
                println("\n Wrap Time: $lt")
                BMs1[idx-1,:,:]=BM_F(model,ss[1],idx-1)
                BMs2[idx-1,:,:]=BM_F(model,ss[2],idx-1)
                BMsinv1[idx-1,:,:]=BMinv_F(model,ss[1],idx-1)
                BMsinv2[idx-1,:,:]=BMinv_F(model,ss[2],idx-1)
                for i in idx-1:-1:1
                    # println("update BL i=",i)
                    BLMs1[i,:,:]=BLMs1[i+1,:,:]*BMs1[i,:,:]
                    BLMs2[i,:,:]=BLMs2[i+1,:,:]*BMs2[i,:,:]
                    BLMs1[i,:,:]=Matrix(qr(BLMs1[i,:,:]').Q)'
                    BLMs2[i,:,:]=Matrix(qr(BLMs2[i,:,:]').Q)'
                end
                for i in idx:size(BLMs1,1)
                    # println("update BR i=",i)
                    BRMs1[i,:,:]=BMs1[i-1,:,:]*BRMs1[i-1,:,:]
                    BRMs2[i,:,:]=BMs2[i-1,:,:]*BRMs2[i-1,:,:]
                    BRMs1[i,:,:]=Matrix(qr(BRMs1[i,:,:]).Q)
                    BRMs2[i,:,:]=Matrix(qr(BRMs2[i,:,:]).Q)
                end
                idx=findfirst(model.nodes .== lt)
                Gt1,G01,Gt01,G0t1=G4(model.nodes,idx,BLMs1,BRMs1,BMs1,BMsinv1)
                Gt2,G02,Gt02,G0t2=G4(model.nodes,idx,BLMs2,BRMs2,BMs2,BMsinv2)
                GM_A=GroverMatrix(G01[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
                gmInv_A=inv(GM_A)
                GM_B=GroverMatrix(G01[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
                gmInv_B=inv(GM_B)
                detg_A=abs2(det(GM_A))
                detg_B=abs2(det(GM_B))

                #####################################################################
                println("--------Test BMs and BMinvs--------")
                for idxx in axes(BMs1,1)
                    println(norm(BM_F(model,ss[1],idxx)-BMs1[idxx,:,:]),",",norm(BMinv_F(model,ss[1],idxx)-BMsinv1[idxx,:,:]))
                    println(norm(BM_F(model,ss[2],idxx)-BMs2[idxx,:,:]),",",norm(BMinv_F(model,ss[2],idxx)-BMsinv2[idxx,:,:]))
                end
                BMs=zeros(ComplexF64,length(model.nodes)-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns
                BMinvs=zeros(ComplexF64,length(model.nodes)-1,model.Ns,model.Ns)  # Number_of_BM*Ns*Ns
                for idxx in axes(BMs,1)
                    BMs[idxx,:,:]=BM_F(model,ss[1],idxx)
                    BMinvs[idxx,:,:]=BMinv_F(model,ss[1],idxx)
                end
                println("--------Test BLMs and BRMs 2--------")
                BLMs=zeros(ComplexF64,length(model.nodes),div(model.Ns,2),model.Ns)
                BRMs=zeros(ComplexF64,length(model.nodes),model.Ns,div(model.Ns,2))
                BLMs[end,:,:]=model.Pt'[:,:]
                BRMs[1,:,:]=model.Pt[:,:]
                for i in axes(BMs,1)
                    BLMs[end-i,:,:]=Matrix(qr( (BLMs[end-i+1,:,:]*BMs[end-i+1,:,:])' ).Q)'
                    BRMs[i+1,:,:]=Matrix(qr( BMs[i,:,:]*BRMs[i,:,:] ).Q)
                    println(norm(BLMs1[end-i,:,:]-BLMs[end-i,:,:]),",",norm(BRMs1[i+1,:,:]-BRMs[i+1,:,:]))
                end
                #####################################################################

            else
                D1=[model.η[x] for x in ss[1][:,lt+1]]
                D2=[model.η[x] for x in ss[2][:,lt+1]]
                Gt1=model.eKinv*diagm(exp.(-1im*model.α.*D1)) *Gt1* diagm(exp.(1im*model.α.*D1))*model.eK 
                Gt2=model.eKinv*diagm(exp.(-1im*model.α.*D2)) *Gt2* diagm(exp.(1im*model.α.*D2))*model.eK
                
                # if lt==div(model.Nt,2)+1
                #     Gt01=model.eKinv*diagm(exp.(-1im*model.α.*D1))*G01
                #     Gt02=model.eKinv*diagm(exp.(-1im*model.α.*D2))*G02
                #     G0t1=-(II-G01)*diagm(exp.(1im*model.α.*D1))*model.eK
                #     G0t2=-(II-G02)*diagm(exp.(1im*model.α.*D2))*model.eK
                # else
                G0t1=G0t1*diagm(exp.(1im*model.α.*D1))*model.eK
                G0t2=G0t2*diagm(exp.(1im*model.α.*D2))*model.eK
                Gt01=model.eKinv*diagm(exp.(-1im*model.α.*D1))*Gt01
                Gt02=model.eKinv*diagm(exp.(-1im*model.α.*D2))*Gt02
            end
            #####################################################################
            Gt1_,G01_,Gt01_,G0t1_=G4_old(model,ss[1],lt,div(model.Nt,2))
            Gt2_,G02_,Gt02_,G0t2_=G4_old(model,ss[2],lt,div(model.Nt,2))
            if norm(Gt1-Gt1_)+norm(Gt2-Gt2_)+norm(Gt01-Gt01_)+norm(Gt02-Gt02_)+norm(G0t1-G0t1_)+norm(G0t2-G0t2_)>1e-3
                println( norm(Gt1-Gt1_),'\n',norm(Gt2-Gt2_),'\n',norm(Gt01-Gt01_),'\n',norm(Gt02-Gt02_),'\n',norm(G0t1-G0t1_),'\n',norm(G0t2-G0t2_) )
                error("WrapTime=$lt ")
            end
            #####################################################################
                for x in 1:model.Ns
                    b_A=transpose(Gt01[x,indexA[:]]) *(2*G02[indexA[:],indexA[:]]-IA)*gmInv_A
                    a_A=G0t1[indexA[:],x]
                    Tau_A=b_A*a_A
                    
                    b_B=transpose(Gt01[x,indexB[:]]) *(2*G02[indexB[:],indexB[:]]-IB)*gmInv_B
                    a_B=G0t1[indexB[:],x]
                    Tau_B=b_B*a_B
                    
                    sp=Random.Sampler(rng,[i for i in elements if i != ss[1][x,lt]])
                    sx1=rand(rng,sp)
                    
                    Δ1=exp(1im*model.α*(model.η[sx1]-model.η[ss[1][x,lt]]))-1
                    r1=1+Δ1*(1-Gt1[x,x])
    
                    p=model.γ[sx1]/model.γ[ss[1][x,lt]]*abs2(r1+Δ1*Tau_A)^λ*abs2(r1+Δ1*Tau_B)^(1-λ)
    
                    if rand(rng)<p
                        rho_A=Δ1/(r1+Tau_A*Δ1)
                        gmInv_A-=rho_A* ( gmInv_A*a_A .* b_A)
                        detg_A*=abs2(1+Δ1/r1*Tau_A)
    
                        rho_B=Δ1/(r1+Tau_B*Δ1)
                        gmInv_B-=rho_B* ( gmInv_B*a_B .* b_B)
                        detg_B*=abs2(1+Δ1/r1*Tau_B)
    
                        G01+=Δ1/r1* (G0t1[:,x] .* transpose(Gt01[x,:]))
                        Gt01+=Δ1/r1* (Gt1[:,x] .* transpose(Gt01[x,:]))
                        G0t1-=Δ1/r1* (G0t1[:,x] .* transpose( (II-Gt1)[x,:] ) )
                        Gt1-=Δ1/r1* (Gt1[:,x] .* transpose( (II-Gt1)[x,:]) )         
                        ss[1][x,lt]=sx1
    
                        #####################################################################
                        # print('-')
                        # Gt1_,G01_,Gt01_,G0t1_=G4(model,ss[1],lt,div(model.Nt,2))
                        # GM_A_=GroverMatrix(G01_[indexA[:],indexA[:]],G02[indexA[:],indexA[:]])
                        # gmInv_A_=inv(GM_A_)
                        # GM_B_=GroverMatrix(G01_[indexB[:],indexB[:]],G02[indexB[:],indexB[:]])
                        # gmInv_B_=inv(GM_B_)
                        # detg_A_=abs2(det(GM_A_))
                        # detg_B_=abs2(det(GM_B_))
    
                        # if norm(Gt1-Gt1_)+norm(G01-G01_)+norm(Gt01-Gt01_)+norm(G0t1-G0t1_)+
                        #    norm(gmInv_A_-gmInv_A)+norm(gmInv_B-gmInv_B_)+abs(detg_A-detg_A_)+abs(detg_B-detg_B_)>1e-3
                        #     println('\n',norm(Gt1-Gt1_),'\n',norm(G01-G01_),'\n',norm(Gt01-Gt01_),'\n',norm(G0t1-G0t1_))
                        #     error("$lt  $x:,,,asdasdasd")
                        # end
                        #####################################################################
                    end
    
    
    
                    b_A=transpose(Gt02[x,indexA[:]]) *gmInv_A
                    a_A=(2*G01[indexA[:],indexA[:]]-IA)*G0t2[indexA[:],x]
                    Tau_A=b_A*a_A
    
                    b_B=transpose(Gt02[x,indexB[:]]) *gmInv_B
                    a_B=(2*G01[indexB[:],indexB[:]]-IB)*G0t2[indexB[:],x]
                    Tau_B=b_B*a_B
    
                    sp=Random.Sampler(rng,[i for i in elements if i != ss[2][x,lt]])
                    sx2=rand(rng,sp)
    
                    Δ2=(exp(1im*model.α*(model.η[sx2]-model.η[ss[2][x,lt]]))-1)
                    r2=(1+Δ2*(1-Gt2[x,x]))
    
                    p=model.γ[sx2]/model.γ[ss[2][x,lt]]*abs2(r2+Δ2*Tau_A)^λ*abs2(r2+Δ2*Tau_B)^(1-λ)
    
                    if rand(rng)<p
                        rho_A=Δ2/(r2+Tau_A*Δ2)
                        gmInv_A-=rho_A* ( gmInv_A*a_A .* b_A)
                        detg_A*=abs2(1+Δ2/r2*Tau_A)
    
                        rho_B=Δ2/(r2+Tau_B*Δ2)
                        gmInv_B-=rho_B* ( gmInv_B*a_B .* b_B)
                        detg_B*=abs2(1+Δ2/r2*Tau_B)
    
                        G02+=Δ2/r2* (G0t2[:,x] .* transpose( Gt02[x,:]))
                        Gt02+=Δ2/r2* (Gt2[:,x] .* transpose( Gt02[x,:]))
                        G0t2-=Δ2/r2* (G0t2[:,x] .* transpose( (II-Gt2)[x,:]))
                        Gt2-=Δ2/r2* (Gt2[:,x] .* transpose( (II-Gt2)[x,:])   )      
                        ss[2][x,lt]=sx2
                        #####################################################################
                        # print('*')
                        # Gt2_,G02_,Gt02_,G0t2_=G4(model,ss[2],lt,div(model.Nt,2))
                        # GM_A_=GroverMatrix(G01[indexA[:],indexA[:]],G02_[indexA[:],indexA[:]])
                        # gmInv_A_=inv(GM_A_)
                        # GM_B_=GroverMatrix(G01[indexB[:],indexB[:]],G02_[indexB[:],indexB[:]])
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
