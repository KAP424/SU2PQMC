

function phy_update(path::String, model::Hubbard_Para_, Sweeps::Int64, s::Array{UInt8, 2},record::Bool=false)
    global LOCK=ReentrantLock()
    ERROR=1e-6

    NN=length(model.nodes)
    Θidx=div(NN,2)+1
    Phy = PhyBuffer(model.Ns, NN) 
    UPD = UpdateBuffer()

    name = if model.Lattice=="SQUARE" "□" 
    elseif model.Lattice=="HoneyComb60" "HC" 
    elseif model.Lattice=="HoneyComb120" "HC120" 
    else error("Lattice: $(model.Lattice) is not allowed !") end  

    Ns=model.Ns
    ns=div(Ns, 2)
    rng = MersenneTwister(Threads.threadid()+time_ns())

    G = Phy.G
    tau = Phy.tau
    ipiv = Phy.ipiv
    BLs = Phy.BLs
    BRs = Phy.BRs
    tmpN = Phy.N
    tmpNN = Phy.NN
    BM = Phy.BM
    tmpNn = Phy.Nn
    tmpnn = Phy.nn
    tmpnN = Phy.nN

    Ek=Eu=CDW0=CDW1=SDW0=SDW1=0.0
    counter = 0

    copyto!(view(BRs,:,:,1) , model.Pt)
    transpose!(view(BLs,:,:,NN) , model.Pt)
    for idx in NN-1:-1:1
        BM_F!(tmpN,tmpNN,BM,model, s, idx)
        mul!(tmpnN,view(BLs,:,:,idx+1), BM)
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau, ns)
        copyto!(view(BLs,:,:,idx) , tmpnN)
    end
    
    idx=1
    get_G!(tmpnn,tmpnN,ipiv,view(BLs,:,:,idx),view(BRs,:,:,idx),G)
    for loop in 1:Sweeps
        # println("\n Sweep: $loop ")
        for lt in 1:model.Nt
            #####################################################################
            # println("lt=",lt-1)
            if norm(G-Gτ_old(model,s,lt-1))>1e-6 
                error(lt-1,"Wrap error:  ",norm(G-Gτ_old(model,s,lt-1)))
            end
            #####################################################################

            @inbounds @simd for iii in 1:Ns
                tmpN[iii] =@fastmath cis( model.α *model.η[s[iii,lt]] ) 
            end
            WrapKV!(tmpNN,model.eK,model.eKinv,tmpN,G,"Forward", "B")

            UpdatePhyLayer!(rng,view(s,:,lt),model,UPD,Phy)
            # ---------------------------------------------------------------------------------------------------------
            # record physical quantities
            if record && 0<=(Θidx-idx)<=1
                tmp=phy_measure(model,lt,s,G,tmpNN,tmpN)
                counter += 1
                Ek+=tmp[1]
                Eu+=tmp[2]
                CDW0+=tmp[3]
                CDW1+=tmp[4]
                SDW0+=tmp[5]
                SDW1+=tmp[6]
            end
            # ---------------------------------------------------------------------------------------------------------
    
            if any(model.nodes .== lt )
                idx+=1
                BM_F!(tmpN,tmpNN,BM,model, s, idx - 1)
                mul!(tmpNn, BM, view(BRs,:,:,idx-1))
                LAPACK.geqrf!(tmpNn, tau)
                LAPACK.orgqr!(tmpNn, tau, ns)
                copyto!(view(BRs,:,:,idx) , tmpNn)

                copyto!(tmpNN , G)

                get_G!(tmpnn,tmpnN,ipiv,view(BLs,:,:,idx),view(BRs,:,:,idx),G)
                #####################################################################
                axpy!(-1.0, G, tmpNN)  
                if norm(tmpNN)>1e-8
                    println("Warning for Batchsize Wrap Error : $(norm(tmpNN))")
                end
                #####################################################################
            end
        end

        for lt in model.Nt:-1:1
            #####################################################################
            print("-")
            if norm(G-Gτ_old(model,s,lt))>1e-6 
                error(lt," Wrap error:  ",norm(G-Gτ_old(model,s,lt)))
            end
            #####################################################################

            UpdatePhyLayer!(rng,view(s,:,lt),model,UPD,Phy)
            # ---------------------------------------------------------------------------------------------------------
            # record physical quantities
            if record && 0<=(idx-Θidx)<=1
                tmp=phy_measure(model,lt,s,G,tmpNN,tmpN)
                counter += 1
                Ek+=tmp[1]
                Eu+=tmp[2]
                CDW0+=tmp[3]
                CDW1+=tmp[4]
                SDW0+=tmp[5]
                SDW1+=tmp[6]
            end
            # ---------------------------------------------------------------------------------------------------------
            @inbounds @simd for iii in 1:Ns
                tmpN[iii] =@fastmath cis(-model.α *model.η[s[iii,lt]] ) 
            end
            WrapKV!(tmpNN,model.eK,model.eKinv,tmpN,G,"Backward", "B")
            
            if any(model.nodes.== (lt-1))
                # println("idx=",idx," lt=",lt-1)
                idx-=1
                BM_F!(tmpN,tmpNN,BM,model, s, idx)
                mul!(tmpnN,view(BLs,:,:,idx+1),BM)
                LAPACK.gerqf!(tmpnN, tau)
                LAPACK.orgrq!(tmpnN, tau, ns)
                copyto!(view(BLs,:,:,idx),tmpnN)

                get_G!(tmpnn,tmpnN,ipiv,view(BLs,:,:,idx),view(BRs,:,:,idx),G)
            end
        end

        if record
            fid = open("$(path)/Phy$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv", "a+")
            writedlm(fid, [Ek Eu CDW0 CDW1 SDW0 SDW1] / counter, ',')
            close(fid)
            Ek=Eu=CDW0=CDW1=SDW0=SDW1=0.0
            counter = 0
        end
    end
    return s
end

"""
    No Return. Overwrite G 
        G = I - BR ⋅ inv(BL ⋅ BR) ⋅ BL 
    ------------------------------------------------------------------------------
"""
function get_G!(tmpnn,tmpnN,ipiv,BL,BR,G)
    # 目标: 计算 G = I - BR * inv(BL * BR) * BL，避免显式求逆 (getri!)
    # 步骤:
    # 1. tmpnn ← M = BL * BR
    # 2. LU 分解 tmpnn 得到 pivot ipiv
    # 3. tmpnN ← BL (右端)，求解 M * X = BL 得 X = inv(M)*BL  (使用 getrs!)
    # 4. G ← BR * X
    # 5. G ← I - G
    # 数值优势: 避免显式逆，提升稳定性与性能，减少 FLOPs。
    mul!(tmpnn, BL, BR)                 # tmpnn = M = BL*BR
    LAPACK.getrf!(tmpnn, ipiv)          # LU 分解 (in-place)
    tmpnN .= BL                         # 右端初始化: RHS = BL
    LAPACK.getrs!('N', tmpnn, ipiv, tmpnN) # 解 M * X = BL, 结果写回 tmpNn
    mul!(G, BR, tmpnN)                  # G = BR * inv(M) * BL
    lmul!(-1.0, G)                      # G = - G
    @inbounds for i in diagind(G)       # G = I - BR * inv(M) * BL
        G[i] += 1.0
    end
end

function WrapKV!(tmpNN,eK,eKinv,D,G,direction,LR)
    if direction=="Forward"
        if LR=="L"
            mul!(tmpNN, eK, G)
            mul!(G,Diagonal(D),tmpNN)
        elseif LR=="R"
            mul!(tmpNN, G,eKinv)
            mul!(G,tmpNN , Diagonal(D))
        elseif LR=="B"
            mul!(tmpNN, eK, G)
            mul!(G,tmpNN,eKinv)
            mul!(tmpNN,Diagonal(D),G)
            conj!(D)
            mul!(G,tmpNN,Diagonal(D))
        end
    elseif direction=="Backward"
        if LR=="L"
            mul!(tmpNN,Diagonal(D),G)
            mul!(G, eKinv, tmpNN)
        elseif LR=="R"
            mul!(tmpNN,G,Diagonal(D))
            mul!(G, tmpNN,eK)
        elseif LR=="B"
            mul!(tmpNN,Diagonal(D),G)
            conj!(D)
            mul!(G,tmpNN,Diagonal(D))
            mul!(tmpNN, eKinv, G)
            mul!(G,tmpNN,eK)
        end
    end
end

"""
    No Return. Overwrite G = G - G · inv(r) ⋅ Δ · (I-G)
    ------------------------------------------------------------------------------
"""
function Gupdate!(Phy::PhyBuffer_,UPD::UpdateBuffer_)
    mul!(Phy.zN,UPD.r,view(Phy.G,UPD.subidx,:))
    lmul!(-1.0,Phy.zN)
    axpy!(1.0,UPD.r,view(Phy.zN,:,UPD.subidx))   # useful for GΘτ,Gτ
    mul!(Phy.NN, view(Phy.G,:,UPD.subidx),Phy.zN)
    axpy!(-1.0, Phy.NN, Phy.G)
end

function phy_measure(model::Hubbard_Para_,lt,s,G,tmpNN,tmpN)
    """
    (Ek,Ev,R0,R1)    
    """
    G0=G[:,:]
    CDW0=CDW1=SDW0=SDW1=0.0
    
    if lt>model.Nt/2
        for t in lt:-1:div(model.Nt,2)+1
            for i in axes(s,1)
                tmpN[i]=cis(-model.α * model.η[s[i,t]])
            end
            WrapKV!(tmpNN,model.eK,model.eKinv,tmpN,G0,"Backward", "B")
        end
    else
        for t in lt+1:div(model.Nt,2)
            for i in axes(s,1)
                tmpN[i]=cis(model.α * model.η[s[i,t]])
            end
            WrapKV!(tmpNN,model.eK,model.eKinv,tmpN,G0,"Forward","B")
        end
    end
    #####################################################################
    # if norm(G0-Gτ_old(model,s,div(model.Nt,2)))>1e-7
    #     error("record error lt=$(lt) : $(norm(G0-Gτ_old(model,s,div(model.Nt,2))))")
    # end
    #####################################################################
    mul!(tmpNN,model.HalfeK,G0)
    mul!(G0,tmpNN,model.HalfeKinv)
    # G0=model.HalfeK* G0 *model.HalfeKinv

    Ek=-2*model.t*real(sum(model.K.*G0))
    Eu=0
    for i in 1:model.Ns
        Eu+=(1-G0[i,i])*adjoint(G0[i,i])
    end
    @assert imag(Eu)<1e-8 "Eu get imaginary part! $(Eu)" 
    Eu=model.U*real(Eu)

    # @assert sum(abs.(imag(diag(G0))))<1e-8 "G0 diag get imaginary part! $(norm(imag(diag(G0))))"

    if occursin("HoneyComb", model.Lattice)
        for rx in 1:model.site[1]
            for ry in 1:model.site[2]
                tmp1=0.0     #<up up> <down down>  
                tmp2=0.0     #<up down> <down up>
                for ix in 1:model.site[1]
                    for iy in 1:model.site[2]
                        idx1=xy_i(model.Lattice,model.site,ix,iy)-1
                        idx2=xy_i(model.Lattice,model.site,mod1(ix+rx,model.site[1]),mod1(iy+ry,model.site[2]))-1
                        
                        delta = idx1==idx2 ? 1 : 0

                        tmp1+=( (1-G0[idx1,idx1]) * (1-G0[idx2,idx2]) + (delta-G0[idx1,idx2])*G0[idx2,idx1] ) 
                            + adjoint( G0[idx1,idx1]*G0[idx2,idx2] + G0[idx1,idx2]*(delta-G0[idx2,idx1]) )  # <i_A j_A>
                        tmp1+=( (1-G0[idx1+1,idx1+1]) * (1-G0[idx2+1,idx2+1]) + (delta-G0[idx1+1,idx2+1])*G0[idx2+1,idx1+1] )
                            + adjoint( G0[idx1+1,idx1+1]*G0[idx2+1,idx2+1] + G0[idx1+1,idx2+1]*(delta-G0[idx2+1,idx1+1]) )  # <i_B j_B>
                        tmp1-=( (1-G0[idx1+1,idx1+1])*(1-G0[idx2,idx2])-G0[idx1+1,idx2]*G0[idx2,idx1+1] )
                            + adjoint( G0[idx1+1,idx1+1]*G0[idx2,idx2]-G0[idx1+1,idx2]*G0[idx2,idx1+1] )  # <i_B j_A>
                        tmp1-=( (1-G0[idx1,idx1])*(1-G0[idx2+1,idx2+1])-G0[idx1,idx2+1]*G0[idx2+1,idx1] )
                            + adjoint( G0[idx1,idx1]*G0[idx2+1,idx2+1]-G0[idx1,idx2+1]*G0[idx2+1,idx1] )  # <i_A j_B>

                        tmp2+= (1-G0[idx1,idx1])*adjoint(G0[idx2,idx2]) 
                            + adjoint(G0[idx1,idx1])*(1-G0[idx2,idx2])    # <i_A j_A>
                        tmp2+= (1-G0[idx1+1,idx1+1])*adjoint(G0[idx2+1,idx2+1])
                            + adjoint(G0[idx1+1,idx1+1])*(1-G0[idx2+1,idx2+1])    # <i_B j_B>
                        tmp2-= (1-G0[idx1,idx1])*adjoint(G0[idx2+1,idx2+1]) 
                            + adjoint(G0[idx1,idx1])*(1-G0[idx2+1,idx2+1])    # <i_A j_B>
                        tmp2-= (1-G0[idx1+1,idx1+1])*adjoint(G0[idx2,idx2]) 
                            + adjoint(G0[idx1+1,idx1+1])*(1-G0[idx2,idx2])    # <i_B j_A>
                    end
                end
                # @assert imag(tmp1)<1e-8 "tmp1 get imaginary part! $(imag(tmp1))"
                # @assert imag(tmp2)<1e-8 "tmp2 get imaginary part! $(imag(tmp2))"
                CDW0+=tmp1+tmp2
                CDW1+=cos(2*π/model.site[1]*rx+2*π/model.site[2]*ry )* (tmp1+tmp2)
                SDW0+=tmp1-tmp2
                SDW1+=cos(2*π/model.site[1]*rx+2*π/model.site[2]*ry )* (tmp1-tmp2)
            end
        end
        
        # @assert imag(CDW0)+imag(CDW1)+imag(SDW0)+imag(SDW1)<1e-8 "struct factor get imaginary part! $(imag(CDW0)+imag(CDW1)+imag(SDW0)+imag(SDW1))"
        CDW0=real(CDW0)/model.Ns^2
        CDW1=real(CDW1)/model.Ns^2
        SDW0=real(SDW0)/model.Ns^2
        SDW1=real(SDW1)/model.Ns^2
    elseif model.Lattice=="SQUARE"
        for rx in 1:model.site[1]
            for ry in 1:model.site[2]
                tmp=0
                for ix in 1:model.site[1]
                    for iy in 1:model.site[2]
                        idx1=ix+(iy-1)*model.site[1]
                        idx2=mod1(rx+ix,model.site[1])+mod((ry+iy-1),model.site[2])*model.site[1]
                        tmp+=(1-G0[idx1,idx1])*(1-G0[idx2,idx2])-G0[idx1,idx2]*G0[idx2,idx1]
                    end
                end
                tmp/=prod(model.site)
                R0+=tmp*cos(π*(rx+ry))
                R1+=cos(π*(rx+ry)+2*π/model.site[1]*rx+2*π/model.site[2]*ry )*tmp
            end
        end
    end
    return Ek,Eu,CDW0,CDW1,SDW0,SDW1
end

function UpdatePhyLayer!(rng,s,model::Hubbard_Para_,UPD::UpdateBuffer_,Phy::PhyBuffer_)
    for i in eachindex(s)
        UPD.subidx = [i]
        sx = rand(rng,  model.samplers_dict[s[i]])
        p = get_r!(UPD,model.α*( model.η[sx] - model.η[s[i]]), Phy.G)
        p*= model.γ[sx] / model.γ[s[i]]
        if rand(rng) < p
            Gupdate!(Phy,UPD)
            s[i] = sx
        end
    end
end