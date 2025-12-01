

function phy_update(path::String, model::_Hubbard_Para, Sweeps::Int64, s::Array{UInt8, 2},record::Bool=false)
    Ns=model.Ns
    ns=div(Ns, 2)
    NN=length(model.nodes)
    Θidx=div(NN,2)+1
    tau = Vector{ComplexF64}(undef, ns)
    ipiv = Vector{LAPACK.BlasInt}(undef, ns)


    name = if model.Lattice=="SQUARE" "□" 
    elseif model.Lattice=="HoneyComb60" "HC" 
    elseif model.Lattice=="HoneyComb120" "HC120" 
    else error("Lattice: $(model.Lattice) is not allowed !") end  

    rng = MersenneTwister(Threads.threadid()+time_ns())
    elements = (1, 2, 3, 4)
    samplers_dict = Dict{UInt8, Random.Sampler}()
    for excluded in elements
        allowed = [i for i in elements if i != excluded]
        samplers_dict[excluded] = Random.Sampler(rng, allowed)
    end

    Ek=Eu=CDW0=CDW1=SDW0=SDW1=0.0
    counter = 0

    G = Matrix{ComplexF64}(undef ,Ns, Ns)

    # 预分配 BL 和 BR
    BLs = Array{ComplexF64}(undef, ns, Ns,NN)
    BRs = Array{ComplexF64}(undef, Ns, ns,NN)

    # 预分配临时数组
    tmpN = Vector{ComplexF64}(undef, Ns)
    tmpNN = Matrix{ComplexF64}(undef, Ns, Ns)
    BM = Matrix{ComplexF64}(undef, Ns, Ns)
    tmpNn = Matrix{ComplexF64}(undef, Ns, ns)
    tmpnn = Matrix{ComplexF64}(undef, ns, ns)
    tmpnN = Matrix{ComplexF64}(undef, ns, Ns)
    tmp1N = Matrix{ComplexF64}(undef, 1, Ns)

    copyto!(view(BRs,:,:,1) , model.Pt)
    transpose!(view(BLs,:,:,NN) , model.Pt)
    for idx in NN-1:-1:1
        BM_F!(tmpN,tmpNN,BM,model, s, idx)
        mul!(tmpnN,view(BLs,:,:,idx+1), BM)
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau, ns)
        copyto!(view(BLs,:,:,idx) , tmpnN)
        # view(BLs,:,:,idx) .= Matrix(qr!(tmpNn).Q)'
    end
    
    idx=1
    get_G!(tmpnn,tmpNn,ipiv,view(BLs,:,:,idx),view(BRs,:,:,idx),G)
    for loop in 1:Sweeps
        # println("\n Sweep: $loop ")
        for lt in 1:model.Nt
            #####################################################################
            # # println("lt=",lt-1)
            # if norm(G-Gτ_old(model,s,lt-1))>1e-6 
            #     error(lt-1,"Wrap error:  ",norm(G-Gτ_old(model,s,lt-1)))
            # end
            #####################################################################

            @inbounds @simd for iii in 1:Ns
                tmpN[iii] =@fastmath cis( model.α *model.η[s[iii,lt]] ) 
            end
            WrapKV!(tmpNN,model.eK,model.eKinv,tmpN,G,"Forward", "B")
            # G= Diagonal(tmp_D) * model.eK * G * model.eKinv * Diagonal(conj(tmp_D))

            @inbounds @simd for x in 1:Ns
                sx = rand(rng,  samplers_dict[s[x, lt]])
                @fastmath Δ = cis( model.α * (model.η[sx] - model.η[s[x, lt]])) - 1
                @fastmath r = 1 + Δ * (1 - G[x, x])

                if rand(rng) < @fastmath model.γ[sx] / model.γ[s[x, lt]] * abs2(r)
                    r=Δ / r
                    Gupdate!(tmpNN, tmp1N, x, r, G)
                    s[x, lt] = sx
                    ####################################################################
                    # if norm(G-Gτ_old(model,s,lt))>1e-6
                    #     error("asd")
                    # end
                    #####################################################################
                end
            end
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

                get_G!(tmpnn,tmpNn,ipiv,view(BLs,:,:,idx),view(BRs,:,:,idx),G)
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
            # print("-")
            # if norm(G-Gτ_old(model,s,lt))>1e-6 
            #     error(lt," Wrap error:  ",norm(G-Gτ_old(model,s,lt)))
            # end
            #####################################################################

            @inbounds @simd for x in 1:Ns
                sx = rand(rng,  samplers_dict[s[x, lt]])
                @fastmath Δ = cis( model.α * (model.η[sx] - model.η[s[x, lt]])) - 1
                @fastmath r = 1 + Δ * (1 - G[x, x])

                if rand(rng) < @fastmath model.γ[sx] / model.γ[s[x, lt]] * abs2(r)
                    r=Δ / r
                    Gupdate!(tmpNN, tmp1N, x, r, G)
                    s[x, lt] = sx
                    ####################################################################
                    # if norm(G-Gτ_old(model,s,lt))>1e-6
                    #     error("asd")
                    # end
                    #####################################################################
                end
            end
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

                get_G!(tmpnn,tmpNn,ipiv,view(BLs,:,:,idx),view(BRs,:,:,idx),G)
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
function get_G!(tmpnn,tmpNn,ipiv,BL,BR,G)
    mul!(tmpnn, BL,BR)
    LAPACK.getrf!(tmpnn,ipiv)
    LAPACK.getri!(tmpnn, ipiv)
    mul!(tmpNn, BR, tmpnn)
    mul!(G, tmpNn, BL)
    lmul!(-1.0,G)
    for i in diagind(G)
        G[i]+=1
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

function Gupdate!(tmpNN::Matrix{ComplexF64},tmp1N::Matrix{ComplexF64},x::Int64,r::ComplexF64,G::Matrix{ComplexF64})
    view(tmp1N,1, :) .= .-view(G,x, :)
    tmp1N[1, x] += 1
    mul!(tmpNN, view(G, :, x), tmp1N)
    axpy!(-r, tmpNN, G)
end


function phy_measure(model::_Hubbard_Para,lt,s,G,tmpNN,tmpN)
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