

function phy_update(path::String, model::_Hubbard_Para, WarmSweeps::Int64, Sweeps::Int64, s::Array{UInt8, 2})
    Ns=model.Ns
    ns=div(Ns, 2)
    NN=length(model.nodes)
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

    mA = mB = nn = R0 = R1 = Ek = C0 = Cmax = 0.0
    counter = 0

    II = Diagonal(ones(ComplexF64,Ns))
    G = Matrix{ComplexF64}(undef ,Ns, Ns)

    # 预分配 BL 和 BR
    BLs = Array{ComplexF64}(undef, ns, Ns,NN)
    BRs = Array{ComplexF64}(undef, Ns, ns,NN)

    # 预分配临时数组
    tmpN = Vector{ComplexF64}(undef, Ns)
    tmpNN = Matrix{ComplexF64}(undef, Ns, Ns)
    DD = Diagonal(tmpN)
    BM = Matrix{ComplexF64}(undef, Ns, Ns)
    tmpNn = Matrix{ComplexF64}(undef, Ns, ns)
    tmpnn = Matrix{ComplexF64}(undef, ns, ns)
    tmpnN = Matrix{ComplexF64}(undef, ns, Ns)
    tmp1N = Matrix{ComplexF64}(undef, 1, Ns)

    view(BRs,:,:,1) .= model.Pt
    view(BLs,:,:,NN) .= model.Pt'
    for idx in NN-1:-1:2
        BM_F!(BM,model, s, idx)
        mul!(tmpnN,view(BLs,:,:,idx+1), BM)
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgqr!(tmpnN, tau, ns)
        view(BLs,:,:,idx) .= tmpnN
        # view(BLs,:,:,idx) .= Matrix(qr!(tmpNn).Q)'
    end
    
    for loop in 1:(Sweeps + WarmSweeps)
        # println("\n Sweep: $loop ")
        
        BM_F!(BM,model, s, 1)
        mul!(tmpnN,view(BLs,:,:,2), BM)
        LAPACK.gerqf!(tmpnN, tau)
        LAPACK.orgrq!(tmpnN, tau, ns)
        view(BLs,:,:,1) .= tmpnN
        
        mul!(tmpnn, view(BLs,:,:,1), view(BRs,:,:,1))
        LAPACK.getrf!(tmpnn,ipiv)
        LAPACK.getri!(tmpnn, ipiv)
        mul!(tmpNn, view(BRs,:,:,1), tmpnn)
        mul!(tmpNN, tmpNn, view(BLs,:,:,1))
        @fastmath G .= II .- tmpNN
        #####################################################################
        # if norm(G-Gτ_old(model,s,0))>1e-6 
        #     error("00000"," Wrap error:  ",norm(G-Gτ_old(model,s,0)))
        # end
        #####################################################################
        idx=1
        for lt in 1:model.Nt
            if any(model.nodes[2:end] .== (lt - 1))
                idx+=1
                BM_F!(BM,model, s, idx - 1)
                mul!(tmpNn, BM, view(BRs,:,:,idx-1))
                LAPACK.geqrf!(tmpNn, tau)
                LAPACK.orgqr!(tmpNn, tau, ns)
                view(BRs,:,:,idx) .= tmpNn
                # BR .= Matrix(qr((BM * BR)).Q)
                mul!(tmpnn, view(BLs,:,:,idx), view(BRs,:,:,idx))
                LAPACK.getrf!(tmpnn,ipiv)
                LAPACK.getri!(tmpnn, ipiv)
                mul!(tmpNn, view(BRs,:,:,idx), tmpnn)
                mul!(tmpNN, tmpNn, view(BLs,:,:,idx))
                @fastmath G .= II .- tmpNN
            end
            s_col = view(model.η,view(s, :, lt))
            @inbounds @simd for iii in 1:Ns
                tmpN[iii] =@fastmath cis( model.α *s_col[iii] ) 
            end
            mul!(tmpNN, model.eK, G)
            mul!(G,tmpNN,model.eKinv)
            DD.diag .= tmpN
            mul!(tmpNN,DD,G)
            conj!(DD)
            mul!(G,tmpNN,DD)
            # G= Diagonal(tmp_D) * model.eK * G * model.eKinv * Diagonal(conj(tmp_D))
            #####################################################################
            # if norm(G-Gτ_old(model,s,lt))>1e-6 
            #     error(lt,"Wrap error:  ",norm(G-Gτ_old(model,s,lt)))
            # end
            #####################################################################

            @inbounds @simd for x in 1:Ns
                sx = rand(rng,  samplers_dict[s[x, lt]])
                @fastmath Δ = cis( model.α * (model.η[sx] - model.η[s[x, lt]])) - 1
                @fastmath r = 1 + Δ * (1 - G[x, x])

                if rand(rng) < @fastmath model.γ[sx] / model.γ[s[x, lt]] * abs2(r)
                    view(tmp1N,1, :) .= -view(G,x, :)
                    tmp1N[1, x] += 1
                    mul!(tmpNN, view(G, :, x), tmp1N)
                    axpy!(-Δ/r, tmpNN, G)
                    s[x, lt] = sx
                    ####################################################################
                    # if norm(G-Gτ_old(model,s,lt))>1e-6
                    #     error("asd")
                    # end
                    #####################################################################
                end
            end

            # ---------------------------------------------------------------------------------------------------------
            # if loop > WarmSweeps && abs(lt - model.Nt / 2) <= model.BatchSize
            #     @views tmpNN .= G
            #     if lt > model.Nt / 2
            #         for i in lt:-1:(div(model.Nt, 2) + 1)
            #             for iii in 1:Ns
            #                 tmp_D[iii] = cis( model.α *model.η[s[iii, i]] ) 
            #             end

            #             tmp_G0= model.eKinv * Diagonal(conj(tmp_D))* tmp_G0 *Diagonal(tmp_D) * model.eK
            
            #         end
            #     else
            #         for i in (lt + 1):div(model.Nt, 2)
            #             for iii in 1:Ns
            #                 tmp_D[iii] = cis( model.α *model.η[s[iii, i]] )
            #             end
            #             tmp_G0= Diagonal(tmp_D) * model.eK * tmp_G0 * model.eKinv * Diagonal(conj(tmp_D))
            #         end
            #     end

            #     @views tmp_G0 .= model.HalfeK * tmp_G0 * model.HalfeKinv

            #     tmp = Magnetism(model, tmp_G0)
            #     mA += tmp[1]
            #     mB += tmp[2]
            #     nn += NN(model, tmp_G0)
            #     Ek += EK(model, tmp_G0)
            #     tmp = CzzofSpin(model, tmp_G0)
            #     R0 += tmp[1]
            #     R1 += tmp[2]
            #     C0 += tmp[3]
            #     Cmax += tmp[4]
            #     counter += 1
            # end
            # ---------------------------------------------------------------------------------------------------------

        end

        BM_F!(BM,model, s, NN-1)
        mul!(tmpNn , BM, view(BRs,:,:,NN-1))
        LAPACK.geqrf!(tmpNn, tau)
        LAPACK.orgqr!(tmpNn, tau, ns)
        view(BRs,:,:,NN) .= tmpNn

        mul!(tmpnn, view(BLs,:,:,NN), view(BRs,:,:,NN))
        LAPACK.getrf!(tmpnn,ipiv)
        LAPACK.getri!(tmpnn, ipiv)
        mul!(tmpNn, view(BRs,:,:,NN), tmpnn)
        mul!(tmpNN, tmpNn, view(BLs,:,:,NN))
        @fastmath G .= II .- tmpNN

        for lt in model.Nt-1:-1:1
            if any(model.nodes.== lt)
                BM_F!(BM,model, s, idx)
                mul!(tmpnN,view(BLs,:,:,idx+1),BM)
                LAPACK.gerqf!(tmpnN, tau)
                LAPACK.orgrq!(tmpnN, tau, ns)
                view(BLs,:,:,idx).=tmpnN
                # BL .= Matrix(qr(( BL * BM )').Q)'

                mul!(tmpnn, view(BLs,:,:,idx), view(BRs,:,:,idx))
                LAPACK.getrf!(tmpnn,ipiv)
                LAPACK.getri!(tmpnn, ipiv)
                mul!(tmpNn, view(BRs,:,:,idx), tmpnn)
                mul!(tmpNN, tmpNn, view(BLs,:,:,idx))
                @fastmath G .= II .- tmpNN

                idx-=1
            else
                s_col = view(model.η,view(s, :, lt+1))
                @inbounds @simd for iii in 1:Ns
                    tmpN[iii] =@fastmath cis( model.α *s_col[iii] ) 
                end
                DD.diag .= tmpN
                mul!(tmpNN,G,DD)
                conj!(DD)
                mul!(G, DD , tmpNN)
                mul!(tmpNN,model.eKinv,G)
                mul!(G,tmpNN,model.eK)
                # G= model.eKinv * Diagonal(conj(tmp_D)) * G * Diagonal(tmp_D) * model.eK
            end
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
                    view(tmp1N,1, :) .= -view(G,x, :)
                    tmp1N[1, x] += 1
                    mul!(tmpNN, view(G, :, x), tmp1N)
                    axpy!(-Δ/r, tmpNN, G)
                    s[x, lt] = sx
                    ####################################################################
                    # if norm(G-Gτ_old(model,s,lt))>1e-6
                    #     error("asd")
                    # end
                    #####################################################################
                end
            end

            # ---------------------------------------------------------------------------------------------------------




            # ---------------------------------------------------------------------------------------------------------

        end

        if loop > WarmSweeps
            fid = open("$(path)Phy$(name)_t$(model.t)U$(model.U)size$(model.site)Δt$(model.Δt)Θ$(model.Θ)BS$(model.BatchSize).csv", "a+")
            writedlm(fid, [Ek model.U * nn nn mA mB R0 R1 C0 Cmax] / counter, ',')
            close(fid)
            mA = mB = nn = R0 = R1 = Ek = C0 = Cmax = 0
            counter = 0
        end
    end
    return s
end


# function Poss(model,s)
#     A=model.Pt[:,:]

#     for i in 1:model.Nt
#         D=[model.η[x] for x in s[:,i]]
#         A=diagm(exp.(1im*model.α.*D))*model.eK*A
#     end
#     A=model.Pt'*A

#     ans=det(A)

#     return abs2(ans)
    
# end
