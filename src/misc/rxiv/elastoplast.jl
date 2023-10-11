@views function deform!(mpD,meD)
    # init mesh quantities to zero
    meD.ΔJ .= 0.0
    # set identity matrix
    ID   = Matrix(I,meD.nD,meD.nD)
    # action
    for p ∈ 1:mpD.nmp
        # get nodal incremental displacement
        iD            = mpD.p2n[p,:]
        Δun           = meD.u[iD,:]
        # compute incremental deformation gradient
        mpD.ΔF[:,:,p].= ID+vcat(mpD.ϕ∂ϕ[p,:,2]'*Δun,mpD.ϕ∂ϕ[p,:,3]'*Δun)'
        mpD.ΔJ[p]     = det(mpD.ΔF[:,:,p])
        # update deformation gradient
        mpD.F[:,:,p] .= mpD.ΔF[:,:,p]*mpD.F[:,:,p]
        # update material point's volume and domain length
        mpD.J[p]      = det(mpD.F[:,:,p])
        mpD.V[p]      = mpD.J[p]*mpD.V0[p]
        mpD.l[p,:]   .= mpD.J[p]^(1/meD.nD).*mpD.l0[p,:] 
        # accumulation
        meD.ΔJ[iD]  .+= mpD.ϕ∂ϕ[p,:].*mpD.m[p].*mpD.ΔJ[p]  
    end 
    # compute nodal deformation determinant
    @threads for no ∈ 1:meD.nno[meD.nD+1]
        if meD.m[no]>0.0 
            meD.ΔJ[no] = meD.ΔJ[no]/meD.m[no]
        end
    end
    # compute determinant Jbar 
    @threads for p ∈ 1:mpD.nmp
        mpD.ΔFbar[:,:,p].= mpD.ΔF[:,:,p].*(((mpD.ϕ∂ϕ[p,:,1]'*meD.ΔJ[mpD.p2n[p,:]])/mpD.ΔJ[p]).^(1/meD.nD))
    end
    return nothing
end
@views function mutate(ϵ,dim,type)
    if type == "tensor"
        ϵ   = [     ϵ[1] 0.5*ϵ[4];
                0.5*ϵ[4]     ϵ[2]]
    elseif type == "voigt"
        ϵ   = vcat(ϵ[1,1],ϵ[2,2],0.0,2*ϵ[1,2])
    end
    return ϵ
end
# For volumetric locking, F-bar method is used, see DOI: 10.1002/nag.3599
@views function elast!(mpD,Del,isΔFbar,fwrkDeform)
    if isΔFbar
        ΔF = copy(mpD.ΔFbar)
    else
        ΔF = copy(mpD.ΔF)
    end
    if fwrkDeform == "finite"
        @threads for p ∈ 1:mpD.nmp
            # compute logarithmic strain tensor
            ϵ          = mutate(mpD.ϵ[:,p],2,"tensor")
            λ,n        = eigvals(ϵ),eigvecs(ϵ)
            b          = n*diagm(exp.(2*λ))*n'
            bt         = ΔF[:,:,p]*b*ΔF[:,:,p]'
            λ,n        = eigvals(bt),eigvecs(bt)
            ϵt         = 0.5.*(n*diagm(log.(λ))*n')
            mpD.ϵ[:,p].= mutate(ϵt,2,"voigt")
            # krichhoff stress tensor
            mpD.τ[:,p].= (Del*mpD.ϵ[:,p]) 
        end
    elseif fwrkDeform == "infinitesimal"
        @threads for p ∈ 1:mpD.nmp
            # calculate elastic strains
            mpD.ϵ[1,p] = ΔF[1,1,p]-1.0
            mpD.ϵ[2,p] = ΔF[2,2,p]-1.0
            mpD.ϵ[4,p] = ΔF[1,2,p]+ΔF[2,1,p]
            mpD.ω[p]   = 0.5*(ΔF[2,1,p]-ΔF[1,2,p])
            # update stresses
            σxx0       = mpD.σ[1,p]
            σyy0       = mpD.σ[2,p]
            σxy0       = mpD.σ[4,p]
            mpD.σ[1,p]+= (Del[1,1]*mpD.ϵ[1,p]+Del[1,2]*mpD.ϵ[2,p]+Del[1,4]*mpD.ϵ[4,p])+2.0*σxy0*mpD.ω[p]
            mpD.σ[2,p]+= (Del[2,1]*mpD.ϵ[1,p]+Del[2,2]*mpD.ϵ[2,p]+Del[2,4]*mpD.ϵ[4,p])-2.0*σxy0*mpD.ω[p]
            mpD.σ[3,p]+= (Del[3,1]*mpD.ϵ[1,p]+Del[3,2]*mpD.ϵ[2,p]+Del[3,4]*mpD.ϵ[4,p])
            mpD.σ[4,p]+= (Del[4,1]*mpD.ϵ[1,p]+Del[4,2]*mpD.ϵ[2,p]+Del[4,4]*mpD.ϵ[4,p])+(σyy0-σxx0)*mpD.ω[p]     
        end        
    end
    return nothing
end
@views function elastoplast!(mpD::NamedTuple,meD::NamedTuple,cmParam::NamedTuple,cmType::String,isΔFbar::Bool,fwrkDeform::String,plastOn::Bool)
    # get def. & logarithmic strains
    deform!(mpD,meD)
    # update kirchoff stresses
    elast!(mpD,cmParam.Del,isΔFbar,fwrkDeform)
    # plastic corrector
    if plastOn 
        ηmax = plast!(mpD,cmParam,cmType,fwrkDeform) 
    else 
        ηmax = 0 
    end
    # get cauchy stresses
    if fwrkDeform == "finite"
        @threads for p ∈ 1:mpD.nmp
            mpD.σ[:,p] .= mpD.τ[:,p]./mpD.J[p]
        end
    end
    return ηmax
end