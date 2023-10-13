@views function deform!(mpD,meD)
    # init mesh quantities to zero
    meD.ΔJ .= 0.0
    # set identity matrix
    ID   = Matrix(1.0I,meD.nD,meD.nD)
    # action
    for p ∈ 1:mpD.nmp
        # get nodal incremental displacement
        iD            = mpD.p2n[p,:]
        # compute incremental deformation gradient
        mpD.ΔF[:,:,p].= ID+(permutedims(mpD.ϕ∂ϕ[p,:,2:end],(2,1))*meD.u[iD,:])'
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
    for p ∈ 1:mpD.nmp
        mpD.ΔFbar[:,:,p].= mpD.ΔF[:,:,p].*((dot(mpD.ϕ∂ϕ[p,:,1],meD.ΔJ[mpD.p2n[p,:]])/mpD.ΔJ[p]).^(1.0/meD.nD))
    end
    return nothing
end
@views function mutate(ϵ,type)
    if type == "tensor"
        ϵmut = [     ϵ[1] 0.5*ϵ[4];
                 0.5*ϵ[4]     ϵ[2]]
    elseif type == "voigt"
        ϵmut = vcat(ϵ[1,1],ϵ[2,2],0.0,2*ϵ[1,2])
    end
    return ϵmut
end
# For volumetric locking, F-bar method is used, see DOI: 10.1002/nag.3599
@views function elast!(mpD,Del,isΔFbar,fwrkDeform)
    if fwrkDeform == "finite"
        @threads for p ∈ 1:mpD.nmp
            # compute logarithmic strain tensor
            λ,n          = eigen(mpD.ϵ[:,:,p])
            mpD.b[:,:,p].= n*diagm(exp.(2.0*λ))*n'
            if isΔFbar
                mpD.bT[:,:,p].= mul!(mpD.bT[:,:,p],mpD.ΔFbar[:,:,p],mpD.b[:,:,p])*mpD.ΔFbar[:,:,p]'
            else
                mpD.bT[:,:,p].= mul!(mpD.bT[:,:,p],mpD.ΔF[:,:,p],mpD.b[:,:,p])*mpD.ΔF[:,:,p]'
            end
            λ,n          = eigen(mpD.bT[:,:,p])
            mpD.ϵ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
            # krichhoff stress tensor
            mul!(mpD.τ[:,p],Del,mutate(mpD.ϵ[:,:,p],"voigt"))
        end
    elseif fwrkDeform == "infinitesimal"
        ID = Matrix(1.0I,size(mpD.ΔF,1),size(mpD.ΔF,2))
        @threads for p ∈ 1:mpD.nmp
            σ0 = zeros(4)
            # calculate elastic strains
            if isΔFbar
                mpD.ϵ[:,:,p].= 0.5.*(mpD.ΔFbar[:,:,p]+mpD.ΔFbar[:,:,p]').-ID
                mpD.ω[p]     = 0.5*(mpD.ΔFbar[1,2,p]-mpD.ΔFbar[2,1,p])
            else
                mpD.ϵ[:,:,p].= 0.5.*(mpD.ΔF[:,:,p]+mpD.ΔF[:,:,p]').-ID
                mpD.ω[p]     = 0.5*(mpD.ΔF[1,2,p]-mpD.ΔF[2,1,p])
            end
            # update cauchy stress tensor
            σ0         .= [2.0*mpD.σ[4,p]*mpD.ω[p],-2.0*mpD.σ[4,p]*mpD.ω[p],0.0,(mpD.σ[2,p]-mpD.σ[1,p])*mpD.ω[p]]
            mpD.σ[:,p].+= Del*mutate(mpD.ϵ[:,:,p],"voigt").+σ0
            
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