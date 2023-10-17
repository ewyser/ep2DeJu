@views function deform!(mpD,meD)
    # init mesh quantities to zero
    meD.ΔJn .= 0.0
    # calculate cst.
    dim  = 1.0/meD.nD
    # action
    iD   = zeros(Int64,meD.nn)
    @simd for p ∈ 1:mpD.nmp
        # get nodal incremental displacement
        iD           .= mpD.p2n[p,:]
        # compute incremental deformation gradient
        mpD.ΔF[:,:,p].= mpD.I+(permutedims(mpD.ϕ∂ϕ[p,:,2:end],(2,1))*meD.Δun[iD,:])'
        mpD.ΔJ[p]     = det(mpD.ΔF[:,:,p])
        # update deformation gradient
        mpD.F[:,:,p] .= mpD.ΔF[:,:,p]*mpD.F[:,:,p]
        # update material point's volume and domain length
        mpD.J[p]      = det(mpD.F[:,:,p])
        mpD.V[p]      = mpD.J[p]*mpD.V0[p]
        mpD.l[p,:]   .= mpD.J[p]^(dim).*mpD.l0[p,:] 
        # accumulation
        meD.ΔJn[iD] .+= mpD.ϕ∂ϕ[p,:].*mpD.m[p].*mpD.ΔJ[p]  
    end 
    # compute nodal deformation determinant
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.mn[n]>0.0 
            meD.ΔJn[n] = meD.ΔJn[n]/meD.mn[n]
        end
    end
    # compute determinant Jbar 
    @threads for p ∈ 1:mpD.nmp
        mpD.ΔFbar[:,:,p].= mpD.ΔF[:,:,p].*((dot(mpD.ϕ∂ϕ[p,:,1],meD.ΔJn[mpD.p2n[p,:]])/mpD.ΔJ[p]).^(dim))
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
@views function finite!(mpD,Del,isΔFbar)
    if isΔFbar
        @threads for p ∈ 1:mpD.nmp
            # update left cauchy-green tensor
            mpD.b[:,:,p].= mpD.ΔFbar[:,:,p]*mpD.b[:,:,p]*mpD.ΔFbar[:,:,p]'
            # compute logarithmic strain tensor
            λ,n          = eigen(mpD.b[:,:,p])
            mpD.ϵ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
            # krichhoff stress tensor
            mul!(mpD.τ[:,p],Del,mutate(mpD.ϵ[:,:,p],"voigt"))
        end
    else
        @threads for p ∈ 1:mpD.nmp
            # update left cauchy-green tensor
            mpD.b[:,:,p].= mpD.ΔF[:,:,p]*mpD.b[:,:,p]*mpD.ΔF[:,:,p]'
            # compute logarithmic strain tensor
            λ,n          = eigen(mpD.b[:,:,p])
            mpD.ϵ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
            # krichhoff stress tensor
            mul!(mpD.τ[:,p],Del,mutate(mpD.ϵ[:,:,p],"voigt"))
        end
    end
    return nothing
end
@views function inifinitesimal!(mpD,Del,isΔFbar)
    if isΔFbar
        @threads for p ∈ 1:mpD.nmp
            # calculate elastic strains
            mpD.ϵ[:,:,p].= 0.5.*(mpD.ΔFbar[:,:,p]+mpD.ΔFbar[:,:,p]')-mpD.I
            mpD.ω[p]     = 0.5.*(mpD.ΔFbar[1,2,p]-mpD.ΔFbar[2,1,p])
            # update cauchy stress tensor
            mpD.σR[:,p].= [2.0*mpD.σ[4,p]*mpD.ω[p],-2.0*mpD.σ[4,p]*mpD.ω[p],0.0,(mpD.σ[2,p]-mpD.σ[1,p])*mpD.ω[p]]
            mpD.σ[:,p].+= Del*mutate(mpD.ϵ[:,:,p],"voigt").+mpD.σR[:,p]
        end   
    else
        @threads for p ∈ 1:mpD.nmp
            # calculate elastic strains
            mpD.ϵ[:,:,p].= 0.5.*(mpD.ΔF[:,:,p]+mpD.ΔF[:,:,p]')-mpD.I
            mpD.ω[p]     = 0.5.*(mpD.ΔF[1,2,p]-mpD.ΔF[2,1,p])
            # update cauchy stress tensor
            mpD.σR[:,p].= [2.0*mpD.σ[4,p]*mpD.ω[p],-2.0*mpD.σ[4,p]*mpD.ω[p],0.0,(mpD.σ[2,p]-mpD.σ[1,p])*mpD.ω[p]]
            mpD.σ[:,p].+= Del*mutate(mpD.ϵ[:,:,p],"voigt").+mpD.σR[:,p]
        end   
    end
    return nothing
end
# For volumetric locking, F-bar method is used, see DOI: 10.1002/nag.3599
@views function elast!(mpD,Del,isΔFbar,fwrkDeform)
    if fwrkDeform == "finite"
        finite!(mpD,Del,isΔFbar) 
    elseif fwrkDeform == "infinitesimal"
        inifinitesimal!(mpD,Del,isΔFbar)
    end
    return nothing
end
@views function elastoplast!(mpD,meD,cmParam,cmType,isΔFbar,fwrkDeform,plastOn)
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
    return ηmax::Int64
end