@views function deform!(mpD,meD)
    # init mesh quantities to zero
    meD.ΔJ .= 0.0
    # set identity matrix
    ID   = Matrix(I,2,2)
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
        mpD.l[p,:]   .= mpD.J[p]^(1/2).*mpD.l0[p,:] 
        # accumulation
        meD.ΔJ[iD]  .+= mpD.ϕ∂ϕ[p,:].*mpD.m[p].*mpD.ΔJ[p]  
    end 
    # compute nodal deformation determinant
    @threads for no ∈ 1:meD.nno[3]
        if meD.m[no]>0.0 
            meD.ΔJ[no] = meD.ΔJ[no]/meD.m[no]
        end
    end
    # compute determinant Jbar 
    @threads for p ∈ 1:mpD.nmp
        mpD.ΔFbar[:,:,p].= mpD.ΔF[:,:,p].*(((mpD.ϕ∂ϕ[p,:,1]'*meD.ΔJ[mpD.p2n[p,:]])/mpD.ΔJ[p]).^(1/2))
    end
    return nothing
end
# For volumetric locking, F-bar method is used, see DOI: 10.1002/nag.3599
@views function elast!(mpD,Del,isΔFbar)
    @threads for p ∈ 1:mpD.nmp
        # compute logarithmic strain tensor
        ϵ          = [    mpD.ϵ[1,p] 0.5*mpD.ϵ[4,p];
                      0.5*mpD.ϵ[4,p]     mpD.ϵ[2,p]]
        λ,n        = eigvals(ϵ),eigvecs(ϵ)
        b          = n*diagm(exp.(2*λ))*n'
        if isΔFbar
            bt = mpD.ΔFbar[:,:,p]*b*mpD.ΔFbar[:,:,p]'
        else
            bt = mpD.ΔF[:,:,p]*b*mpD.ΔF[:,:,p]'
        end
        λ,n        = eigvals(bt),eigvecs(bt)
        ϵt         = 0.5.*(n*diagm(log.(λ))*n')
        mpD.ϵ[:,p].= vcat(ϵt[1,1],ϵt[2,2],0.0,2*ϵt[1,2])
        # krichhoff stress tensor
        mpD.τ[:,p].= (Del*mpD.ϵ[:,p]) 
    end
    return nothing
end
@views function elastoplast!(mpD,meD,K,Del,Hp,cr,isΔFbar,cmType,plastOn)
    # get def. & logarithmic strains
    deform!(mpD,meD)
    # update kirchoff stresses
    elast!(mpD,Del,isΔFbar)
    # plastic corrector
    if plastOn 
        ηmax = plast!(mpD,K,Del,Hp,cr,cmType) 
    else 
        ηmax=0 
    end
    # get cauchy stresses
    @threads for p ∈ 1:mpD.nmp
        mpD.σ[:,p] .= mpD.τ[:,p]./mpD.J[p]
    end
    return ηmax
end