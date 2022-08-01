@views function elast!(τ,ϵ,J,v,v0,l,l0,F,un,∂ϕx,∂ϕz,p2n,nmp,Del)
    I   = [1.0 0.0;0.0 1.0]
    @threads for p in 1:nmp
        # compute incremental deformation gradient
        id     = p2n[p,:]
        Δun    = hcat(un[id,1],un[id,2])
        # compute incremental deformation gradient
        ΔF     = I+vcat(∂ϕx[p,:]'*Δun,∂ϕz[p,:]'*Δun)'
        Ft     = ΔF*reshape(F[p,:],2,2)'
        F[p,:].= vec(Ft')
        # compute logarithmic strain tensor
        ϵt     = [    ϵ[1,p] 0.5*ϵ[4,p];
                  0.5*ϵ[4,p]     ϵ[2,p]]
        λ      = eigvals(ϵt)
        n      = eigvecs(ϵt)
        b      = n*diagm(exp.(2*λ))*n'
        bt     = ΔF*b*ΔF'
        λ      = eigvals(bt)
        n      = eigvecs(bt)
        ϵt     = 0.5.*(n*diagm(log.(λ))*n')
        ϵ[:,p].= vcat(ϵt[1,1],ϵt[2,2],0.0,2*ϵt[1,2])
        # krichhoff stress tensor
        τ[:,p].= (Del*ϵ[:,p]) 
        # update material point's volume and domain length
        J[p]   = det(Ft)
        v[p]   = J[p]*v0[p]
        l[p,:].= J[p]^(1/2).*l0[p,:] 
    end
end