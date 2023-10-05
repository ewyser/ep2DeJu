
# For volumetric locking, F-bar method is used, see DOI: 10.1002/nag.3599
@views function deform!(ΔJ,m,ΔF,ΔFbar,mn,un,ΔJn,ϕ∂ϕ,p2n,nmp)
    # init mesh quantities to zero
    ΔJn .= 0.0
    # set identity matrix
    ID   = Matrix(I,2,2)
    # action
    @threads for p in 1:nmp
        # get nodal incremental displacement
        iD        = p2n[p,:]
        Δun       = un[iD,:]
        # compute incremental deformation gradient
        ΔF[:,:,p].= ID+vcat(ϕ∂ϕ[p,:,2]'*Δun,ϕ∂ϕ[p,:,3]'*Δun)'
        ΔJ[p]     = det(ΔF[:,:,p])
        # accumulation
        ΔJn[iD]  .+= ϕ∂ϕ[p,:].*m[p].*ΔJ[p]  
    end 
    # compute nodal deformation determinant
    @threads for no in eachindex(ΔJn)
        if mn[no]>0.0 
            ΔJn[no] = ΔJn[no]/mn[no]
        end
    end
    # compute determinant Jbar 
    @threads for p in 1:nmp
        ΔFbar[:,:,p].= ΔF[:,:,p].*(((ϕ∂ϕ[p,:,1]'*ΔJn[p2n[p,:]])/ΔJ[p]).^(1/2))
    end
end
@views function elast!(τ,ϵ,ΔJ,J,v,v0,l,l0,ΔF,ΔFbar,F,nmp,Del)
    @threads for p in 1:nmp
        # update deformation gradient
        F[:,:,p].= ΔF[:,:,p]*F[:,:,p]
        # update material point's volume and domain length
        J[p]     = det(F[:,:,p])
        v[p]     = J[p]*v0[p]
        l[p,:]  .= J[p]^(1/2).*l0[p,:] 
        # compute logarithmic strain tensor
        ϵt       = [    ϵ[1,p] 0.5*ϵ[4,p];
                    0.5*ϵ[4,p]     ϵ[2,p]]
        λ,n      = eigvals(ϵt),eigvecs(ϵt)
        b        = n*diagm(exp.(2*λ))*n'
        bt       = ΔFbar[:,:,p]*b*ΔFbar[:,:,p]'
        λ,n      = eigvals(bt),eigvecs(bt)
        ϵt       = 0.5.*(n*diagm(log.(λ))*n')
        ϵ[:,p]  .= vcat(ϵt[1,1],ϵt[2,2],0.0,2*ϵt[1,2])
        # krichhoff stress tensor
        τ[:,p]  .= (Del*ϵ[:,p]) 
    end
end