
@views function deform!(τ,ϵ,ΔJ,J,Jbar,Jn,v,Vn,v0,l,l0,F,un,ϕ,∂ϕx,∂ϕz,p2n,nmp,Del)
    # init mesh quantities to zero
    Jn .= 0.0
    Vn .= 0.0
    # set identity matrix
    ID   = Matrix(I,2,2)
    # action
    @threads for p in 1:nmp
        # compute incremental deformation gradient
        id     = p2n[p,:]
        Δun    = hcat(un[id,1],un[id,2])
        # compute incremental deformation gradient
        ΔF     = ID+vcat(∂ϕx[p,:]'*Δun,∂ϕz[p,:]'*Δun)'
        ΔJ[p]  = det(ΔF)
        Jbar[p]= det(ΔF*reshape(F[p,:],2,2)')
        # accumulation
        Jn[id] .+= ϕ[p,:].*v[p].*(J[p].*ΔJ[p])
        Vn[id] .+= ϕ[p,:].*v[p] 
    end 
    # compute nodal deformed volume
    Vn .= Jn./Vn 
    # compute determinant Jbar 
    @threads for p in 1:nmp
        Jbar[p] = (ϕ[p,:]'*Vn[p2n[p,:]])'
    end
end
#==#
@views function elast!(τ,ϵ,J,v,v0,l,l0,F,un,∂ϕx,∂ϕz,p2n,nmp,Del)
    ID   = Matrix(I,2,2)
    @threads for p in 1:nmp
        # compute incremental deformation gradient
        id     = p2n[p,:]
        Δun    = hcat(un[id,1],un[id,2])
        # compute incremental deformation gradient
        ΔF     = ID+vcat(∂ϕx[p,:]'*Δun,∂ϕz[p,:]'*Δun)'
        Ft     = ΔF*reshape(F[p,:],2,2)'
        F[p,:].= vec(Ft')
        # update material point's volume and domain length
        J[p]   = det(Ft)
        v[p]   = J[p]*v0[p]
        l[p,:].= J[p]^(1/2).*l0[p,:] 
        # compute logarithmic strain tensor
        ϵt     = [    ϵ[1,p] 0.5*ϵ[4,p];
                  0.5*ϵ[4,p]     ϵ[2,p]]
        λ,n    = eigvals(ϵt),eigvecs(ϵt)
        b      = n*diagm(exp.(2*λ))*n'
        bt     = ΔF*b*ΔF'
        λ,n    = eigvals(bt),eigvecs(bt)
        ϵt     = 0.5.*(n*diagm(log.(λ))*n')
        ϵ[:,p].= vcat(ϵt[1,1],ϵt[2,2],0.0,2*ϵt[1,2])
        # krichhoff stress tensor
        τ[:,p].= (Del*ϵ[:,p]) 

    end
end