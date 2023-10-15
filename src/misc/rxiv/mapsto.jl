@views function accum!(mpD::NamedTuple,meD::NamedTuple,g::Matrix{Float64})
    # initialize nodal quantities
    meD.m   .= 0.0
    meD.p   .= 0.0
    meD.fext.= 0.0
    meD.fint.= 0.0
    # accumulate material point contributions
    for n ∈ 1:meD.nn
        for p ∈ 1:mpD.nmp
            # index & buffer
            iD             = mpD.p2n[p,n]
            buff           = mpD.ϕ∂ϕ[p,n,1]*mpD.m[p]
            # accumulation
            meD.m[iD  ]   += buff
            for dim ∈ 1:meD.nD
                meD.p[iD,dim]   += buff*mpD.v[p,dim]
                meD.p[iD,dim]   += buff*mpD.v[p,dim]
                meD.fext[iD,dim]+= g[1]
                meD.fint[iD,dim]+= mpD.V[p].*reshape(mpD.B[:,:,p]'*mpD.σ[:,p],meD.nD,meD.nn)' 
                meD.fint[iD,dim]+= mpD.V[p]*(mpD.∂ϕx[p,n,meD.nD]*mpD.σ[1,p]+mpD.∂ϕz[p,n,meD.nD+1]*mpD.σ[4,p])
            end
        end
    end
    return nothing
end
@views function flipDM!(mpD::NamedTuple,meD::NamedTuple,Δt::Float64)
    # flip update
    for p ∈ 1:mpD.nmp
        dvpx = dvpy = dxp = dyp = 0.0
        for n ∈ 1:nn
            iD    = mpD.p2n[p,n]
            dvpx += mpD.ϕ∂ϕ[p,n,1]*meD.a[iD,1]
            dvpy += mpD.ϕ∂ϕ[p,n,1]*meD.a[iD,2]
            dxp  += mpD.ϕ∂ϕ[p,n,1]*meD.v[iD,1]
            dyp  += mpD.ϕ∂ϕ[p,n,1]*meD.v[iD,2]
        end
        mpD.v[p,1] += Δt*dvpx
        mpD.v[p,2] += Δt*dvpy
        mpD.x[p,1] += Δt*dxp
        mpD.x[p,2] += Δt*dyp
    end
    # initialize for DM + BCs procedure
    meD.p.= 0.0
    meD.u.= 0.0
    # accumulate material point contributions
    for n ∈ 1:meD.nn
        for p ∈ 1:mpD.nmp
            # index & buffer
            iD           = mpD.p2n[p,n]
            buff         = mpD.ϕ∂ϕ[p,n,1]*mpD.m[p]
            # accumulation
            meD.p[iD,:] += buff*mpD.v[p,:]
            meD.p[iD,:] += buff*mpD.v[p,:]
        end
    end 
    # solve for nodal incremental displacement
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.m[n]>0.0
            m          = 1.0/meD.m[n]
            meD.u[n,1] = (Δt*meD.p[n,1]*m)*meD.bc[n,1]
            meD.u[n,2] = (Δt*meD.p[n,2]*m)*meD.bc[n,2]
        end
    end
    # update material point's displacement
    for p in 1:mpD.nmp
        dupx = dupz = 0.0
        for n in 1:meD.nn
            iD    = mpD.p2n[p,n]
            dupx += mpD.ϕ∂ϕ[p,n,1]*meD.u[iD,1]
            dupz += mpD.ϕ∂ϕ[p,n,1]*meD.u[iD,2]
        end
        mpD.u[p,1] += Δt*dupx
        mpD.u[p,2] += Δt*dupy
    end
    return nothing
end
@views function mapsto!(mpD::NamedTuple,meD::NamedTuple,g::Matrix{Float64},Δt::Float64,whereto::String)
    if whereto == "p->N"
        accum!(mpD,meD,g)
    elseif whereto == "p<-N"
        flipDM!(mpD,meD,Δt)
    end
    return nothing
end