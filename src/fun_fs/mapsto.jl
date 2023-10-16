@views function accum!(mpD,meD,g)
    # initialize nodal quantities
    meD.m   .= 0.0
    meD.p   .= 0.0
    meD.fext.= 0.0
    meD.fint.= 0.0
    # accumulate material point contributions
    iD,ϕmp  = zeros(Int64,meD.nn),zeros(Float64,meD.nn)
    for p ∈ 1:mpD.nmp
        # index & buffer
        iD .= mpD.p2n[p,:]
        ϕmp.= mpD.ϕ∂ϕ[p,:,1]*mpD.m[p]
        # accumulation
        meD.m[iD].+= ϕmp
        for nD in 1:meD.nD
            meD.p[   iD,nD].+= (ϕmp.*mpD.v[p,nD])
            meD.fext[iD,nD].+= (ϕmp.*g[nD])
            meD.fint[iD,nD].+= (mpD.V[p].*(mpD.B[:,nD:meD.nD:end,p]'*mpD.σ[:,p]))
        end
    end
    return nothing
end
@views function flipDM!(mpD,meD,Δt)
    # flip update
    @threads for p ∈ 1:mpD.nmp
        # mapping back to mp's
        for nD in 1:meD.nD
            mpD.v[p,nD]+= Δt*(mpD.ϕ∂ϕ[p,:,1]'*meD.a[mpD.p2n[p,:],nD])
            mpD.x[p,nD]+= Δt*(mpD.ϕ∂ϕ[p,:,1]'*meD.v[mpD.p2n[p,:],nD])
        end          
    end
    # initialize for DM + BCs procedure
    meD.p.= 0.0
    meD.u.= 0.0
    # accumulate material point contributions
    iD  = zeros(Int64,meD.nn)
    for p ∈ 1:mpD.nmp
        # index & buffer
        iD .= mpD.p2n[p,:]
        # accumulation
        for nD in 1:meD.nD
            meD.p[iD,nD].+= (mpD.ϕ∂ϕ[p,:,1].*mpD.m[p].*mpD.v[p,nD])
        end
    end    
    # solve for nodal incremental displacement
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.m[n]>0.0
            for nD in 1:meD.nD
                meD.u[n,nD]= (Δt*meD.p[n,nD]*(1.0/meD.m[n])*meD.bc[n,nD])
            end    
        end
    end
    # update material point's displacement
    @threads for p ∈ 1:mpD.nmp
        for nD in 1:meD.nD
            mpD.u[p,nD]+= (mpD.ϕ∂ϕ[p,:,1]'*meD.u[mpD.p2n[p,:],nD])
        end        
    end
    return nothing
end
@views function mapsto!(mpD,meD,g,Δt,whereto)
    if whereto == "p->N"
        accum!(mpD,meD,g)
    elseif whereto == "p<-N"
        flipDM!(mpD,meD,Δt)
    end
    return nothing
end