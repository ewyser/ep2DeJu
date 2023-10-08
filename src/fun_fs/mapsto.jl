@views function accum!(mpD,meD,g)
    # initialize nodal quantities
    meD.m   .= 0.0
    meD.p   .= 0.0
    meD.fext.= 0.0
    meD.fint.= 0.0
    # accumulate material point contributions
    iD  = zeros(Int64  ,meD.nn)
    buff= zeros(Float64,meD.nn)
    for p ∈ 1:mpD.nmp
        # index & buffer
        iD            .= mpD.p2n[p,:]
        buff          .= mpD.ϕ∂ϕ[p,:,1].*mpD.m[p]
        # accumulation
        mpD.σ[:,p]     .= mpD.τ[:,p]./mpD.J[p]
        meD.m[iD  ]   .+= buff
        meD.p[iD,:]   .+= repeat(buff,1,meD.nD).*repeat(mpD.v[p,:]',meD.nn,1) 
        meD.fext[iD,2].-= buff.*g
        meD.fint[iD,:].+= mpD.V[p].*reshape(mpD.B[:,:,p]'*mpD.σ[:,p],meD.nD,meD.nn)' 
    end
    return nothing
end
@views function flipDM!(mpD,meD,Δt)
    # init.
    iD = zeros(Int64,meD.nn)
    # flip update
    for p ∈ 1:mpD.nmp
        # index & buffer
        iD          .= mpD.p2n[p,:]
        # mapping back to mp's
        mpD.v[p,:].+= Δt.*(mpD.ϕ∂ϕ[p,:,1]'*meD.a[iD,:])'
        mpD.x[p,:].+= Δt.*(mpD.ϕ∂ϕ[p,:,1]'*meD.v[iD,:])'
    end
    # initialize for DM + BCs procedure
    meD.p.= 0.0
    meD.u.= 0.0
    # accumulate material point contributions
    for p ∈ 1:mpD.nmp
        # index & buffer
        iD          .= mpD.p2n[p,:]
        # accumulation
        meD.p[iD,:].+= repeat(mpD.ϕ∂ϕ[p,:,1].*mpD.m[p],1,meD.nD).*repeat(mpD.v[p,:]',meD.nn,1) 
    end    
    # solve for nodal incremental displacement
    @threads for n ∈ 1:meD.nno[3]
        if meD.m[n]>0.0
            mnT        = [1.0/meD.m[n];1.0/meD.m[n]]
            meD.u[n,:].= Δt.*meD.p[n,:].*mnT.*meD.bc[n,:]
        end
    end
    # update material point's displacement
    for p ∈ 1:mpD.nmp
        mpD.u[p,:].+= (mpD.ϕ∂ϕ[p,:,1]'*meD.u[mpD.p2n[p,:],:])'
    end
    return nothing
end
@views function mapsto!(mpD,meD,g,Δt,to)
    if to == "p->N"
        accum!(mpD,meD,g)
    elseif to == "p<-N"
        flipDM!(mpD,meD,Δt)
    end
    return nothing
end