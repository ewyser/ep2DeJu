@views function accum!(mpD,meD,g)
    # initialize nodal quantities
    meD.mn .= 0.0
    meD.pn .= 0.0
    meD.fen.= 0.0
    meD.fin.= 0.0
    # accumulate material point contributions
    iD  = zeros(Int64  ,meD.nn)
    buff= zeros(Float64,meD.nn)
    for p ∈ 1:mpD.nmp
        # index & buffer
        iD            .= mpD.p2n[p,:]
        buff          .= mpD.ϕ∂ϕ[p,:,1].*mpD.mp[p]
        # accumulation
        mpD.σ[:,p]    .= mpD.τ[:,p]./mpD.J[p]
        meD.mn[iD  ] .+= buff
        meD.pn[iD,:] .+= repeat(buff,1,meD.dof).*repeat(mpD.vp[p,:]',meD.nn,1) 
        meD.fen[iD,2].-= buff.*g
        meD.fin[iD,:].+= mpD.v[p].*reshape(mpD.B[:,:,p]'*mpD.σ[:,p],meD.dof,meD.nn)' 
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
        mpD.vp[p,:].+= Δt.*(mpD.ϕ∂ϕ[p,:,1]'*meD.an[iD,:])'
        mpD.xp[p,:].+= Δt.*(mpD.ϕ∂ϕ[p,:,1]'*meD.vn[iD,:])'
    end
    # initialize for DM + BCs procedure
    meD.pn.= 0.0
    meD.un.= 0.0
    # accumulate material point contributions
    for p ∈ 1:mpD.nmp
        # index & buffer
        iD           .= mpD.p2n[p,:]
        # accumulation
        meD.pn[iD,:].+= repeat(mpD.ϕ∂ϕ[p,:,1].*mpD.mp[p],1,meD.dof).*repeat(mpD.vp[p,:]',meD.nn,1) 
    end    
    # solve for nodal incremental displacement
    @threads for n ∈ 1:meD.nno[3]
        if meD.mn[n]>0.0
            mnT         = [1.0/meD.mn[n];1.0/meD.mn[n]]
            meD.un[n,:].= Δt.*meD.pn[n,:].*mnT.*meD.bc[n,:]
        end
    end
    # update material point's displacement
    for p ∈ 1:mpD.nmp
        mpD.up[p,:].+= (mpD.ϕ∂ϕ[p,:,1]'*meD.un[mpD.p2n[p,:],:])'
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