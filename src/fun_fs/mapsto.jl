@views function accum!(mpD::NamedTuple,meD::NamedTuple,g::Matrix{Float64})
    # initialize nodal quantities
    meD.m   .= 0.0
    meD.p   .= 0.0
    meD.fext.= 0.0
    meD.fint.= 0.0
    # accumulate material point contributions
    iD   = zeros(Int64  ,meD.nn)
    buff = zeros(Float64,meD.nn)
    for p ∈ 1:mpD.nmp
        # index & buffer
        iD             .= mpD.p2n[p,:]
        buff           .= mpD.ϕ∂ϕ[p,:,1].*mpD.m[p]
        # accumulation
        meD.m[iD  ]   .+= buff
        meD.p[iD,:]   .+= repeat(buff,1,meD.nD).*repeat(mpD.v[p,:]',meD.nn,1) 
        meD.fext[iD,:].-= buff.*g
        meD.fint[iD,:].+= mpD.V[p].*reshape(mpD.B[:,:,p]'*mpD.σ[:,p],meD.nD,meD.nn)' 
    end
    return nothing
end
@views function flipDM!(mpD::NamedTuple,meD::NamedTuple,Δt::Float64)
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
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.m[n]>0.0
            mnT        = fill(1.0/meD.m[n],meD.nD)
            meD.u[n,:].= Δt.*meD.p[n,:].*mnT.*meD.bc[n,:]
        end
    end
    # update material point's displacement
    for p ∈ 1:mpD.nmp
        mpD.u[p,:].+= (mpD.ϕ∂ϕ[p,:,1]'*meD.u[mpD.p2n[p,:],:])'
    end
    return nothing
end
@views function mapsto!(mpD::NamedTuple,meD::NamedTuple,g::Matrix{Float64},Δt::Float64,to::String)
    if to == "p->N"
        accum!(mpD,meD,g)
    elseif to == "p<-N"
        flipDM!(mpD,meD,Δt)
    end
    return nothing
end