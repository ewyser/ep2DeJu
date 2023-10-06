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
        meD.pn[iD,:] .+= repeat(buff,1,2).*repeat(mpD.vp[p,:]',meD.nn,1) 
        meD.fen[iD,2].-= buff.*g
        meD.fin[iD,:].+= mpD.v[p].*reshape(mpD.B[:,:,p]'*mpD.σ[:,p],2,meD.nn)' 
    end
end