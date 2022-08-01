@views function accum!(mn,pn,fen,fin,σ,τ,J,vp,v,mp,ϕ,B,p2n,g,nmp,nn)
    # initialize nodal quantities
    mn .= 0.0
    pn .= 0.0
    fen.= 0.0
    fin.= 0.0
    # accumulate material point contributions
    iD  = zeros(Int64  ,nn)
    buff= zeros(Float64,nn)
    for p in 1:nmp
        # index & buffer
        iD        .= p2n[p,:]
        buff      .= ϕ[p,:].*mp[p]
        # accumulation
        σ[:,p]    .= τ[:,p]./J[p]
        mn[iD  ] .+= buff
        pn[iD,:] .+= repeat(buff,1,2).*repeat(vp[p,:]',nn,1) 
        fen[iD,2].-= buff.*g
        fin[iD,:].+= v[p].*reshape(B[:,:,p]'*σ[:,p],2,nn)' 
    end
end