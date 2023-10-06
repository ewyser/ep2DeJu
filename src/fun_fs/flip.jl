@views function flipDM!(mpD,meD,bcx,bcz,Δt)
    # init.
    iD = zeros(Int64,meD.nn)
    # flip update
    @threads for p in 1:mpD.nmp
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
    @threads for p in 1:mpD.nmp
        # index & buffer
        iD           .= mpD.p2n[p,:]
        # accumulation
        meD.pn[iD,:].+= repeat(mpD.ϕ∂ϕ[p,:,1].*mpD.mp[p],1,2).*repeat(mpD.vp[p,:]',meD.nn,1) 
    end    
    # solve for nodal incremental displacement
    @threads for n in 1:meD.nno[3]
        if(meD.mn[n]>0.0)
            mnT         = [1.0/meD.mn[n] 1.0/meD.mn[n]]
            meD.un[n,:].= reshape(Δt.*meD.pn[n,:]'.*mnT.*[bcx[n] bcz[n]],2)
        end
    end
    # update material point's displacement
    @threads for p in 1:mpD.nmp
        mpD.up[p,:].+= (mpD.ϕ∂ϕ[p,:,1]'*meD.un[mpD.p2n[p,:],:])'
    end
end