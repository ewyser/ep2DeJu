@views function DMBC!(up,pn,un,mn,ϕ,vp,mp,p2n,bcx,bcz,nmp,nn,nno,Δt)
    # initialize
    pn.= 0.0
    un.= 0.0
    # accumulate material point contributions
    iD = zeros(Int64,nn)
    for p in 1:nmp
        # index & buffer
        iD       .= p2n[p,:]
        # accumulation
        pn[iD,:].+= repeat(ϕ[p,:].*mp[p],1,2).*repeat(vp[p,:]',nn,1) 
    end    
    # solve for nodal incremental displacement
    @threads for n in 1:nno[3]
        if(mn[n]>0.0)
            mnT     = [1.0/mn[n] 1.0/mn[n]]
            un[n,:].= reshape(Δt.*pn[n,:]'.*mnT.*[bcx[n] bcz[n]],2)
        end
    end
    # update material point's displacement
    @threads for p in 1:nmp
        up[p,:].+= (ϕ[p,:]'*un[p2n[p,:],:])'
    end
end