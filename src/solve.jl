@views function solve!(fn,an,pn,vn,mn,fen,fin,bcx,bcz,nno,Δt)
    D   = 0.1
    # initialize
    fn .= 0.0
    an .= 0.0
    vn .= 0.0
    # solve momentum equation on the mesh
    @threads for n in 1:nno[3]
        if(mn[n]>0.0)
            mnT      = [1.0/mn[n];1.0/mn[n]] #(2,)
            fnT      = fen[n,:].-fin[n,:]    #(2,)
            vnT      = pn[n,:] .*mnT         #(2,)

            η        = sqrt(fnT[1]^2+fnT[2]^2)      #()
            fnT      = fnT .- D.*η.*sign.(vnT)      #(2,)
            
            an[n,:] .= fnT.*mnT.*[bcx[n];bcz[n]]
            vn[n,:] .= (pn[n,:].+Δt.*fnT).*mnT.*[bcx[n];bcz[n]]
        end
    end
end