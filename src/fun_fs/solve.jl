@views function solve!(meD,bc,Δt)
    D   = 0.1
    # initialize
    meD.fn .= 0.0
    meD.an .= 0.0
    meD.vn .= 0.0
    # solve momentum equation on the mesh
    @threads for n ∈ 1:meD.nno[3]
        if meD.mn[n]>0.0 
            mnT          = [1.0/meD.mn[n];1.0/meD.mn[n]] #(2,)
            fnT          = meD.fen[n,:].-meD.fin[n,:]    #(2,)
            vnT          = meD.pn[n,:] .*mnT         #(2,)
            η            = sqrt(fnT[1]^2+fnT[2]^2)      #()
            fnT          = fnT .- D.*η.*sign.(vnT)      #(2,)
            meD.an[n,:] .= fnT.*mnT.*[bc.x[n];bc.z[n]]
            meD.vn[n,:] .= (meD.pn[n,:].+Δt.*fnT).*mnT.*[bc.x[n];bc.z[n]]
        end
    end
    return nothing
end