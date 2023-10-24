@views function solve!(meD,Δt)
    # viscous damping
    η   = 0.1
    # initialize
    meD.vn .= 0.0
    meD.Δun.= 0.0
    # solve momentum equation on the mesh
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.mn[n]>0.0 
            m             = (1.0/meD.mn[n]).*meD.bc[n,:]                   #(2,)
            meD.vn[n,:]  .= meD.pn[n,:].*m 
            meD.vn[n,:] .+= (Δt.*meD.oobf[n,:]).*m            #(2,)
            meD.Δun[n,:] .= Δt.*meD.vn[n,:]
        end
    end
    return nothing
end