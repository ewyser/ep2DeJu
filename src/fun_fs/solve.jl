@views function solve!(meD,Δt)
    # viscous damping
    η       = 0.1
    # initialize
    meD.vn .= 0.0
    # solve momentum equation on the mesh
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.mn[n]>0.0 
            m           = (1.0/meD.mn[n]).*meD.bc[n,:]                   #(2,)
            meD.Dn[n,:].= η.*norm(meD.oobf[n,:]).*sign.(meD.pn[n,:].*m)  #(2,)
            meD.fn[n,:].= meD.oobf[n,:].-meD.Dn[n,:]                     #(2,)
            meD.vn[n,:].= (meD.pn[n,:].+Δt.*meD.fn[n,:]).*m              #(2,)
        end
    end
    return nothing
end