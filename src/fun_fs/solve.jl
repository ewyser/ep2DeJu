@views function solve!(meD,Δt)
    # viscous damping
    η   = 0.1
    # initialize
    meD.f .= 0.0
    meD.a .= 0.0
    meD.v .= 0.0
    # solve momentum equation on the mesh
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.m[n]>0.0 
            m          = (1.0/meD.m[n]).*meD.bc[n,:]                   #(2,)
            meD.f[n,:].= meD.fext[n,:].-meD.fint[n,:]                  #(2,)
            meD.D[n,:].= η.*norm(meD.f[n,:]).*sign.(meD.p[n,:].*m)     #(2,)
            meD.f[n,:].= meD.f[n,:].-meD.D[n,:]                        #(2,)
            meD.a[n,:].= meD.f[n,:].*m                                 #(2,)
            meD.v[n,:].= (meD.p[n,:].+Δt.*meD.f[n,:]).*m               #(2,)
        end
    end
    return nothing
end