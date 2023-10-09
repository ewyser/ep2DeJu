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
            m           = fill(1.0/meD.m[n],meD.nD)                     #(2,)
            f           = meD.fext[n,:].-meD.fint[n,:]                  #(2,)
            D           = η.*sqrt(f[1]^2+f[2]^2).*sign.(meD.p[n,:].*m)  #(2,)
            f           = f.-D                                          #(2,)
            meD.a[n,:] .= f.*m.*meD.bc[n,:]                             #(2,)
            meD.v[n,:] .= (meD.p[n,:].+Δt.*f).*m.*meD.bc[n,:]           #(2,)
        end
    end
    return nothing
end