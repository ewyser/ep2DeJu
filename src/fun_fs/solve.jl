@views function solve!(meD,Δt)
    D   = 0.1
    # initialize
    meD.f .= 0.0
    meD.a .= 0.0
    meD.v .= 0.0
    # solve momentum equation on the mesh
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.m[n]>0.0 
            mnT          = fill(1.0/meD.m[n],meD.nD) #(2,)
            fnT          = meD.fext[n,:].-meD.fint[n,:]    #(2,)
            vnT          = meD.p[n,:] .*mnT         #(2,)
            η            = sqrt(fnT[1]^2+fnT[2]^2)      #()
            fnT          = fnT .- D.*η.*sign.(vnT)      #(2,)
            meD.a[n,:] .= fnT.*mnT.*meD.bc[n,:]
            meD.v[n,:] .= (meD.p[n,:].+Δt.*fnT).*mnT.*meD.bc[n,:]
        end
    end
    return nothing
end