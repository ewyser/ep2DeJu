@views function topol!(p2e,p2n,e2n,xn,zn,xp,h,nel,nmp,nn)
    xmin = minimum(xn)
    zmin = minimum(zn)
    Δx::Float64 = 1.0/h[1]
    Δz::Float64 = 1.0/h[2]
    nez::Int64  = nel[2]
    id::Int64   = 0
    for p ∈ 1:nmp
        id = (floor(Int64,(xp[p,2]-zmin)*Δz)+1::Int64)+(nez)*floor(Int64,(xp[p,1]-xmin)*Δx)
        for n ∈ 1:nn
            p2n[p,n] = e2n[id,n]
        end
        p2e[p] = id
    end
end