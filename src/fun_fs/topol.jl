@views function topol!(mpD,meD)
    xmin = minimum(meD.xn)
    zmin = minimum(meD.zn)
    Δx::Float64 = 1.0/meD.h[1]
    Δz::Float64 = 1.0/meD.h[2]
    nez::Int64  = meD.nel[2]
    id::Int64   = 0
    @threads for p ∈ 1:mpD.nmp
        id = (floor(Int64,(mpD.xp[p,2]-zmin)*Δz)+1::Int64)+(nez)*floor(Int64,(mpD.xp[p,1]-xmin)*Δx)
        for n ∈ 1:meD.nn
            mpD.p2n[p,n] = meD.e2n[id,n]
        end
        mpD.p2e[p] = id
    end
end
#@threads