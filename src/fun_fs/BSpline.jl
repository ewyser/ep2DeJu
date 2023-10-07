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
    return nothing
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function whichType(xn::Float64,xB::Vector{Float64},Δx::Float64)
    if xn==xB[1] ||  xn==xB[2] 
        type = 1::Int64
    elseif xn<(xB[1]+1.1*Δx) && xn>(xB[1]+0.9*Δx) 
        type = 2::Int64
    elseif xn>(xB[1]+1.5*Δx) && xn<(xB[2]-1.5*Δx) 
        type = 3::Int64
    elseif xn<(xB[2]-0.9*Δx) && xn>(xB[2]-1.1*Δx)
        type = 4::Int64
    else
        type = 0::Int64
    end
    return type::Int64
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function ϕ∇ϕ(ξ::Float64,type::Int64,Δx::Float64)
    ϕ = 0.0
    ∂ϕ= 0.0
    if type==1 
        if -2<=ξ && ξ<=-1 
            ϕ = 1/6     *ξ^3+     ξ^2   +2*ξ    +4/3
            ∂ϕ= 3/(6*Δx)*ξ^2+2/Δx*ξ     +2/Δx
        elseif -1<=ξ && ξ<=-0 
            ϕ = -1/6     *ξ^3           +  ξ    +1
            ∂ϕ= -3/(6*Δx)*ξ^2           +  1/Δx
        elseif  0<=ξ && ξ<= 1 
            ϕ =  1/6     *ξ^3           -  ξ    +1
            ∂ϕ=  3/(6*Δx)*ξ^2           -  1/Δx
        elseif  1<=ξ && ξ<= 2 
            ϕ = -1/6     *ξ^3+     ξ^2  -2*ξ    +4/3
            ∂ϕ= -3/(6*Δx)*ξ^2+2/Δx*ξ    -2/Δx
        end    
    elseif type==2 
        if -1<=ξ && ξ<=0 
            ϕ = -1/3 *ξ^3-     ξ^2    +2/3
            ∂ϕ= -1/Δx*ξ^2-2/Δx*ξ
        elseif 0<=ξ && ξ<=1 
            ϕ =  1/2     *ξ^3-     ξ^2    +2/3
            ∂ϕ=  3/(2*Δx)*ξ^2-2/Δx*ξ
        elseif 1<=ξ && ξ<=2 
            ϕ = -1/6     *ξ^3+     ξ^2-2*ξ+4/3
            ∂ϕ= -3/(6*Δx)*ξ^2+2/Δx*ξ  -2/Δx
        end
    elseif type==3 
        if -2<=ξ && ξ <=-1 
            ϕ =  1/6     *ξ^3+     ξ^2+2*ξ+4/3
            ∂ϕ=  3/(6*Δx)*ξ^2+2/Δx*ξ  +2/Δx
        elseif -1<=ξ && ξ<=0 
            ϕ = -1/2     *ξ^3-     ξ^2    +2/3
            ∂ϕ= -3/(2*Δx)*ξ^2-2/Δx*ξ
        elseif  0<=ξ && ξ<=1
            ϕ =  1/2     *ξ^3-     ξ^2    +2/3
            ∂ϕ=  3/(2*Δx)*ξ^2-2/Δx*ξ
        elseif  1<=ξ && ξ<=2    
            ϕ = -1/6     *ξ^3+     ξ^2-2*ξ+4/3
            ∂ϕ= -3/(6*Δx)*ξ^2+2/Δx*ξ  -2/Δx
        end
    elseif type==4
        if -2<=ξ && ξ<=-1
            ϕ =  1/6     *ξ^3+     ξ^2+2*ξ+4/3
            ∂ϕ=  3/(6*Δx)*ξ^2+2/Δx*ξ  +2/Δx 
        elseif -1<=ξ && ξ<=0
            ϕ = -1/2     *ξ^3-     ξ^2    +2/3
            ∂ϕ= -3/(2*Δx)*ξ^2-2/Δx*ξ      
        elseif 0<=ξ && ξ<=1
            ϕ =  1/3     *ξ^3-     ξ^2    +2/3
            ∂ϕ=  3/(3*Δx)*ξ^2-2/Δx*ξ      
        end
    else
        ϕ = 0.0
        ∂ϕ= 0.0
    end    
    return ϕ,∂ϕ
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
@views function ϕ∂ϕ!(mpD,meD)
    topol!(mpD,meD)
    #preprocessing
    xb = copy(meD.xB[1:2])
    zb = copy(meD.xB[3:4])
    Δx = meD.h[1]
    Δz = meD.h[2]
    #action
    @threads for mp ∈ 1:mpD.nmp
        for nn ∈ 1:meD.nn
            # compute basis functions
            id     = mpD.p2n[mp,nn]
            ξ      = (mpD.xp[mp,1] - meD.xn[id])/Δx 
            type   = whichType(meD.xn[id],xb,Δx)
            ϕx,dϕx = ϕ∇ϕ(ξ,type,Δx)
            η      = (mpD.xp[mp,2] - meD.zn[id])/Δz
            type   = whichType(meD.zn[id],zb,Δz)
            ϕz,dϕz = ϕ∇ϕ(η,type,Δz)
            # convolution of basis function
            mpD.ϕ∂ϕ[mp,nn,1] =  ϕx*  ϕz                                        
            mpD.ϕ∂ϕ[mp,nn,2] = dϕx*  ϕz                                        
            mpD.ϕ∂ϕ[mp,nn,3] =  ϕx* dϕz
        end
        # B-matrix assembly
        mpD.B[1,1:2:end,mp].= mpD.ϕ∂ϕ[mp,:,2]
        mpD.B[2,2:2:end,mp].= mpD.ϕ∂ϕ[mp,:,3]
        mpD.B[4,1:2:end,mp].= mpD.ϕ∂ϕ[mp,:,3]
        mpD.B[4,2:2:end,mp].= mpD.ϕ∂ϕ[mp,:,2]
    end
    return nothing
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------