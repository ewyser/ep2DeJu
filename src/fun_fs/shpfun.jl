@views function topol!(mpD,meD)
    xmin,zmin = meD.minC[1],meD.minC[2]
    Δx,Δz     = 1.0/meD.h[1],1.0/meD.h[2]
    nez       = meD.nel[2]
    @threads for p ∈ 1:mpD.nmp
        mpD.p2e[p]   = (floor(Int64,(mpD.x[p,2]-zmin)*Δz)+1)+(nez)*floor(Int64,(mpD.x[p,1]-xmin)*Δx)
        mpD.p2n[:,p].= meD.e2n[:,mpD.p2e[p]]
    end
    return nothing
end
function NdN(δx::Float64,h::Float64,lp::Float64)                                                         
    if abs(δx) < lp                       
        ϕ  = 1.0-((4.0*δx^2+(2.0*lp)^2)/(8.0*h*lp))                                   
        ∂ϕ = -((8.0*δx)/(8.0*h*lp))                                     
    elseif (abs(δx)>=   lp ) && (abs(δx)<=(h-lp))
        ϕ  = 1.0-(abs(δx)/h)                                                       
        ∂ϕ = sign(δx)*(-1.0/h)                                                   
    elseif (abs(δx)>=(h-lp)) && (abs(δx)< (h+lp))
        ϕ  = ((h+lp-abs(δx))^2)/(4.0*h*lp)                                       
        ∂ϕ = -sign(δx)*(h+lp-abs(δx))/(2.0*h*lp)
    else
        ϕ  = 0.0                                                                 
        ∂ϕ = 0.0                                  
    end
    return ϕ,∂ϕ    
end
@views function whichType(xn,xB,Δx)
    type = -404
    if xn==xB[1] || xn==xB[2] 
        type = 1
    elseif (xB[1]+0.9*Δx)<xn<(xB[1]+1.1*Δx)
        type = 2
    elseif (xB[1]+1.5*Δx)<xn<(xB[2]-1.5*Δx) 
        type = 3
    elseif (xB[2]-1.1*Δx)<xn<(xB[2]-0.9*Δx)
        type = 4
    end
    return type
end
function ϕ∇ϕ(ξ,type,Δx)
    ϕ,∂ϕ = 0.0,0.0
    if type == 1 
        if -2.0<=ξ<=-1.0 
            ϕ = 1.0/6.0     *ξ^3+     ξ^2   +2.0*ξ    +4.0/3.0
            ∂ϕ= 3.0/6.0     *ξ^2+2.0 *ξ     +2.0
        elseif -1.0<=ξ<=0.0 
            ϕ = -1.0/6.0     *ξ^3           +    ξ    +1.0
            ∂ϕ= -3.0/6.0     *ξ^2           +  1.0
        elseif  0.0<=ξ<= 1.0 
            ϕ =  1.0/6.0     *ξ^3           -    ξ    +1.0
            ∂ϕ=  3.0/6.0     *ξ^2           -  1.0
        elseif  1.0<=ξ<= 2.0 
            ϕ = -1.0/6.0     *ξ^3+     ξ^2  -2.0*ξ    +4.0/3.0
            ∂ϕ= -3.0/6.0     *ξ^2+2.0 *ξ    -2.0
        end    
    elseif type == 2 
        if -1.0<=ξ<=0.0 
            ϕ = -1.0/3.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ= -3.0/3.0     *ξ^2-2.0 *ξ
        elseif 0.0<=ξ<=1.0 
            ϕ =  1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ=  3.0/2.0     *ξ^2-2.0 *ξ
        elseif 1.0<=ξ<=2.0 
            ϕ = -1.0/6.0     *ξ^3+     ξ^2-2.0*ξ+4.0/3.0
            ∂ϕ= -3.0/6.0     *ξ^2+2.0 *ξ  -2.0
        end
    elseif type == 3 
        if -2.0<=ξ<=-1.0 
            ϕ =  1.0/6.0     *ξ^3+     ξ^2+2.0*ξ+4.0/3.0
            ∂ϕ=  3.0/6.0     *ξ^2+2.0 *ξ  +2.0
        elseif -1.0<=ξ<=0.0 
            ϕ = -1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ= -3.0/2.0     *ξ^2-2.0 *ξ
        elseif  0.0<=ξ<=1.0
            ϕ =  1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ=  3.0/2.0     *ξ^2-2.0 *ξ
        elseif  1.0<=ξ<=2.0    
            ϕ = -1.0/6.0     *ξ^3+     ξ^2-2.0*ξ+4.0/3.0
            ∂ϕ= -3.0/6.0     *ξ^2+2.0 *ξ  -2.0
        end
    elseif type == 4
        if -2.0<=ξ<=-1.0
            ϕ =  1.0/6.0     *ξ^3+     ξ^2+2.0*ξ+4.0/3.0
            ∂ϕ=  3.0/6.0     *ξ^2+2.0 *ξ  +2.0 
        elseif -1.0<=ξ<=0.0
            ϕ = -1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ= -3.0/2.0     *ξ^2-2.0 *ξ      
        elseif 0.0<=ξ<=1.0
            ϕ =  1.0/3.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ=  3.0/3.0     *ξ^2-2.0 *ξ      
        end
    end    
    ∂ϕ/=Δx
    return ϕ,∂ϕ
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
@views function ϕ∂ϕ!(mpD,meD,ϕ∂ϕType)
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    topol!(mpD,meD)
    # calculate shape functions
    if ϕ∂ϕType == :bsmpm
        @threads for mp ∈ 1:mpD.nmp
            @simd for nn ∈ 1:meD.nn
                # compute basis functions
                id     = mpD.p2n[nn,mp]
                ξ      = (mpD.x[mp,1]-meD.xn[id,1])/meD.h[1] 
                type   = whichType(meD.xn[id,1],meD.xB[1:2],meD.h[1])
                ϕx,dϕx = ϕ∇ϕ(ξ,type,meD.h[1])
                η      = (mpD.x[mp,2]-meD.xn[id,2])/meD.h[2]
                type   = whichType(meD.xn[id,2],meD.xB[3:4],meD.h[2])
                ϕz,dϕz = ϕ∇ϕ(η,type,meD.h[2])
                # convolution of basis function
                mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕz        
            end
            # B-matrix assembly
            mpD.B[1:meD.nD:end,1,mp].= mpD.ϕ∂ϕ[:,mp,2]
            mpD.B[2:meD.nD:end,2,mp].= mpD.ϕ∂ϕ[:,mp,3]
            mpD.B[1:meD.nD:end,4,mp].= mpD.ϕ∂ϕ[:,mp,3]
            mpD.B[2:meD.nD:end,4,mp].= mpD.ϕ∂ϕ[:,mp,2]
        end
    elseif ϕ∂ϕType == :gimpm
        @threads for mp in 1:mpD.nmp
            @simd for nn in 1:meD.nn
                # compute basis functions
                id     = mpD.p2n[nn,mp]
                ξ      = (mpD.x[mp,1]-meD.xn[id,1])
                η      = (mpD.x[mp,2]-meD.xn[id,2])
                ϕx,dϕx = NdN(ξ,meD.h[1],mpD.l[mp,1])
                ϕz,dϕz = NdN(η,meD.h[2],mpD.l[mp,2])
                # convolution of basis function
                mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕz
            end
            # B-matrix assembly
            mpD.B[1:meD.nD:end,1,mp].= mpD.ϕ∂ϕ[:,mp,2]
            mpD.B[2:meD.nD:end,2,mp].= mpD.ϕ∂ϕ[:,mp,3]
            mpD.B[1:meD.nD:end,4,mp].= mpD.ϕ∂ϕ[:,mp,3]
            mpD.B[2:meD.nD:end,4,mp].= mpD.ϕ∂ϕ[:,mp,2]
        end
    elseif ϕ∂ϕType == :smpm
        @threads for mp in 1:mpD.nmp
            @simd for nn in 1:meD.nn
                # compute basis functions
                id     = mpD.p2n[nn,mp]
                ξ      = (mpD.x[mp,1]-meD.xn[id,1])
                η      = (mpD.x[mp,2]-meD.xn[id,2])
                ϕx,dϕx = NdN(ξ,meD.h[1],0.0)
                ϕz,dϕz = NdN(η,meD.h[2],0.0)
                # convolution of basis function
                mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕz
            end
            # B-matrix assembly
            mpD.B[1:meD.nD:end,1,mp].= mpD.ϕ∂ϕ[:,mp,2]
            mpD.B[2:meD.nD:end,2,mp].= mpD.ϕ∂ϕ[:,mp,3]
            mpD.B[1:meD.nD:end,4,mp].= mpD.ϕ∂ϕ[:,mp,3]
            mpD.B[2:meD.nD:end,4,mp].= mpD.ϕ∂ϕ[:,mp,2]
        end        
    end
    return nothing
end



























#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
@views function topol3D!(mpD,meD)
    xmin,ymin,zmin = meD.minC[1],meD.minC[2],meD.minC[3]
    Δx,Δy,Δz       = 1.0/meD.h[1],1.0/meD.h[2],1.0/meD.h[3]
    nex,ney,nez    = meD.nel[1],meD.nel[2],meD.nel[3]
    @threads for p ∈ 1:mpD.nmp
        mpD.p2e[p  ] = (floor(Int64,(mpD.x[p,3]-zmin)*Δz)+1)+(nez)*floor(Int64,(mpD.x[p,1]-xmin)*Δx)+(nez*nex)*floor(Int64,(mpD.x[p,2]-ymin)*Δy)
        mpD.p2n[:,p].= meD.e2n[:,mpD.p2e[p]]
    end
    return nothing
end
@views function ϕ∂ϕ3D!(mpD,meD,ϕ∂ϕType)
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    topol3D!(mpD,meD)
    # calculate shape functions
    @threads for mp ∈ 1:mpD.nmp
        @simd for nn ∈ 1:meD.nn
            # compute basis functions
            id     = mpD.p2n[nn,mp]
            ξ      = (mpD.x[mp,1]-meD.xn[id,1])/meD.h[1] 
            type   = whichType(meD.xn[id,1],meD.xB[1:2],meD.h[1])
            ϕx,dϕx = ϕ∇ϕ(ξ,type,meD.h[1])
            η      = (mpD.x[mp,2]-meD.xn[id,2])/meD.h[2]
            type   = whichType(meD.xn[id,2],meD.xB[3:4],meD.h[2])
            ϕy,dϕy = ϕ∇ϕ(η,type,meD.h[2])
            ζ      = (mpD.x[mp,3]-meD.xn[id,3])/meD.h[3]
            type   = whichType(meD.xn[id,3],meD.xB[5:6],meD.h[3])
            ϕz,dϕz = ϕ∇ϕ(ζ,type,meD.h[3])
            # convolution of basis function
            mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕy*  ϕz                                                                                
            mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕy*  ϕz                                                                                
            mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕy*  ϕz                                   
            mpD.ϕ∂ϕ[nn,mp,4] =  ϕx*  ϕy* dϕz       
        end
        # B-matrix assembly
        mpD.B[1:meD.nD:end,1,mp].= mpD.ϕ∂ϕ[:,mp,2]
        mpD.B[2:meD.nD:end,2,mp].= mpD.ϕ∂ϕ[:,mp,3]
        mpD.B[3:meD.nD:end,3,mp].= mpD.ϕ∂ϕ[:,mp,4]
        mpD.B[2:meD.nD:end,4,mp].= mpD.ϕ∂ϕ[:,mp,4]
        mpD.B[3:meD.nD:end,4,mp].= mpD.ϕ∂ϕ[:,mp,3]
        mpD.B[1:meD.nD:end,5,mp].= mpD.ϕ∂ϕ[:,mp,4]
        mpD.B[3:meD.nD:end,5,mp].= mpD.ϕ∂ϕ[:,mp,2]
        mpD.B[1:meD.nD:end,6,mp].= mpD.ϕ∂ϕ[:,mp,3]
        mpD.B[2:meD.nD:end,6,mp].= mpD.ϕ∂ϕ[:,mp,2]
    end
end
#=
1 xx
2 yy
3 zz
4 yz
5 xz
6 xy
=#