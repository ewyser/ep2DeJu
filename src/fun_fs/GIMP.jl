function NdN(δx::Float64,h::Float64,lp::Float64)                                                         
    if     abs(δx) <    lp                       
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
    return(ϕ,∂ϕ)    
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function ϕ∂ϕ!(B,ϕ,∂ϕx,∂ϕz,σ,τ,J,xp,l,xn,zn,p2n,h,xB,nn,nmp)
    Threads.@threads for mp in 1:nmp
        for nn in 1:nn
            # compute basis functions
            id     = p2n[mp,nn]
            ξ      = xp[mp,1] - xn[id] 
            η      = xp[mp,2] - zn[id]
            ϕx,dϕx = NdN(ξ,h[1],l[mp,1])
            ϕz,dϕz = NdN(η,h[2],l[mp,2])
            # convolution of basis function
            ϕ[mp,nn]   =  ϕx*  ϕz                                        
            ∂ϕx[mp,nn] = dϕx*  ϕz                                        
            ∂ϕz[mp,nn] =  ϕx* dϕz
        end
        # B-matrix assembly
        B[1,1:2:end,mp] = ∂ϕx[mp,:]
        B[2,2:2:end,mp] = ∂ϕz[mp,:]
        B[4,1:2:end,mp] = ∂ϕz[mp,:]
        B[4,2:2:end,mp] = ∂ϕx[mp,:]
        σ[:,mp]         = τ[:,mp]/J[mp]
    end
end