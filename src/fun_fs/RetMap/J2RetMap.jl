@views function getInitial(σ0,nstr)
    if nstr == 3
        P  = (σ0[1]+σ0[2])/2.0
        ξ  = σ0.-[P,P,0.0]
        J2 = 0.5*(ξ[1]^2+ξ[2]^2+2.0*ξ[3]^2) # Borja (2013), p.33
        ξn = sqrt(2.0*J2) 
    elseif nstr == 6
        P  = (σ0[1]+σ0[2]+σ0[3])/3.0
        ξ  = σ0.-[P,P,P,0.0,0.0,0.0]
        J2 = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2+2.0*ξ[5]^2+2.0*ξ[6]^2) # Borja (2013), p.33
        ξn = sqrt(2.0*J2) 
    end
    return ξ,ξn,J2
end
@views function J2RetMap!(mpD,cmParam,fwrkDeform) # Borja (1990); De Souza Neto (2008)
    ftol,ηtol,ηmax = 1e-9,1e4,20
    Hp,χ = 0.35*cmParam.Hp,3.0/2.0
    # create an alias
    if fwrkDeform == :finite
        σ,nstr = mpD.τ,size(mpD.τ,1)
    elseif fwrkDeform == :infinitesimal
        σ,nstr = mpD.σ,size(mpD.σ,1)
    end
    @threads for p in 1:mpD.nmp
        κ = 2.5*mpD.c0[p]+cmParam.Hp*mpD.ϵpII[p]
        if κ <= mpD.cr[p] κ = mpD.cr[p] end
        ξ,ξn,J2 = getInitial(σ[:,p],nstr)
        n,q,f   = ξ./ξn,sqrt(χ)*ξn,ξn-κ
        if f>0.0 
            γ0,σ0 = copy(mpD.ϵpII[p]),copy(σ[:,p])
            ηit  = 1
            while abs(f)>ftol && ηit<ηmax
                ∂f∂σ    = n
                Δλ      = f/(∂f∂σ'*cmParam.Del*∂f∂σ)
                Δσ      = (Δλ*cmParam.Del*∂f∂σ)        
                σ0    .-= Δσ 
                γ0     += Δλ
                ξ,ξn,J2 = getInitial(σ0,nstr)
                κ       = 2.5*mpD.c0[p]+cmParam.Hp*mpD.ϵpII[p]
                if κ <= mpD.cr[p] κ = mpD.cr[p] end
                n,q,f   = ξ./ξn,sqrt(χ)*ξn,ξn-κ
                ηit +=1
            end
            mpD.ϵpII[p] = γ0
            σ[:,p]     .= σ0
            if fwrkDeform == :finite
                # update strain tensor
                mpD.ϵ[:,:,p].= mutate(cmParam.Del\σ[:,p],0.5,:tensor)
                # update left cauchy green tensor
                λ,n          = eigen(mpD.ϵ[:,:,p],sortby=nothing)
                mpD.b[:,:,p].= n*diagm(exp.(2.0.*λ))*n'
            end
            ηmax = max(ηit,ηmax)
        end
    end
    return ηmax::Int64
end

