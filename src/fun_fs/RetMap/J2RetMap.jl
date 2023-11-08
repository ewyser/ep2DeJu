@views function J2plast!(mpD,Del,Kc,Hp,fwrkDeform) # Borja (1990); De Souza Neto (2008)
    ftol,ηtol,ηit,ηmax = 1e-6,1e4,0,0
    Hp,χ = 0.35*Hp,3.0/2.0
    # create an alias
    if fwrkDeform == :finite
        σ = mpD.τ
    elseif fwrkDeform == :infinitesimal
        σ = mpD.σ
    end
    @threads for mp in 1:mpD.nmp
        κ = 2.5*mpD.coh[mp]+Hp*mpD.ϵpII[mp]
        cr= mpD.cohr[mp]
        if κ <= cr κ = cr end
        p    = (σ[1,mp]+σ[2,mp]+σ[3,mp])/3.0
        ξ    = σ[:,mp].-[p;p;p;0.0]
        J2   = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2) # Borja (2013), p.33
        ξn   = sqrt(2.0*J2) 
        n    = ξ./ξn
        q    = sqrt(χ)*ξn
        f    = ξn - κ        
        if f>0.0 
            Δλ  = 0.0
            γ   = copy(mpD.ϵpII[mp])
            σ0  = copy(σ[:,mp])
            ϵ0  = copy(mpD.ϵ[:,mp])
            ηit = 1
            while abs(f)>1e-9 && ηit < 20
                ∂f∂σ = n
                Δλ   = f/(∂f∂σ'*Del*∂f∂σ)
                Δσ   = (Δλ*Del*∂f∂σ)        
                σ0 .-= Δσ 
                ϵ0 .-= Del\Δσ
                γ   += Δλ
                p    = (σ0[1]+σ0[2]+σ0[3])/3.0
                ξ    = σ0.-[p;p;p;0.0]
                J2   = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2)
                ξn   = sqrt(2.0*J2)
                n    = ξ./ξn
                q    = sqrt(χ)*ξn
                κ    = 2.5*mpD.coh[mp]+Hp*γ
                if κ <= cr κ = cr end
                f    = ξn - κ
                ηit +=1
            end
            σ[:,mp]     .= σ0
            mpD.ϵ[:,mp] .= ϵ0
            mpD.ϵpII[mp] = γ
            ηmax         = max(ηit,ηmax)
        end
    end
    return ηmax::Int64
end