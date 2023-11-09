@views function getCamCParam(σ0,nstr)
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
    return P,ξ,ξn,J2
end
@views function ∂f(∂f∂p,∂f∂q,n,χ,nstr)
    if nstr == 3
        ∂f∂σ = [∂f∂p*1.0/3.0+sqrt(χ)*∂f∂q*n[1];
                ∂f∂p*1.0/3.0+sqrt(χ)*∂f∂q*n[2];
                ∂f∂p*0.0/3.0+sqrt(χ)*∂f∂q*n[3]]
    elseif nstr == 6

    end
    return ∂f∂σ
end
@views function camCRetMap!(mpD,cmParam,fwrkDeform) # Borja (1990); De Souza Neto (2008)
    ηmax = 20
    χ   = 3.0/2.0
    Pt  = cmParam.Kc/10.0
    a0  = Pt+cmParam.Kc/10.0
    β   = 1.0/1.0
    Pc  = β*a0+(a0-Pt)     
    ϕcs = 20.0*pi/180.0
    M   = 6.0*sin(ϕcs)/(3.0-sin(ϕcs))
    ζ   = 0.0

    # create an alias
    if fwrkDeform == :finite
        σ,nstr = mpD.τ,size(mpD.τ,1)
    elseif fwrkDeform == :infinitesimal
        σ,nstr = mpD.σ,size(mpD.σ,1)
    end

    Ps = zeros(mpD.nmp)
    Qs = zeros(mpD.nmp)
    @threads for p in 1:mpD.nmp
        ϕcs,a     = mpD.ϕ[p],a0*(exp(-ζ*mpD.ϵpV[p]))
        P,ξ,ξn,J2 = getCamCParam(σ[:,p],nstr)
        n,q       = ξ./ξn,sqrt(χ)*ξn
        if P<(Pt-a) b = β else b = 1.0 end
        f         = (1.0/b^2)*(P-Pt+a)^2+(q/M)^2-a^2 # De Souza Neto (2008)
        if f>0.0 
            σ0    = σ[:,p]
            ϵpV,γ = mpD.ϵpV[p],mpD.ϵpII[p]
            Δλ,η  = 0.0,1
            while abs(f)>1e-6 && η < ηmax
                ∂f∂p = (2.0/b^2)*(P-Pt+a)
                ∂f∂q = 2.0*q/M^2
                ∂f∂σ = ∂f(∂f∂p,∂f∂q,n,χ,nstr)   
                Δλ   = f/(∂f∂σ'*cmParam.Del*∂f∂σ)        
                σ0 .-= (Δλ*cmParam.Del*∂f∂σ)  
                ϵpV += Δλ*∂f∂p
                γ   += Δλ*∂f∂q
                a    = a0*(exp(-ζ*ϵpV))
                P,ξ,ξn,J2 = getCamCParam(σ0,nstr)
                n,q       = ξ./ξn,sqrt(χ)*ξn
                if P<(Pt-a) b = β else b = 1.0 end
                f         = (1.0/b^2)*(P-Pt+a)^2+(q/M)^2-a^2 # De Souza Neto (2008)
                η   +=1
            end
        end
        Ps[p]=P
        Qs[p]=q
    end
    println(size(Ps),size(Qs))
    gr() # We will continue onward using the GR backend
    tit = ""
    plot(Ps./(Pc), Qs./abs(Pc), show=true, markershape=:circle,markersize=1.0, color = :blue, seriestype = :scatter, title = tit,xlabel=L"p/p_c",ylabel=L"q/p_c",aspect_ratio=1,xlim=(-1.5,0.5),ylim=(0.0,1.1),)

    return ηmax::Int64
end