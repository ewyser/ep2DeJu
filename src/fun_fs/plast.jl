@views function MCplast!(mpD,Del,Hp)
    ftol,ηtol,ηmax = 1e-6,1e4,0
    ψ              = 0.5*pi/180.0
    @threads for p ∈ 1:mpD.nmp
        ϕ,H,ϵII0 = mpD.phi[p],cos(mpD.phi[p])*Hp,mpD.ϵpII[p]
        c0       = mpD.coh[p]+Hp*ϵII0
        cr       = mpD.cohr[p]
        if c0<cr c0 = cr end
        σm       = 0.5*(mpD.τ[1,p]+mpD.τ[2,p])
        τII      = sqrt(0.25*(mpD.τ[1,p]-mpD.τ[2,p])^2+mpD.τ[4,p]^2)
        f        = τII+σm*sin(ϕ)-c0*cos(ϕ)    
        if f>0.0
            ϵII = ϵII0
            Δϵ  = zeros(Float64,4)
            ηit = 0
            while abs(f)>ftol
                ηit+= 1
                ∂σf = [ (mpD.τ[1,p]-mpD.τ[2,p])/(4*τII)+sin(ϕ)/2;
                       -(mpD.τ[1,p]-mpD.τ[2,p])/(4*τII)+sin(ϕ)/2;
                        0.0                                     ;
                        mpD.τ[4,p]/τII                      ]
                ∂σg = [ (mpD.τ[1,p]-mpD.τ[2,p])/(4*τII)+sin(ψ)/2;
                       -(mpD.τ[1,p]-mpD.τ[2,p])/(4*τII)+sin(ψ)/2;
                        0.0                                     ;
                        mpD.τ[4,p]/τII                          ] 

                Δγ  = f/(H+∂σf'*Del*∂σg)
                Δσ  = Δγ*Del*∂σg
                Δϵ  = Δϵ+Del\Δσ
                ϵII = ϵII0+sqrt(2/3*(Δϵ[1]^2+Δϵ[2]^2+Δϵ[3]^2+2*Δϵ[4]^2))
                c0  = mpD.coh[p]+Hp*ϵII
                if c0<cr c0 = cr end
                mpD.τ[:,p] .-= Δσ
                σm           = 0.5*(mpD.τ[1,p]+mpD.τ[2,p])
                τII          = sqrt(0.25*(mpD.τ[1,p]-mpD.τ[2,p])^2+mpD.τ[4,p]^2)
                f            = τII+σm*sin(ϕ)-c0*cos(ϕ)
                if ηit>ηtol
                    @printf("\nCPA: max(η_it)>%d",ηtol)
                    @printf("\n     f = %.6f",f)
                    @printf("\n     program killed...")
                    exit(1)
                end
                ηmax = max(ηit,ηmax)
            end
            mpD.ϵ[:,p] .-= Δϵ
            mpD.ϵpII[p]  = ϵII 
        end        
    end
    return ηmax
end
@views function J2plast!(mpD,Del,Kc,Hp) # Borja (1990); De Souza Neto (2008)
    ftol,ηtol,ηit,ηmax = 1e-6,1e4,0,0
    Hp,χ = 0.35*Hp,3.0/2.0
    @threads for mp in 1:mpD.nmp
        κ = 2.5*mpD.coh[mp]+Hp*mpD.ϵpII[mp]
        cr= mpD.cohr[mp]
        if κ <= cr κ = cr end
        p    = (mpD.τ[1,mp]+mpD.τ[2,mp]+mpD.τ[3,mp])/3.0
        ξ    = mpD.τ[:,mp].-[p;p;p;0.0]
        J2   = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2) # Borja (2013), p.33
        ξn   = sqrt(2.0*J2) 
        n    = ξ./ξn
        q    = sqrt(χ)*ξn
        f    = ξn - κ        
        if f>0.0 
            γ  = mpD.ϵpII[mp]
            Δλ = 0.0
            ηit= 1
            τ0 = mpD.τ[:,mp]
            while abs(f)>1e-9 && ηit < 20
                ∂f∂τ = n
                Δλ   = f/(∂f∂τ'*Del*∂f∂τ)        
                τ0 .-= (Δλ*Del*∂f∂τ)  
                γ   += Δλ
                p    = (τ0[1]+τ0[2]+τ0[3])/3.0
                ξ    = τ0[:].-[p;p;p;0.0]
                J2   = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2)
                ξn   = sqrt(2.0*J2)
                n    = ξ./ξn
                q    = sqrt(χ)*ξn
                κ    = 2.5*mpD.coh[mp]+Hp*γ
                if κ <= cr κ = cr end
                f    = ξn - κ
                ηit +=1
            end
            mpD.τ[:,mp] .= τ0
            mpD.ϵpII[mp] = γ
            ηmax         = max(ηit,ηmax)
        end
    end
    return ηmax
end
function plast!(mpD,cmParam,cmType)
    if cmType == "mohr"
        ηmax = MCplast!(mpD,cmParam.Del,cmParam.Hp)
    elseif cmType == "J2"
        ηmax = J2plast!(mpD,cmParam.Del,cmParam.Kc,cmParam.Hp)
    elseif cmType == "camC"
        
    else
        @error "invalid plastic model --"*string(cmType)*"--"
        exit(1) 
    end
    return ηmax
end









































#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function bEplast!(τ,ϵ,epII,coh,phi,nmp,Del,Hp,cr)
    ψ    = 0.5*pi/180.0
    Hp   = 0.0
    ftol = 1e-6
    ηtol = 1e4
    I    = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    @threads for p ∈ 1:nmp
        ϕ   = phi[p]
        H   = cos(ϕ)*Hp
        ϵII0= ϵpII[p]
        c0  = coh[p]+Hp*ϵII0
        if c0<cr
            c0 = cr
        end
        σm  = 0.5*(τ[1,p]+τ[2,p])
        τII = sqrt(0.25*(τ[1,p]-τ[2,p])^2+τ[4,p]^2)
        f   = τII+σm*sin(ϕ)-c0*cos(ϕ)
        e   = [0 f]
        σ0  = τ[:,p]
        σ   = σ0
        if f>0.0
            δϵ  = zeros(Float64,4,1)
            δϵII= 0.0
            δγ  = 0.0
            η   = 1
            while(abs(f)>ftol && η<10)
                ∂σf = [ (σ[1]-σ[2])/(4*τII)+sin(ϕ)/2;
                       -(σ[1]-σ[2])/(4*τII)+sin(ϕ)/2;
                        0.0                         ;
                        σ[4]/τII                    ]
                ∂∂σf= [ 1/(4*τII)                   ;
                       -1/(4*τII)                   ;
                        0.0                         ;
                        1/(1*τII)                   ]
                ∂σg = [ (σ[1]-σ[2])/(4*τII)+sin(ψ)/2;
                       -(σ[1]-σ[2])/(4*τII)+sin(ψ)/2;
                        0.0                         ;
                        σ[4]/τII                    ] 
                ∂∂σg= ∂∂σf
                
                Rσ    = σ-σ0+δγ*Del*∂σg
                #=
                Rf    = f
                R     = [Rσ;Rf]
                ∂R∂σ  = I+δγ*Del*∂∂σg
                ∂R∂γ  = Del*(∂σg+δγ*∂∂σg)
                ∂Rf∂σ = ∂σf'
                ∂Rf∂γ = -h
                B     = [σ;δγ]-[∂R∂σ ∂R∂γ;∂Rf∂σ ∂Rf∂γ]\R
                
                σ     = B[1:4]
                δγ    = B[5]
                δϵ    = Del\(σ0 - σ)
                δϵII  = sqrt(2/3*(δϵ[1]^2+δϵ[2]^2+δϵ[3]^2+2*δϵ[4]^2))
                c0    = coh[p]+Hp*ϵII0
                if c0<cr
                    c0 = cr
                end
                σm  = 0.5*(σ[1]+σ[2])
                τII = sqrt(0.25*(σ[1]-σ[2])^2+σ[4]^2)
                f   = τII+σm*sin(ϕ)-c0*cos(ϕ)
                =#
                η += 1 
            end
            ϵ[:,p] -= δϵ
            ϵpII[p] = δϵII
            τ[:,p]  = σ 
        end
    end
end
function MCnoCPAplast!(τ,ϵ,epII,coh,phi,nmp,Del,Hp,cr)
    for mp ∈ 1:nmp
        c   = coh[mp]+Hp*epII[mp]
        if c<cr
            c = cr
        end
        ϕ  = phi[mp]
        σn = τ[:,mp]
        Δσ = σn[1]-σn[2]
        τt = sqrt(0.25*Δσ^2+σn[4]^2)
        σt = 0.5*(σn[1]+σn[2])
        f  = τt+σt*sin(ϕ)-c*cos(ϕ)
        
        if f>0.0
            if σt <= (c/tan(ϕ))
                β = abs(c*cos(ϕ)-sin(ϕ)*σt)/τt
                σn[1] = σt+β*Δσ/2
                σn[2] = σt-β*Δσ/2
                σn[4] = β*σn[4]
            else
                σn[1] = c/tan(ϕ)
                σn[2] = c/tan(ϕ)
                σn[4] = 0.0
            end
        end
        Δσ       = σn-τ[:,mp]
        τ[:,mp]  = σn
        ϵp       = Del\Δσ
        ϵ[:,mp]  = ϵ[:,mp]+ϵp
        epII[mp]+= sqrt(2/3*(ϵp[1]^2+ϵp[2]^2+ϵp[3]^2)+2*ϵp[4]^2) 
    end
end