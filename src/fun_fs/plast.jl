@views function σTr(σ,nstr)
    σ0 = copy(σ)
    if nstr == 3
        P   = (σ0[1]+σ0[2])/2.0
        τ0  = σ0.-[P,P,0.0]
        τII = sqrt(0.5*(τ0[1]^2+τ0[2]^2)+τ0[4]^2)
    elseif nstr == 6
        P  = (σ0[1]+σ0[2]+σ0[3])/3.0
        τ0     = σ0.-[P,P,P,0.0,0.0,0.0]
        τII    = sqrt(0.5*(τ0[1]^2+τ0[2]^2+τ0[3]^2)+τ0[4]^2+τ0[5]^2+τ0[6]^2)
    end
    return P,σ0,τ0,τII
end
@views function materialParam(ϕ,ψ,c,nstr)
    if nstr == 3
        η   = 3.0*tan(ϕ)/(sqrt(9.0+12.0*tan(ϕ)*tan(ϕ)))
        ηB  = 3.0*tan(ψ)/(sqrt(9.0+12.0*tan(ψ)*tan(ψ)))
        ξ   = 3.0*c     /(sqrt(9.0+12.0*tan(ϕ)*tan(ϕ)))
    elseif nstr == 6
        η   = 6.0  *sin(ϕ)/(sqrt(3.0)*(3.0+sin(ϕ)))
        ηB  = 6.0  *sin(ψ)/(sqrt(3.0)*(3.0+sin(ψ)))
        ξ   = 6.0*c*cos(ϕ)/(sqrt(3.0)*(3.0+sin(ϕ))) 
    end
    return η,ηB,ξ
end
@views function σn(Pn,τ0,τn,τII,nstr)
    if nstr == 3
        σ = τ0.*(τn/τII).+[Pn,Pn,0.0]
    elseif nstr == 6
        σ = τ0.*(τn/τII).+[Pn,Pn,Pn,0.0,0.0,0.0]
    end
    return σ 
end
@views function DPplast!(mpD,cmParam,fwrkDeform)
    ψ,nstr   = 0.0*π/180.0,size(mpD.σ,1)
    # create an alias for stress tensor
    if fwrkDeform == :finite
        σ = mpD.τ
    elseif fwrkDeform == :infinitesimal
        σ = mpD.σ
    end
    # closed-form solution return-mapping for D-P
    for p ∈ 1:mpD.nmp
        c   = mpD.c0[p]+cmParam.Hp*mpD.ϵpII[p]
        if c<mpD.cr[p] c = mpD.cr[p] end
        P,σ0,τ0,τII = σTr(σ[:,p],nstr)
        η,ηB,ξ      = materialParam(mpD.ϕ[p],ψ,c,nstr)
        σm,τP       = ξ/η,ξ-η*(ξ/η)
        fs,ft       = τII+η*P-ξ,P-σm         
        αP,h        = sqrt(1.0+η^2)-η,τII-τP-(sqrt(1.0+η^2))*(P-σm)  
        if fs>0.0 && P<σm || h>0.0 && P>=σm
            Δλ          = fs/(cmParam.Gc+cmParam.Kc*η*ηB)
            Pn,τn       = P-cmParam.Kc*ηB*Δλ,ξ-η*(P-cmParam.Kc*ηB*Δλ)
            σ[:,p]     .= σn(Pn,τ0,τn,τII,nstr)
            mpD.ϵpII[p]+= Δλ*sqrt(1/3+2/9*ηB^2)
            if fwrkDeform == :finite
                mpD.ϵ[:,:,p].= mutate(cmParam.Del\σ[:,p],0.5,:tensor)
                # update left cauchy green tensor
                λ,n          = eigen(mpD.ϵ[:,:,p],sortby=nothing)
                mpD.b[:,:,p].= n*diagm(exp.(2.0.*λ))*n'
            end
        end
        if h<=0.0 && P>=σm
            Δλ          = (P-σm)/cmParam.Kc
            Pn          = σm-P
            σ[:,p]     .= σn(Pn,τ0,0.0,τII,nstr)
            mpD.ϵpII[p]+= sqrt(2.0)*Δλ/3.0
            if fwrkDeform == :finite
                mpD.ϵ[:,:,p].= mutate(cmParam.Del\σ[:,p],0.5,:tensor)
                # update left cauchy green tensor
                λ,n          = eigen(mpD.ϵ[:,:,p],sortby=nothing)
                mpD.b[:,:,p].= n*diagm(exp.(2.0.*λ))*n'
            end
        end
    end
    return 0
end
@views function MCplast!(mpD,cmParam,fwrkDeform)
    ftol,ηtol,ηmax = 1e-6,1e4,0
    ψ              = 0.5*π/180.0
    # create an alias
    if fwrkDeform == :finite
        σ = mpD.τ
    elseif fwrkDeform == :infinitesimal
        σ = mpD.σ
    end
    @threads for p ∈ 1:mpD.nmp
        ϕ,H,ϵII0 = mpD.ϕ[p],cos(mpD.ϕ[p])*cmParam.Hp,mpD.ϵpII[p]
        c0,cr    = mpD.c0[p]+cmParam.Hp*ϵII0,mpD.cr[p]
        if c0<cr c0 = cr end
        σm,τII   = 0.5*(σ[1,p]+σ[2,p]),sqrt(0.25*(σ[1,p]-σ[2,p])^2+σ[3,p]^2)
        f        = τII+σm*sin(ϕ)-c0*cos(ϕ)    
        if f>0.0
            ϵII = ϵII0
            Δϵ  = zeros(Float64,3)
            ηit = 0
            while abs(f)>ftol
                ηit+= 1
                ∂σf = [ (σ[1,p]-σ[2,p])/(4*τII)+sin(ϕ)/2;
                       -(σ[1,p]-σ[2,p])/(4*τII)+sin(ϕ)/2;
                         σ[3,p]/τII                     ]
                ∂σg = [ (σ[1,p]-σ[2,p])/(4*τII)+sin(ψ)/2;
                       -(σ[1,p]-σ[2,p])/(4*τII)+sin(ψ)/2;
                         σ[3,p]/τII                     ] 

                Δγ  = f/(H+∂σf'*cmParam.Del*∂σg)
                Δσ  = Δγ*cmParam.Del*∂σg
                Δϵ.+= cmParam.Del\Δσ
                ϵII = ϵII0+sqrt(2/3*(Δϵ[1]^2+Δϵ[2]^2+2*Δϵ[3]^2))
                c0  = mpD.c0[p]+cmParam.Hp*ϵII
                if c0<cr c0 = cr end
                σ[:,p].-= Δσ
                σm,τII  = 0.5*(σ[1,p]+σ[2,p]),sqrt(0.25*(σ[1,p]-σ[2,p])^2+σ[3,p]^2)
                f       = τII+σm*sin(ϕ)-c0*cos(ϕ)
                if ηit>ηtol
                    @error "CPA: η_it>$(ηit): program killed..."
                    exit(1)
                end
                ηmax = max(ηit,ηmax)
            end
            mpD.ϵ[:,:,p].= mutate(cmParam.Del\σ[:,p],0.5,:tensor)
            mpD.ϵpII[p]  = ϵII 
            if fwrkDeform == :finite
                # update left cauchy green tensor
                λ,n           = eigen(mpD.ϵ[:,:,p],sortby=nothing)
                mpD.b[:,:,p] .= n*diagm(exp.(2.0.*λ))*n'
            end
        end        
    end
    return ηmax::Int64
end
function plast!(mpD,cmParam,cmType,fwrkDeform)
    if cmType == "MC"
        ηmax = MCplast!(mpD,cmParam,fwrkDeform)
    elseif cmType == "DP"        
        ηmax = DPplast!(mpD,cmParam,fwrkDeform)
    else
        err_msg = "$(cmType): invalid plastic model"
        throw(error(err_msg))
    end
    return ηmax::Int64
end