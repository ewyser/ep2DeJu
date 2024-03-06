@kernel inbounds = true function kernel_ΔJn(mpD,meD,arg)
    k = @index(Global)
    if arg == :p2n && k≤mpD.nmp 
        # accumulation
        for nn ∈ 1:meD.nn
            @atom meD.ΔJn[mpD.p2n[nn,k]]+= mpD.ϕ∂ϕ[nn,k,1]*(mpD.m[k]*mpD.ΔJ[k])  
        end
    elseif arg == :solve && k≤meD.nno[end] 
        # solve
        meD.mn[k]>0.0 ? meD.ΔJn[k]/= meD.mn[k] : meD.ΔJn[k] = 0.0
    elseif arg == :n2p && k≤mpD.nmp 
        # mapping back to mp's
        @views mpD.ΔF[:,:,k].*= (dot(mpD.ϕ∂ϕ[:,k,1],meD.ΔJn[mpD.p2n[:,k]])/mpD.ΔJ[k]).^(1.0/meD.nD)
    end
end
function ΔFbar!(mpD,meD)
    # init mesh quantities to zero
    meD.ΔJn.= 0.0
    # calculate dimensional cst.
    dim     = 1.0/meD.nD
    # action
    @isdefined(ΔJn!) ? nothing : ΔJn! = kernel_ΔJn(CPU())
    # mapping to mesh
    ΔJn!(mpD,meD,:p2n; ndrange=mpD.nmp);sync(CPU())
    # compute nodal determinant of incremental deformation 
    ΔJn!(mpD,meD,:solve; ndrange=meD.nno[end]);sync(CPU())
    # compute determinant Jbar 
    ΔJn!(mpD,meD,:n2p; ndrange=mpD.nmp);sync(CPU())
    return nothing
end
@views function domainUpd!(mpD)
    for p ∈ 1:mpD.nmp
        # update material point's domain length using symmetric material stretch tensor U
        λ,n        = eigen(mpD.F[:,:,p]'*mpD.F[:,:,p],sortby=nothing)
        U          = (n*diagm(sqrt.(λ))*n')
        mpD.l[p,:].= U*mpD.l0[p,:]
    end
    return nothing
end
@views function deform!(mpD,meD,Δt)
    for p ∈ 1:mpD.nmp
        # compute velocity & displacement gradients
        mpD.∇v[:,:,p].= (permutedims(mpD.ϕ∂ϕ[:,p,2:end],(2,1))*meD.vn[mpD.p2n[:,p],:])'
        mpD.∇u[:,:,p].= Δt.*mpD.∇v[:,:,p]
        # compute incremental deformation gradient
        mpD.ΔF[:,:,p].= mpD.I.+mpD.∇u[:,:,p]
        mpD.ΔJ[p]     = det(mpD.ΔF[:,:,p])
        # update deformation gradient
        mpD.F[:,:,p] .= mpD.ΔF[:,:,p]*mpD.F[:,:,p]
        # update material point's volume
        mpD.J[p]      = det(mpD.F[:,:,p])
        mpD.V[p]      = mpD.J[p]*mpD.V0[p]
    end
    return nothing
end
@views function mutate(ϵ,Χ,type)
    if type == :tensor # α = 1/2 when ϵ := strain, α = 1.0 when ϵ := stress
        if size(ϵ) == (3,)
            ϵmut = [  ϵ[1] Χ*ϵ[3];
                    Χ*ϵ[3]   ϵ[2]]
        elseif size(ϵ) == (6,)
            ϵmut = [  ϵ[1] Χ*ϵ[6] Χ*ϵ[5];
                    Χ*ϵ[6]   ϵ[2] Χ*ϵ[4];
                    Χ*ϵ[5] Χ*ϵ[4]   ϵ[3]]
        end
    elseif type == :voigt # α = 2.0 when ϵ := strain, α = 1.0 when ϵ := stress
        if size(ϵ) == (2,2)
            ϵmut = vcat(ϵ[1,1],ϵ[2,2],Χ*ϵ[1,2]) #xx,yy,zz,xy
        elseif size(ϵ) == (3,3)
            ϵmut = vcat(ϵ[1,1],ϵ[2,2],ϵ[3,3],Χ*ϵ[2,3],Χ*ϵ[1,3],Χ*ϵ[1,2]) #xx,yy,zz,yz,xz,xy
        end
    end
    return ϵmut
end
@kernel inbounds = true function kernel_elast(mpD,Del,fwrkDeform)
    p = @index(Global)
    # deformation framework dispatcher
    if fwrkDeform == :finite
        if p ≤ mpD.nmp 
            # update left cauchy-green tensor
            mpD.b[:,:,p].= mpD.ΔF[:,:,p]*mpD.b[:,:,p]*mpD.ΔF[:,:,p]'
            # compute logarithmic strain tensor
            λ,n          = eigen(mpD.b[:,:,p],sortby=nothing)
            mpD.ϵ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
            # krichhoff stress tensor
            mpD.τ[:,p]   = Del*mutate(mpD.ϵ[:,:,p],2.0,:voigt)
        end
    elseif fwrkDeform == :infinitesimal
        if p ≤ mpD.nmp 
            # calculate elastic strains & spin(s)
            mpD.ϵ[:,:,p] .= 0.5.*(mpD.ΔF[:,:,p]+mpD.ΔF[:,:,p]').-mpD.I
            mpD.ω[:,:,p] .= 0.5.*(mpD.ΔF[:,:,p]-mpD.ΔF[:,:,p]')
            # update cauchy stress tensor
            mpD.σJ[:,:,p].= mutate(mpD.σ[:,p],1.0,:tensor)
            mpD.σJ[:,:,p].= mpD.σJ[:,:,p]*mpD.ω[:,:,p]'+mpD.σJ[:,:,p]'*mpD.ω[:,:,p]
            mpD.σ[:,p]  .+= Del*mutate(mpD.ϵ[:,:,p],2.0,:voigt).+mutate(mpD.σJ[:,:,p],1.0,:voigt)
        end   
    end
end
function plast!(mpD,cmParam,cmType,fwrkDeform)
    # plastic return-mapping dispatcher
    if cmType == "MC"
        ηmax = MCRetMap!(mpD,cmParam,fwrkDeform)
    elseif cmType == "DP"        
        ηmax = DPRetMap!(mpD,cmParam,fwrkDeform)
    elseif cmType == "J2"
        ηmax = J2RetMap!(mpD,cmParam,fwrkDeform)
    elseif cmType == "camC"
        ηmax = camCRetMap!(mpD,cmParam,fwrkDeform)
    else
        err_msg = "$(cmType): invalid return mapping for plastic correction"
        throw(error(err_msg))
    end
    return ηmax::Int64
end
@views function elastoplast!(mpD,meD,cmParam,cmType,Δt,ϕ∂ϕType,isΔFbar,fwrkDeform,plastOn)
    # get incremental deformation tensor & strains
    deform!(mpD,meD,Δt)
    # update material point's domain
    ϕ∂ϕType == :gimpm ? domainUpd!(mpD) : nothing
    # volumetric locking correction
    isΔFbar ? ΔFbar!(mpD,meD) : nothing
    # update kirchoff/cauchy stresses
    @isdefined(elastK!) ? nothing : elastK! = kernel_elast(CPU())
    elastK!(mpD,cmParam.Del,fwrkDeform; ndrange=mpD.nmp);sync(CPU())
    # plastic corrector
    if plastOn 
        ηmax = plast!(mpD,cmParam,cmType,fwrkDeform) 
    else 
        ηmax = 0 
    end
    # get cauchy stresses
    if fwrkDeform == :finite
        for p ∈ 1:mpD.nmp
            mpD.σ[:,p] .= mpD.τ[:,p]./mpD.J[p]
        end
    end
    return ηmax::Int64
end