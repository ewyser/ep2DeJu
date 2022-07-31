function RFS(xp,zp,coh0,cohr,phi0,phir)
    # parameters
    θx  = 20.0
    θz  = 2.0
    β   = 45.0*pi/180
    μc  = coh0
    σc  = μc/5.0
    μϕ  = phi0
    σϕ  = μϕ/10.0
    # vector format
    xp  = vec(xp)
    zp  = vec(zp)
    nmp = length(xp)
    # relative distance
    Δx  = (xp.-xp')
    Δz  = zp.-zp'
    # exponential covariance matrix
    if β != 0.0
        C = real.(exp.(-sqrt.(complex.((( Δx.*cos(β).+Δz.*sin(β))./θx).^2+((-Δx.*cos(β).+Δz.*sin(β))./θz).^2))))
    else
        C = real.(exp.(-sqrt.(complex.((Δx./θx).^2+(Δz./θz).^2))))  
    end
    C[diagind(C)].= 1.0    
    cϕ   = cholesky(C).L*randn(Float64,nmp,2)
    p    = 0.5
    R    = [1.0 0.0;p sqrt(1.0-p^2)]
    cϕ   = R*cϕ'
    c    = μc.+σc.*cϕ[1,:]
    ϕ    = μϕ.+σϕ.*cϕ[2,:]

    p    = findall(x->x<=cohr, c)
    c[p].= cohr
	return(c,ϕ)
end

    #=
    #ρ  = exp.(-sqrt.(complex.((  Δx                     ./θx).^2+(  Δz                     ./θz).^2)))
    ρ  = exp.(-(complex.((Δx./θx).^2+(Δz./θz).^2)))
    C  = real.(ρ)
    Q  = eigvecs(C)
    Λ  = diagm(eigvals(C))
    c  = (Q*Λ.^(0.5)*randn(Float64,nmp)).+μ    
    =#