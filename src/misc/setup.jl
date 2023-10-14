function meshSetup(nel,L,nD,typeD)
    # geometry                                               
    if nD == 2
        L   = [L[1],ceil(L[2])]
        h   = [L[1]/nel,L[1]/nel]
    elseif nD == 3
        L   = [Lx,Ly,ceil(Lz)]
        h   = [L[1]/nel[1],L[1]/nel[1],L[1]/nel[1]]
    end
    # mesh 
    if nD == 2
        xn  = (0.0-2*h[1]):h[1]:(L[1]+2.0*h[1])
        zn  = (0.0-2*h[2]):h[2]:(L[2]+2.0*h[2])
        zn  = reverse(zn)
        nno = [length(xn) length(zn) length(xn)*length(zn)] 
        nel = [nno[1]-1 nno[2]-1 (nno[1]-1)*(nno[2]-1)]
        nn  = 16
        xn  = (xn'.*ones(typeD,nno[2],1     ))     
        zn  = (     ones(typeD,nno[1],1     )'.*zn)
        xn  = vec(xn)
        zn  = vec(zn)
        x   = hcat(xn,zn)
    elseif nD == 3
        xn  = (0.0-2*h[1]):h[1]:(L[1]+2.0*h[1])
        yn  = (0.0-2*h[2]):h[2]:(L[2]+2.0*h[2])
        zn  = (0.0-2*h[3]):h[3]:(L[3]+2.0*h[3])
        zn  = reverse(zn)
        nno = [length(xn),length(yn),length(zn),length(xn)*length(yn)*length(zn)] 
        nel = [nno[1]-1,nno[2]-1,nno[3]-1,(nno[1]-1)*(nno[2]-1)*(nno[3]-1)]
        nn  = 64
        xn  = (xn'.*ones(typeD,nno[3],1     ))     .*ones(typeD,1,1,nno[2])
        zn  = (     ones(typeD,nno[1],1     )'.*zn).*ones(typeD,1,1,nno[2])
        yn  = (     ones(typeD,nno[3],nno[1]))     .*reshape(yn,1,1,nno[2])
        xn  = vec(xn)
        yn  = vec(yn)
        zn  = vec(zn)
        x   = hcat(xn,yn,zn)
    end

    # boundary conditions
    if nD == 2
        xB  = [minimum(xn)+2*h[1],maximum(xn)-2*h[1],0.0,Inf]                                    
        bcx = vcat(findall(x->x<=xB[1], xn),findall(x->x>=xB[2], xn))
        bcz = findall(x->x<=xB[3], zn)
        bcX = ones(Int64,nno[nD+1],1)
        bcX[bcx] .= 0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0
        #bcX[bcz] .= 0
        bc   = hcat(bcX,bcZ)
    elseif nD == 3
        xB  = [minimum(xn)+2*h[1],maximum(xn)-2*h[1],minimum(yn)+2*h[1],maximum(yn)-2*h[1],0.0,Inf]                                    
        bcx = vcat(findall(x->x<=xB[1], xn),findall(x->x>=xB[2], xn))
        bcy = vcat(findall(x->x<=xB[3], yn),findall(x->x>=xB[4], yn))
        bcz = findall(x->x<=xB[5], zn)
        bcX = ones(Int64,nno[nD+1],1)
        bcX[bcx] .= 0
        bcY = ones(nno[nD+1],1)
        bcY[bcy] .= 0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0
        bc   = hcat(bcX,bcY,bcZ)
    end
    # push to named-Tuple
    meD = (
        nD   = nD,
        nel  = nel,
        nno  = nno,
        nn   = nn,
        L    = L,
        h    = h,
        x    = x,
        # nodal quantities
        m    = zeros(typeD,nno[nD+1],1), 
        fext = zeros(typeD,nno[nD+1],nD), 
        fint = zeros(typeD,nno[nD+1],nD),
        D    = zeros(typeD,nno[nD+1],nD),
        f    = zeros(typeD,nno[nD+1],nD),
        a    = zeros(typeD,nno[nD+1],nD),
        p    = zeros(typeD,nno[nD+1],nD),
        v    = zeros(typeD,nno[nD+1],nD),
        u    = zeros(typeD,nno[nD+1],nD),
        pel  = zeros(typeD,nno[nD+1],nD),
        ΔJ   = zeros(typeD,nno[nD+1],nD),
        # mesh-to-node topology
        e2n = e2N(nD,nno,nel,nn),
        xB   = xB,
        bc   = bc,
    )
    return meD
end
function continuum(xp,zp,x,z,c,wl,a,nD)
    xlt = Float64[]
    zlt = Float64[]
    clt = Float64[]
    pos = Float64 
    for mp ∈ 1:length(xp)
        for p ∈ 1:length(z)
            Δx = xp[mp]-x[p]
            Δz = zp[mp]-z[p]
            nx = a
            nz = -1.0
            s  = Δx*nx+Δz*nz        
            if s>0
                pos = 1
            else
                pos = 0
            end
            if zp[mp]<wl 
                pos = 1
            end
        end
        if pos==1
            push!(xlt, xp[mp]) # push!(inArray, What), incremental construction of an array of arbitrary size
            push!(zlt, zp[mp]) # push!(inArray, What), incremental construction of an array of arbitrary size
            push!(clt, c[mp])
        end
    end
    xp = if nD == 2 hcat(xlt,zlt) elseif nD == 3 hcat(xlt,ylt,zlt) end
    return xp,clt
end
function pointSetup(meD,ni,lz,coh0,cohr,phi0,phir,rho0,nstr,typeD)
    # mpm initialization
    xL          = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
    zL          = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:lz-0.5*meD.h[2]/ni
    npx,npz     = length(xL),length(zL)
    xp,zp       = ((xL'.*ones(npz,1  )      )),((     ones(npx,1  )'.*zL )) 
    c           = GRFS_gauss(xp,coh0,cohr,ni,meD.h[1])
    xp,zp,c     = vec(xp),vec(zp),vec(c)
    wl          = 0.15*lz
    x           = LinRange(minimum(xp),maximum(xp),200)
    a           = -1.25
    x,z         = x.+0.5.*meD.L[1],a.*x
    xp,clt      = continuum(xp,zp,x,z,c,wl,a,meD.nD)
    # material point's quantities
    # scalars & vectors
    nmp         = size(xp,1)
    l0          = ones(typeD,nmp,2).*0.5.*(meD.h[1]./ni)
    l           = ones(typeD,nmp,2).*0.5.*(meD.h[1]./ni)
    v0          = ones(typeD,nmp,1).*(2.0.*l0[:,1].*2.0.*l0[:,2])
    v           = ones(typeD,nmp,1).*(2.0.*l[:,1].*2.0.*l[:,2])
    m           = rho0.*v0
    coh         =  ones(typeD,nmp,1).*coh0#clt
    #coh  =  clt
    #coh,phi  = RFS(xp[:,1],xp[:,2],coh0,cohr,phi0,phir)
    cohr        =  ones(typeD,nmp,1).*cohr
    phi         =  ones(typeD,nmp,1).*phi0
    p           = findall(x->x<=2*wl, xp[:,2])
    phi[p]     .= phir
    # push/init. to mpD()::NamedTuple data structure 
    mpD = (
        nmp  = nmp,
        x    = xp,
        u    = zeros(typeD,nmp,meD.nD), 
        v    = zeros(typeD,nmp,meD.nD),
        p    = zeros(typeD,nmp,meD.nD),
        l0   = l0,
        l    = l,
        V0   = v0,
        V    = v,
        m    = m,
        coh  = coh,
        cohr = cohr,
        phi  = phi,
        ϵpII = zeros(typeD,nmp,1),
        ϵpV  = zeros(typeD,nmp,1), 
        ΔJ   = ones(typeD,nmp,1),
        J    = ones(typeD,nmp,1),
        # tensor in matrix notation
        ΔF   = zeros(typeD,meD.nD,meD.nD,nmp),
        ΔFbar= zeros(typeD,meD.nD,meD.nD,nmp),
        F    = repeat(Matrix(1.0I,meD.nD,meD.nD),1,1,nmp),
        ϵ    = zeros(typeD,meD.nD,meD.nD,nmp),
        b    = zeros(typeD,meD.nD,meD.nD,nmp),
        bT   = zeros(typeD,meD.nD,meD.nD,nmp),
        # tensor in voigt notation
        ω    = zeros(typeD,1,nmp),
        σR   = zeros(typeD,nstr,nmp),
        σ    = zeros(typeD,nstr,nmp),
        τ    = zeros(typeD,nstr,nmp),
        dev  = zeros(typeD,nstr,nmp),
        ep   = zeros(typeD,nstr,nmp),
        # additional quantities
        ϕ∂ϕ  = zeros(typeD,nmp ,meD.nn,meD.nD+1   ),
        B    = zeros(typeD,nstr,meD.nn.*meD.nD,nmp),
        # connectivity
        p2e  = zeros(Int64,nmp,1),
        p2n  = zeros(Int64,nmp,meD.nn),
    )
    return mpD 
end
function e2N(nD,nno,nel,nn)
	e2n  = zeros(nel[nD+1],nn)
    if nD == 2
        gnum = reverse(reshape(1:(nno[3]),nno[2],nno[1]),dims=1)
        iel  = 1
        for i ∈ 1:nel[1]#nelx
            for j ∈ 1:nel[2]#nelz
                if i>1 && i<nel[1] && j>1 && j<nel[2]
                    e2n[iel,1 ] = gnum[j-1,i-1]
                    e2n[iel,2 ] = gnum[j-0,i-1]
                    e2n[iel,3 ] = gnum[j+1,i-1]
                    e2n[iel,4 ] = gnum[j+2,i-1]

                    e2n[iel,5 ] = gnum[j-1,i  ]
                    e2n[iel,6 ] = gnum[j-0,i  ]
                    e2n[iel,7 ] = gnum[j+1,i  ]
                    e2n[iel,8 ] = gnum[j+2,i  ]

                    e2n[iel,9 ] = gnum[j-1,i+1]
                    e2n[iel,10] = gnum[j-0,i+1]
                    e2n[iel,11] = gnum[j+1,i+1]
                    e2n[iel,12] = gnum[j+2,i+1]

                    e2n[iel,13] = gnum[j-1,i+2]
                    e2n[iel,14] = gnum[j-0,i+2]
                    e2n[iel,15] = gnum[j+1,i+2]
                    e2n[iel,16] = gnum[j+2,i+2]
                end
                iel = iel+1;
            end
        end
    elseif nD == 3
        gnum = reverse(reshape(1:(nno[4]),nno[3],nno[1],nno[2]),dims=1)
        iel  = 1
        for k ∈ 1:nel[2]#nely
            for i ∈ 1:nel[1]#nelx
                for j ∈ 1:nel[3]#nelz
                    if i>1 && i<nel[1] && j>1 && j<nel[3] && k>1 && k<nel[2]
                        e2n[iel,1 ] = gnum[j-1,i-1,k-1]
                        e2n[iel,2 ] = gnum[j-0,i-1,k-1]
                        e2n[iel,3 ] = gnum[j+1,i-1,k-1]
                        e2n[iel,4 ] = gnum[j+2,i-1,k-1]
                        e2n[iel,5 ] = gnum[j-1,i  ,k-1]
                        e2n[iel,6 ] = gnum[j-0,i  ,k-1]
                        e2n[iel,7 ] = gnum[j+1,i  ,k-1]
                        e2n[iel,8 ] = gnum[j+2,i  ,k-1]
                        e2n[iel,9 ] = gnum[j-1,i+1,k-1]
                        e2n[iel,10] = gnum[j-0,i+1,k-1]
                        e2n[iel,11] = gnum[j+1,i+1,k-1]
                        e2n[iel,12] = gnum[j+2,i+1,k-1]
                        e2n[iel,13] = gnum[j-1,i+2,k-1]
                        e2n[iel,14] = gnum[j-0,i+2,k-1]
                        e2n[iel,15] = gnum[j+1,i+2,k-1]
                        e2n[iel,16] = gnum[j+2,i+2,k-1]
                        
                        e2n[iel,17] = gnum[j-1,i-1,k  ]
                        e2n[iel,18] = gnum[j-0,i-1,k  ]
                        e2n[iel,19] = gnum[j+1,i-1,k  ]
                        e2n[iel,20] = gnum[j+2,i-1,k  ]
                        e2n[iel,21] = gnum[j-1,i  ,k  ]
                        e2n[iel,22] = gnum[j-0,i  ,k  ]
                        e2n[iel,23] = gnum[j+1,i  ,k  ]
                        e2n[iel,24] = gnum[j+2,i  ,k  ]
                        e2n[iel,25] = gnum[j-1,i+1,k  ]
                        e2n[iel,26] = gnum[j-0,i+1,k  ]
                        e2n[iel,27] = gnum[j+1,i+1,k  ]
                        e2n[iel,28] = gnum[j+2,i+1,k  ]
                        e2n[iel,29] = gnum[j-1,i+2,k  ]
                        e2n[iel,30] = gnum[j-0,i+2,k  ]
                        e2n[iel,31] = gnum[j+1,i+2,k  ]
                        e2n[iel,32] = gnum[j+2,i+2,k  ]
                        
                        e2n[iel,33] = gnum[j-1,i-1,k+1]
                        e2n[iel,34] = gnum[j-0,i-1,k+1]
                        e2n[iel,35] = gnum[j+1,i-1,k+1]
                        e2n[iel,36] = gnum[j+2,i-1,k+1]
                        e2n[iel,37] = gnum[j-1,i  ,k+1]
                        e2n[iel,38] = gnum[j-0,i  ,k+1]
                        e2n[iel,39] = gnum[j+1,i  ,k+1]
                        e2n[iel,40] = gnum[j+2,i  ,k+1]
                        e2n[iel,41] = gnum[j-1,i+1,k+1]
                        e2n[iel,42] = gnum[j-0,i+1,k+1]
                        e2n[iel,43] = gnum[j+1,i+1,k+1]
                        e2n[iel,44] = gnum[j+2,i+1,k+1]
                        e2n[iel,45] = gnum[j-1,i+2,k+1]
                        e2n[iel,46] = gnum[j-0,i+2,k+1]
                        e2n[iel,47] = gnum[j+1,i+2,k+1]
                        e2n[iel,48] = gnum[j+2,i+2,k+1]
                            
                        e2n[iel,49] = gnum[j-1,i-1,k+2]
                        e2n[iel,50] = gnum[j-0,i-1,k+2]
                        e2n[iel,51] = gnum[j+1,i-1,k+2]
                        e2n[iel,52] = gnum[j+2,i-1,k+2]
                        e2n[iel,53] = gnum[j-1,i  ,k+2]
                        e2n[iel,54] = gnum[j-0,i  ,k+2]
                        e2n[iel,55] = gnum[j+1,i  ,k+2]
                        e2n[iel,56] = gnum[j+2,i  ,k+2]
                        e2n[iel,57] = gnum[j-1,i+1,k+2]
                        e2n[iel,58] = gnum[j-0,i+1,k+2]
                        e2n[iel,59] = gnum[j+1,i+1,k+2]
                        e2n[iel,60] = gnum[j+2,i+1,k+2]
                        e2n[iel,61] = gnum[j-1,i+2,k+2]
                        e2n[iel,62] = gnum[j-0,i+2,k+2]
                        e2n[iel,63] = gnum[j+1,i+2,k+2]
                        e2n[iel,64] = gnum[j+2,i+2,k+2]
                    end
                    iel = iel+1;
                end
            end
        end
    end
	return convert(Array{Int64},e2n)
end