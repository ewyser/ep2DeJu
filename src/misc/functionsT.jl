function kwargsOut(kwargs)
    if isempty(kwargs)
        # ϵp2De(40,"P","mohr")
        ϕ∂ϕType,fwrkDeform,isΔFbar = "bsmpm","finite",true
    else
        #ϵp2De(40,"P","mohr";shpfun="bsmpm",fwrk="finite",vollock=true)
        shpfun,fwrk,vollock = kwargs
        ϕ∂ϕType,fwrkDeform,isΔFbar = shpfun[2],fwrk[2],vollock[2]
    end
    return ϕ∂ϕType,fwrkDeform,isΔFbar
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
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
    # nodal quantities
    mn  = zeros(typeD,nno[nD+1],1) 
    fen = zeros(typeD,nno[nD+1],nD) 
    fin = zeros(typeD,nno[nD+1],nD)
    fn  = zeros(typeD,nno[nD+1],nD)
    an  = zeros(typeD,nno[nD+1],nD)
    pn  = zeros(typeD,nno[nD+1],nD)
    vn  = zeros(typeD,nno[nD+1],nD)
    un  = zeros(typeD,nno[nD+1],nD)
    pel = zeros(typeD,nno[nD+1],nD)
    ΔJn = zeros(typeD,nno[nD+1],nD)
    # mesh-to-node topology
    e2n = e2N(nD,nno,nel,nn)
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
        m    = mn,
        fext = fen,
        fint = fin,
        f    = fn,
        a    = an,
        p    = pn,
        v    = vn,
        u    = un,
        pel  = pel,
        ΔJ   = ΔJn,
        e2n  = e2n,
        xB   = xB,
        bc   = bc,
    )
    return meD
end
function pointSetup(meD,ni,lz,coh0,cohr,phi0,phir,rho0,nstr,typeD)
    # mpm initialization
    xL  = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
    zL  = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:lz-0.5*meD.h[2]/ni
    npx = length(xL)
    npz = length(zL)
    xp  = ((xL'.*ones(npz,1  )      ))
    zp  = ((     ones(npx,1  )'.*zL ))
    c   = GRFS_gauss(xp,coh0,cohr,ni,meD.h[1])

    xp  = vec(xp)
    zp  = vec(zp)
    c   = vec(c)

    wl  = 0.15*lz
    x   = LinRange(minimum(xp),maximum(xp),200)
    a   = -1.25
    z   = a.*x
    x   = x.+0.5.*meD.L[1]

    nx  = size(xp,2) 
    nz  = size(zp,1)
  
    xlt = Float64[]
    zlt = Float64[]
    clt = Float64[]
    pos = Float64 
    for mp in 1:length(xp)
        for p in 1:length(z)
            Δx = xp[mp]-x[p]
            Δz = zp[mp]-z[p]
            nx = a
            nz = -1.0
            s  = Δx*nx+Δz*nz        
            if(s>0)
                pos = 1
            else
                pos = 0
            end
            if(zp[mp]<wl) 
                pos = 1
            end
        end
        if(pos==1)
            push!(xlt, xp[mp]) # push!(inArray, What), incremental construction of an array of arbitrary size
            push!(zlt, zp[mp]) # push!(inArray, What), incremental construction of an array of arbitrary size
            push!(clt, c[mp])
        end
    end
    #scatter!(xlt,zlt,markershape=:circle,label="",show=true,aspect_ratio=1)
    # material point's quantities
    # scalars & vectors
    nmp  = length(xlt)
    l0   =  ones(typeD,nmp,2).*0.5.*(meD.h[1]./ni)
    l    =  ones(typeD,nmp,2).*0.5.*(meD.h[1]./ni)
    v0   =  ones(typeD,nmp,1).*(2.0.*l0[:,1].*2.0.*l0[:,2])
    v    =  ones(typeD,nmp,1).*(2.0.*l[:,1].*2.0.*l[:,2])
    m    = rho0.*v0
    xp   = hcat(xlt,zlt)
    up   = zeros(typeD,nmp,2) 
    vp   = zeros(typeD,nmp,2)
    pp   = zeros(typeD,nmp,2)
    coh  =  ones(typeD,nmp,1).*coh0#clt
    #coh  =  clt
    #coh,phi  = RFS(xp[:,1],xp[:,2],coh0,cohr,phi0,phir)
    cohr =  ones(typeD,nmp,1).*cohr
    phi  =  ones(typeD,nmp,1).*phi0
    p    = findall(x->x<=2*wl, xp[:,2])
    phi[p] .= phir

    epII = zeros(typeD,nmp,1)
    epV  = zeros(typeD,nmp,1)
    ΔJ   = ones(typeD,nmp,1)
    J    = ones(typeD,nmp,1)
    ΔJbar= ones(typeD,nmp,1)
    Jbar = ones(typeD,nmp,1)
    # tensors
    dF   = zeros(2,2,nmp)
    dFbar= zeros(2,2,nmp)
    F    = ones(2,2,nmp)
    F[1,2,:] .= F[2,1,:] .= 0.0
    Fbar = F
    b    = zeros(2,2,nmp)
    bT   = zeros(2,2,nmp)
    e    = zeros(typeD,nstr,nmp)
    ome  = zeros(typeD,1,nmp)
    s    = zeros(typeD,nstr,nmp)
    τ    = zeros(typeD,nstr,nmp)
    dev  = zeros(typeD,nstr,nmp)
    ep   = zeros(typeD,nstr,nmp)
    # additional quantities
    nn   = convert(UInt64,meD.nn)
    ϕ    = zeros(typeD,nmp ,nn       )
    ∂ϕx  = zeros(typeD,nmp ,nn       )
    ∂ϕz  = zeros(typeD,nmp ,nn       )
    ϕ∂ϕ  = zeros(typeD,nmp ,nn,3     )
    B    = zeros(typeD,nstr,nn.*2,nmp)
    # connectivity
    p2e  = zeros(UInt64,nmp,1)
    p2n  = zeros(UInt64,nmp,nn)
    # push to named-Tuple
    mpD = (
        nmp = nmp,
        l0  = l0,
        l   = l,
        V0  = v0,
        V   = v,
        m   = m,
        x   = xp,
        u   = up,
        v   = vp,
        p   = pp,
        coh = coh,
        cohr= cohr,
        phi = phi,
        ϵpII= epII,
        ϵpV = epV, 
        ΔJ  = ΔJ,
        J   = J,
        ΔF  = dF,
        ΔFbar= dFbar,
        F    = F,
        Fbar = Fbar,
        b    = b,
        bT   = bT,
        ϵ    = e,
        ω    = ome,
        σ    = s,
        τ    = τ,
        dev  = dev,
        ep   = ep,
        ϕ∂ϕ  = ϕ∂ϕ,
        B    = B,
        p2e  = p2e,
        p2n  = p2n,
    )
    return mpD 
end
function e2N(nD,nno,nel,nn)
	e2n  = zeros(nel[nD+1],nn)
    if nD == 2
        gnum = reverse(reshape(1:(nno[3]),nno[2],nno[1]),dims=1)
        iel  = 1
        for i in 1:nel[1]#nelx
            for j in 1:nel[2]#nelz
                if(i>1 && i<nel[1] && j>1 && j<nel[2])
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
        for k in 1:nel[2]#nely
            for i in 1:nel[1]#nelx
                for j in 1:nel[3]#nelz
                    if(i>1 && i<nel[1] && j>1 && j<nel[3] && k>1 && k<nel[2])
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
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
@views function get_vals(meD,mpD,it,ηmax,ηtot,cmpl,symb)
    # completion [%]
    cmpl = round(100.0*cmpl,digits=1)
    # save vals
    vals = [("[nel,np]",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),
            ("iteration(s)",it),
            ("ηmax,ηtot",(ηmax,ηtot)),
            (symb*" t/T",cmpl)]
    return vals
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
@views function D(E,ν)
    Gc = E/(2.0*(1.0+ν))                                                   # shear modulus               [Pa]
    Kc = E/(3.0*(1.0-2.0*ν))                                               # bulk modulus                [Pa]
    D  = [ Kc+4/3*Gc Kc-2/3*Gc Kc-2/3*Gc 0.0 ;
             Kc-2/3*Gc Kc+4/3*Gc Kc-2/3*Gc 0.0 ;
             Kc-2/3*Gc Kc-2/3*Gc Kc+4/3*Gc 0.0 ;
             0.0       0.0       0.0       Gc  ]
    return Kc,Gc,D
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
@views function get_Δt(vp,h,yd)
    Δx   = h[1]
    Δz   = h[2]
    vmax = [abs.(vp[:,1]) abs.(vp[:,2])]
    cmax = [maximum(vmax[:,1]) maximum(vmax[:,2])]
    cmax = [Δx/(cmax[1]+yd) Δz/(cmax[2]+yd)]
    Δt   = 0.5*maximum(cmax)
    return Δt
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function get_g(tw::Float64,tg::Float64,nD::Int64)
    g = 0.0
    if tw<=tg 
        g = 9.81*tw/tg
    else
        g = 9.81
    end
    return if nD == 2 g = [0.0 g] elseif nD == 3 g = [0.0 0.0 g] end
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
default(
    fontfamily="Computer Modern",
    titlefont=12, 
    guidefont=12,  
    tickfont=10, 
    legendfont=10,
    linewidth=2,
    framestyle=:box,
    label=nothing,
    grid=false
    )
@views function plot_coh(xp,coh,phi,ϕ0)
    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    scatter(xp[:,1],xp[:,2],zcolor=coh./1e3,
    markershape=:circle,
    label="",
    show=true,
    aspect_ratio=1,
    c=:vik,
    clims=(10.0,30.0),
    markerstrokecolor=:auto,
    markerstrokewidth=0,
    ylim=(-10,20),
    )
    savefig(path_plot*"coh0.png")
    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    scatter(xp[:,1],xp[:,2],zcolor=phi,
    markershape=:circle,
    label="",
    show=true,
    aspect_ratio=1,
    c=:vik,
    clims=(ϕ0-ϕ0/5,ϕ0+ϕ0/5),
    markerstrokecolor=:auto,
    markerstrokewidth=0,
    ylim=(-10,20),
    )
    savefig(path_plot*"phi0.png")
end
@views function __plotStuff(mpD,type,ctr)
    xlab,ylab = L"$x-$direction",L"$z-$direction"
    gr(size=(2*250,2*125),legend=true,markersize=2.5,markershape=:circle,markerstrokewidth=0.75,)#markerstrokecolor=:match,)
    if type == "P"
        p = -(mpD.σ[1,:]+mpD.σ[2,:]+mpD.σ[3,:])/3/1e3
        scatter(mpD.x[:,1],mpD.x[:,2],zcolor=p,
            xlabel = xlab,
            ylabel = ylab,
            label=L"$p=-\dfrac{1}{3}\left(\sigma_{xx,p}+\sigma_{yy,p}+\sigma_{zz,p}\right)$",
            aspect_ratio=1,
            c=:viridis,
            ylim=(-10.0,20.0),
            title="Pressure",
            show=true,
            )  
    elseif type == "epII"
        scatter(mpD.x[:,1],mpD.x[:,2],zcolor=mpD.ϵpII,
            xlabel = xlab,
            ylabel = ylab,    
            label=L"$\epsilon_{\mathrm{II}}^{\mathrm{acc}}$",
            aspect_ratio=1,
            c=:viridis,
            clims=(0.0,2.0),
            ylim=(-10.0,20.0),
            title="Plastic strain",
            show=true,
            ) 
    elseif type == "du"
        scatter(mpD.x[:,1],mpD.x[:,2],zcolor=sqrt.(mpD.u[:,1].^2+mpD.u[:,2].^2),
            markershape=:circle,
            xlabel = xlab,
            ylabel = ylab,
            label=L"$\Delta u$",
            aspect_ratio=1,
            c=:viridis,
            ylim=(-10.0,20.0),
            title="Displacement",
            show=true,
            )
    else
        @error "field --"*string(type)*"-- not found as a valid arg. for plot"
        exit(1)
    end
    return ctr+=1
end