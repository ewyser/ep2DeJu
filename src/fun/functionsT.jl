#----------------------------------------------------------------------------------------------------------
function meshSetup(nel,Lx,Lz,typeD)
    # geometry                                               
    L   = [Lx,ceil(Lz)]
    h   = [L[1]/nel,L[1]/nel]
    # mesh 
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
    # nodal quantities
    mn  = zeros(typeD,nno[3],1) 
    fen = zeros(typeD,nno[3],2) 
    fin = zeros(typeD,nno[3],2)
    fn  = zeros(typeD,nno[3],2)
    an  = zeros(typeD,nno[3],2)
    pn  = zeros(typeD,nno[3],2)
    vn  = zeros(typeD,nno[3],2)
    un  = zeros(typeD,nno[3],2)
    pel = zeros(typeD,nno[3],2)
    # mesh-to-node topology
    e2n = e2N(nno,nel,nn)
    # boundary conditions
    xB  = [minimum(xn)+2*h[1],maximum(xn)-2*h[1],0.0,Inf]                                    
    bcx = vcat(findall(x->x<=xB[1], xn),findall(x->x>=xB[2], xn))
    bcz = findall(x->x<=xB[3], zn)
    # push to struct
    meD = mesh(nel,nno,nn,L,h,xn,zn,mn,fen,fin,fn,an,pn,vn,un,pel,e2n,xB)
    bc  = boundary(bcx,bcz)
    return(meD,bc)
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function pointSetup(meD,ni,lz,coh0,cohr,phi0,phir,rho0,nstr,typeD)
    # mpm initialization
    xL  = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
    zL  = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:lz-0.5*meD.h[2]/ni
    npx = length(xL)
    npz = length(zL)
    xp  = ((xL'.*ones(npz,1  )      ))
    zp  = ((     ones(npx,1  )'.*zL ))



    xp  = vec(xp)
    zp  = vec(zp)
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
    coh  =  clt
    coh,phi  = RFS(xp[:,1],xp[:,2],coh0,cohr,phi0,phir)
    
    cohr =  ones(typeD,nmp,1).*cohr
    phi  =  ones(typeD,nmp,1).*phi0
    p    = findall(x->x<=2*wl, xp[:,2])
    phi[p] .= phir

    epII = zeros(typeD,nmp,1)
    J    = ones(typeD,nmp,1)
    # tensors
    dF   = zeros(typeD,nmp,4) 
    F    = repeat([1 0 0 1],nmp,1)
    b    = zeros(typeD,nmp,4)
    bT   = zeros(typeD,nmp,4)
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
    B    = zeros(typeD,nstr,nn.*2,nmp)
    # connectivity
    p2e  = zeros(UInt64,nmp,1)
    p2n  = zeros(UInt64,nmp,nn)
    # push to struct
    mpD  = point(nmp,l0,l,v0,v,m,xp,up,vp,pp,coh,cohr,phi,epII,J,
                 dF,F,b,bT,e,ome,s,τ,dev,ep,ϕ,∂ϕx,∂ϕz,B,p2e,p2n)
    return(mpD)
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function e2N(nno,nel,nn)
	gnum = reverse(reshape(1:(nno[3]),nno[2],nno[1]),dims=1)
	e2n  = zeros(nel[3],nn)
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
	return(convert(Array{Int64},e2n))
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function save2txt(meD,mpD,bc)
    xn = [meD.xn meD.zn]
    writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/xn.txt",vec(xn))
    writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/e2n.txt",vec(meD.e2n.-1))
    
    writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/cohp.txt",vec(mpD.coh))
    writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/phip.txt",vec(mpD.phi))

    p = [mpD.nmp meD.nn meD.nno[3] meD.h[1] meD.h[2] minimum(meD.xn) minimum(meD.zn) meD.nno[1] meD.nno[2]]
    writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/param.txt",vec(p))    

    writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/mp.txt" ,vec(mpD.mp) )    
    writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/xp.txt" ,vec(mpD.xp))    
    writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/vol.txt",vec(mpD.v) )    
    writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/lp.txt" ,vec(mpD.l)) 


    bcx = bc.x.+(0*meD.nno[3])
    bcz = bc.z.+(1*meD.nno[3])
    BC  = ones(Int64,meD.nno[3]*2,1)
    BC[vcat(bcx,bcz)].= 0
    writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/bcs.txt" ,vec(BC)) 

    return("done")
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
    return(Kc,Gc,D)
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
    return(Δt)
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function get_g(tw::Float64,tg::Float64)
    g = 0.0
    if(tw<=tg)
        g = 9.81*tw/tg
    else
        g = 9.81
    end
    return(g)
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
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
@views function plot_Δu(xp,up)
            Δu = sqrt.(up[:,1].^2+up[:,2].^2)
            gr(size=(2*250,2*125),legend=true,markersize=2.5)
            scatter(xp[:,1],xp[:,2],zcolor=Δu,
                markershape=:circle,
                label="",
                show=true,
                aspect_ratio=1,
                c=:viridis,
                ylim=(-10.0,20.0),
                )     
end
@views function plot_P(xp,σ)
            p = -(σ[1,:]+σ[2,:]+σ[3,:])/3/1e3
            gr(size=(2*250,2*125),legend=true,markersize=2.5)
            scatter(xp[:,1],xp[:,2],zcolor=p,
                markershape=:circle,
                label="",
                show=true,
                aspect_ratio=1,
                c=:viridis,
                ylim=(-10.0,20.0),
                )     
end
@views function plot_Δϵp(xp,epII)
            gr(size=(2*250,2*125),legend=true,markersize=2.5)
            scatter(xp[:,1],xp[:,2],zcolor=epII,
                markershape=:circle,
                label="",
                show=true,
                aspect_ratio=1,
                c=:viridis,
                clims=(0.0,2.0),
                ylim=(-10.0,20.0),
                )     
end
