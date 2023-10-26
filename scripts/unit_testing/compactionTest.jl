# include("./scripts/unit_testing/compactionTest.jl")
# include dependencies
include("../../src/superInclude.jl")
# main program
function materialGeomCompact(meD,lz,wl,coh0,cohr,ni)
    xL          = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
    zL          = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:lz-0.5*meD.h[2]/ni
    npx,npz     = length(xL),length(zL)
    xp,zp       = ((xL'.*ones(npz,1  )      )),((     ones(npx,1  )'.*zL )) 
    xp,zp       = vec(xp),vec(zp)
    id          = shuffle(collect(1:size(xp,1)))
    return hcat(xp[id,:],zp[id,:])
end
function pointSetup(meD,L,coh0,cohr,phi0,phir,rho0,typeD)
    # non-dimensional constant                                                   
    ni,nstr = 2,4                                                               # number of material point along 1d, number of stresses
    # material geometry
    lz     = L[end]
    wl     = 0.15*lz
    xp     = materialGeomCompact(meD,lz,wl,coh0,cohr,ni)
    # scalars & vectors
    nmp    = size(xp,1)
    l0,l   = ones(typeD,nmp,2).*0.5.*(meD.h[1]./ni),ones(typeD,nmp,2).*0.5.*(meD.h[1]./ni)
    v0,v   = ones(typeD,nmp  ).*(2.0.*l0[:,1].*2.0.*l0[:,2]),ones(typeD,nmp  ).*(2.0.*l[:,1].*2.0.*l[:,2])
    m      = rho0.*v0
    coh    = ones(typeD,nmp ).*coh0
    cohr   = ones(typeD,nmp).*cohr
    phi    = ones(typeD,nmp).*phi0
    phi[xp[:,2].<=2*wl] .= phir
    # constructor
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
        ϵpII = zeros(typeD,nmp),
        ϵpV  = zeros(typeD,nmp), 
        ΔJ   = ones(typeD,nmp),
        J    = ones(typeD,nmp),
        # tensor in matrix notation
        I    = Matrix(1.0I,meD.nD,meD.nD    ),
        ΔF   = zeros(typeD,meD.nD,meD.nD,nmp),
        F    = repeat(Matrix(1.0I,meD.nD,meD.nD),1,1,nmp),
        ∇v   = zeros(typeD,meD.nD,meD.nD,nmp),
        ϵ    = zeros(typeD,meD.nD,meD.nD,nmp),
        b    = repeat(Matrix(1.0I,meD.nD,meD.nD),1,1,nmp),
        # tensor in voigt notation
        ω    = zeros(typeD,nmp),
        σR   = zeros(typeD,nstr,nmp),
        σ    = zeros(typeD,nstr,nmp),
        τ    = zeros(typeD,nstr,nmp),
        dev  = zeros(typeD,nstr,nmp),
        ep   = zeros(typeD,nstr,nmp),
        # additional quantities
        ϕ∂ϕ  = zeros(typeD,meD.nn,nmp ,meD.nD+1   ),
        B    = zeros(typeD,meD.nn.*meD.nD,nstr,nmp),
        # connectivity
        p2e  = zeros(Int64,nmp),
        p2n  = zeros(Int64,meD.nn,nmp),
    )
    return mpD 
end
@views function plotStuff(mpD,t,type,ctr,title)
    xlab,ylab = L"$x-$direction",L"$z-$direction"
    gr(size=(2*250,2*125),legend=true,markersize=2.5,markershape=:circle,markerstrokewidth=0.0,markerstrokecolor=:match,)
    temp = title
    if type == "P"
        p = -mpD.σ[2,:]/1e3
        scatter(mpD.x[:,1],mpD.x[:,2],zcolor=p,
            xlabel = xlab,
            ylabel = ylab,
            label=L"$\sigma_{zz}$",
            aspect_ratio=1,
            c=:viridis,
            ylim=(0.0,50.0),
            title=temp,
            show=true,
            )  
    end
    return ctr+=1
end

@views function compactTest(nel,varPlot,ν,E,ρ0,l0; kwargs...)
    cmType = "MC"
    ϕ∂ϕType,fwrkDeform,trsfrAp,isΔFbar = getKwargs(kwargs)
    @info "** ϵp2De v$(getVersion()): compaction of a two-dimensional column under self weight **"
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    K,G,Del = D(E,ν)                                                      # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
    t,te,tg = 10.0,10.0,10.0                                                    # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    L       = [10.0,l0]                                                        # domain geometry
    meD     = meshSetup(nel,L,typeD)                                            # mesh geometry setup
    mpD     = pointSetup(meD,L,c0,cr,ϕ0,ϕr,ρ0,typeD)                            # material point geometry setup
    Hp      = -60.0e3*meD.h[1]                                                  # softening modulus
    # constitutive model param.
    cmParam = (Kc = K, Gc = G, Del = Del, Hp = Hp,)
    @info "mesh & mp feature(s):" dim=meD.nD nel=Int64(meD.nel[end]) nno=meD.nno[end] nmp=mpD.nmp
    # plot & time stepping parameters
    tw,tC,it,ctr,ηmax,ηtot = 0.0,1.0,0,0,0,0    
    # action
    g = [0.0,0.0]
    @info "launch $(ϕ∂ϕType) calculation cycle..."
    prog  = ProgressUnknown("working hard:", spinner=true,showspeed=true)
    while tw<=t
        # plot/save
        if tw >= ctr*tC ctr = plotStuff(mpD,tw,varPlot,ctr,L"$g = $"*string(round(g[end],digits=2))*L" [m.s$^{-2}$]") end
        # set clock on/off
        tic = time_ns()
        # adaptative Δt & linear increase in gravity
        Δt,g  = get_Δt(mpD.v,meD.h,yd),get_g(tw,tg,meD.nD)
        # bsmpm cycle
        ϕ∂ϕ!(mpD,meD,ϕ∂ϕType)
        mapsto!(mpD,meD,g,Δt,trsfrAp,"p->n")                  
        solve!(meD,Δt)
        mapsto!(mpD,meD,g,Δt,trsfrAp,"p<-n")
        ηmax = elastoplast!(mpD,meD,cmParam,cmType,Δt,isΔFbar,fwrkDeform,tw>te)
        # update sim time
        tw,it,toc,ηtot = tw+Δt,it+1,((time_ns()-tic)),max(ηmax,ηtot)
        # update progress bas
        next!(prog;showvalues = getVals(meD,mpD,it,ηmax,ηtot,tw/t,"(✗)"))
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = getVals(meD,mpD,it,ηmax,ηtot,1.0,"(✓)"))
    ctr     = plotStuff(mpD,tw,varPlot,ctr,L"$g = $"*string(round(g[end],digits=2))*L" [m.s$^{-2}$]")
    sleep(2.5)
    savefig(path_plot*"$(varPlot)_compaction_self_weight_test.png")
    @info "Figs saved in" path_plot
    return msg("(✓) Done! exiting...")
end
# initial parameters
nel,l0 = 5,50.0
ν,E,ρ0 = 0.3,1.0e6,80.0
#action
compactTest(nel,"P",ν,E,ρ0,l0;shpfun=:bsmpm,fwrk=:finite,vollock=true)
