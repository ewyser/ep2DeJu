# include("./scripts/unit_testing/compactionTest.jl")
# include dependencies
include("../../src/superInclude.jl")
using Test
# main program
function meshGeom(L,nel)
    nD = length(L)
    h   = [L[1],L[2]/nel]
    return L,h,nD
end
function meshCoord(nD,L,h)
    xn  = collect((0.0-2*h[1]):h[1]:(L[1]+2.0*h[1])) 
    zn  = reverse(collect((0.0-2*h[2]):h[2]:(L[2]+2.0*h[2])))
    nno = [length(xn),length(zn),length(xn)*length(zn)] 
    nel = [nno[1]-1,nno[2]-1,(nno[1]-1)*(nno[2]-1)]
    nn  = 16
    xn  = (xn'.*ones(typeD,nno[2],1     ))     
    zn  = (     ones(typeD,nno[1],1     )'.*zn)
    x   = hcat(vec(xn),vec(zn))
    return x,nn,nel,nno
end
function meshBCs(xn,h,nno,nD)
    if nD == 2
        xB  = [minimum(xn[:,1])+2*h[1],maximum(xn[:,1])-2*h[1],0.0,Inf]                                    
        bcx = vcat(findall(x->x<=xB[1], xn[:,1]),findall(x->x>=xB[2], xn[:,1]))
        bcz = findall(x->x<=xB[3], xn[:,2])
        bcX = ones(Int64,nno[nD+1],1)
        bcX[bcx] .= 0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0
        #bcX[bcz] .= 0
        bc   = hcat(bcX,bcZ)
    elseif nD == 3
        xB  = [minimum(xn[:,1])+2*h[1],maximum(xn[:,1])-2*h[1],minimum(xn[:,2])+2*h[1],maximum(xn[:,2])-2*h[1],0.0,Inf]                                    
        bcx = vcat(findall(x->x<=xB[1], xn[:,1]),findall(x->x>=xB[2], xn[:,1]))
        bcy = vcat(findall(x->x<=xB[3], xn[:,2]),findall(x->x>=xB[4], xn[:,2]))
        bcz = findall(x->x<=xB[5], xn[:,3])
        bcX = ones(Int64,nno[nD+1],1)
        bcX[bcx] .= 0
        bcY = ones(nno[nD+1],1)
        bcY[bcy] .= 0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0
        bc   = hcat(bcX,bcY,bcZ)
    end
    return bc,xB
end
function meshSetup(nel,L,typeD)
    # geometry                                               
    L,h,nD       = meshGeom(L,nel)
    # mesh 
    x,nn,nel,nno = meshCoord(nD,L,h)
    # boundary conditions
    bc,xB        = meshBCs(x,h,nno,nD)
    # constructor
    meD = (
        nD   = nD,
        nel  = nel,
        nno  = nno,
        nn   = nn,
        L    = L,
        h    = h,
        minC = minimum(x,dims=2),
        # nodal quantities
        xn   = x,
        mn   = zeros(typeD,nno[nD+1]             ), # lumped mass vector
        Mn   = zeros(typeD,nno[nD+1],nno[nD+1]   ), # consistent mass matrix
        fext = zeros(typeD,nno[nD+1],nD          ), 
        fint = zeros(typeD,nno[nD+1],nD          ),
        oobf = zeros(typeD,nno[nD+1],nD          ),
        Dn   = zeros(typeD,nno[nD+1],nD          ),
        fn   = zeros(typeD,nno[nD+1],nD          ),
        an   = zeros(typeD,nno[nD+1],nD          ),
        pn   = zeros(typeD,nno[nD+1],nD          ),
        vn   = zeros(typeD,nno[nD+1],nD          ),
        Δun  = zeros(typeD,nno[nD+1],nD          ),
        ΔJn  = zeros(typeD,nno[nD+1],nD          ),
        bn   = zeros(typeD,nD       ,nD,nno[nD+1]),
        # mesh-to-node topology
        e2n  = e2N(nD,nno,nel,nn),
        xB   = xB,
        # mesh boundary conditions
        bc   = bc,
    )
    return meD
end
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
@views function solve!(meD,Δt)
    # viscous damping
    η   = 0.05
    # initialize
    meD.fn .= 0.0
    meD.an .= 0.0
    meD.vn .= 0.0
    # solve momentum equation on the mesh
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.mn[n]>0.0 
            m           = (1.0/meD.mn[n]).*meD.bc[n,:]                   #(2,)
            meD.Dn[n,:].= η.*norm(meD.oobf[n,:]).*sign.(meD.pn[n,:].*m)  #(2,)
            meD.fn[n,:].= meD.oobf[n,:].-meD.Dn[n,:]                     #(2,)
            meD.an[n,:].= meD.fn[n,:].*m                                 #(2,)
            meD.vn[n,:].= (meD.pn[n,:].+Δt.*meD.fn[n,:]).*m              #(2,)
        end
    end
    return nothing
end
@views function getVals(it,cmpl,symb)
    # completion [%]
    cmpl = round(100.0*cmpl,digits=1)
    # save vals
    vals = [("iteration(s)",it),
            (symb*" t/T",cmpl)]
    return vals
end

@views function compactTest(nel,varPlot,ν,E,ρ0,l0; kwargs...)
    cmType = "MC"
    ϕ∂ϕType,fwrkDeform,trsfrAp,isΔFbar = getKwargs(kwargs)
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    K,G,Del = D(E,ν)                                                            # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
    tg      = ceil((1.0/yd)*(2.0*l0)*40.0)
    t,te    = 1.25*tg,1.25*tg
    # mesh & mp setup
    L       = [l0/nel,l0]                                                       # domain geometry
    meD     = meshSetup(nel,L,typeD)                                            # mesh geometry setup
    mpD     = pointSetup(meD,L,c0,cr,ϕ0,ϕr,ρ0,typeD)                            # material point geometry setup
    z0      = copy(mpD.x[:,end])
    Hp      = -60.0e3*meD.h[1]                                                  # softening modulus
    # constitutive model param.
    cmParam = (Kc = K, Gc = G, Del = Del, Hp = Hp,)
    @info "mesh & mp feature(s):" nel=Int64(meD.nel[2]-4)
    # plot & time stepping parameters
    tw,tC,it,ctr,ηmax,ηtot = 0.0,1.0,0,0,0,0    
    # action
    g = [0.0,0.0]
    prog  = ProgressUnknown("working hard:", spinner=true,showspeed=true)
    while tw<=t
        # set clock on/off
        tic = time_ns()
        # adaptative Δt & linear increase in gravity
        Δt,g  = get_Δt(mpD.v,meD.h,yd),get_g(tw,tg,meD.nD)
        # bsmpm cycle
        ϕ∂ϕ!(mpD,meD,ϕ∂ϕType)
        mapsto!(mpD,meD,g,Δt,trsfrAp,"p->n")                  
        solve!(meD,Δt)
        mapsto!(mpD,meD,g,Δt,trsfrAp,"p<-n")
        ηmax = elastoplast!(mpD,meD,cmParam,cmType,Δt,ϕ∂ϕType,isΔFbar,fwrkDeform,tw>te)
        # update sim time
        tw,it,toc,ηtot = tw+Δt,it+1,((time_ns()-tic)),max(ηmax,ηtot)
        next!(prog;showvalues = getVals(it,tw/t,"(✗)"))
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = getVals(it,1.0,"(✓)"))
    ctr     = plotStuff(mpD,tw,varPlot,ctr,L"$g = $"*string(round(g[end],digits=2))*L" [m.s$^{-2}$]")
    savefig(path_plot*"$(varPlot)_compaction_self_weight_test_$(ϕ∂ϕType)_$(fwrkDeform).png")
    # analytics
    xN,yN = abs.(mpD.σ[2,:]),z0
    xA,yA = abs.(ρ0.*g[end].*(l0.-z0)),z0
    err   = sum(sqrt.((xN.-xA).^2).*mpD.V0)/(abs(g[end])*ρ0*l0*sum(mpD.V0))
    return (xN,yN,xA,yA),meD.h,err
end
@views function compacTest()
    ϕ∂ϕType    = :gimpm
    fwrkDeform = :finite
    @info "** ϵp2De v$(getVersion()): compaction of a two-dimensional column under self weight **"
    store,H,error = [],[],[]
    try
        @testset "convergence using $(ϕ∂ϕType), $(fwrkDeform) deformation" begin
            # geometry
            n         = [0,1,2,3,4,5]
            nel       = 2.0.^n
            # initial parameters 
            l0,ν,E,ρ0 = 50.0,0.0,1.0e4,80.0
            # init error
            ϵ         = 1.0
            for (it,nel) in enumerate(nel)
                #action
                DAT,h,err = compactTest(nel,"P",ν,E,ρ0,l0;shpfun=ϕ∂ϕType,fwrk=fwrkDeform,vollock=true)
                push!(store,DAT )
                push!(H ,h[end])
                push!(error,err)
                # test
                @test (err < ϵ )
                ϵ = err
            end 
        end
        gr(size=(2.0*250,2*125),legend=false,markersize=2.25,markerstrokecolor=:auto)
        p1 = plot(1.0./H,error,seriestype=:scatter, label="convergence",xlabel=L"$1/h$ [m$^{-1}$]",ylabel="error",xaxis=:log,yaxis=:log) 
        display(plot(p1; layout=(1,1), size=(450,250)))
        savefig(path_plot*"convergence_pass_compacTest_$(ϕ∂ϕType)_$(fwrkDeform).png")
    catch
        gr(size=(2.0*250,2*125),legend=false,markersize=2.25,markerstrokecolor=:auto)
        p1 = plot(1.0./H,error,seriestype=:scatter, label="convergence",xlabel=L"$1/h$ [m$^{-1}$]",ylabel="error",xaxis=:log,yaxis=:log) 
        display(plot(p1; layout=(1,1), size=(450,250)))
        savefig(path_plot*"convergence_fail_compacTest_$(ϕ∂ϕType)_$(fwrkDeform).png")
    end
    xN,yN,xA,yA = store[end]
    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    p1 = plot(xN.*1e-3,yN,seriestype=:scatter, label="numerical approximation")
    p1 = plot!(xA.*1e-3,yA,label="analytical solution",xlabel=L"$\sigma_{yy}$ [kPa]",ylabel=L"$y-$position [m]") 
    display(plot(p1; layout=(1,1), size=(450,250)))
    savefig(path_plot*"numericVsAnalytic_compacTest_$(ϕ∂ϕType)_$(fwrkDeform).png")
end
compacTest()








































#=
xN,yN,xA,yA,err,h = store[k]
nely = nel[k]
gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
p1 = plot(xN.*1e-3,yN,seriestype=:scatter, label="numerical approximation")
p1 = plot!(xA.*1e-3,yA,label="analytical solution",xlabel=L"$\sigma_{yy}$ [kPa]",ylabel=L"$y-$position [m]") 
display(plot(p1; layout=(1,1), size=(450,250)))
savefig(path_plot*"numericVsAnalytic_compaction_self_weight_test_nel_$(nely)_$(ϕ∂ϕType)_$(fwrkDeform).png")
=#