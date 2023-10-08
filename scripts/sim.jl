# julia -i -O3 -t auto --check-bounds=no --project=.
# include("./scripts/sim.jl")
# ϵp2De(80,"P","mohr",true)

# include dependencies
include("../src/superInclude.jl")
# arithmetic precision (double=Float64 or single=Float32)
typeD = Float64  
# relative path for figs & data
path_plot = "./out/"
if isdir(path_plot)==false mkdir(path_plot) end

@views function ϵp2De(nel::Int64,varPlot::String,cmType::String,isΔFbar::Bool)
    @info "** ϵp2-3De v1.0: finite strain formulation **"
    # non-dimensional constant                                                   
    ni,ndim,nstr = 2,2,4                                                        # number of material point along 1d, number of stresses
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    K,G,Del = D(1.0e6,0.3)                                                      # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    ρ0      = 2700.0                                                            # density [kg/m^3]
    yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*pi/180,7.5*pi/180,0.0                                        # friction angle [Rad], dilation angle [Rad]                                                              
    t,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    lx,lz   = 64.1584,12.80                                                     # domain geometry
    meD     = meshSetup(nel,(lx,lz),ndim,typeD)                                 # mesh geometry setup
    mpD     = pointSetup(meD,ni,lz,c0,cr,ϕ0,ϕr,ρ0,nstr,typeD)                   # material point geometry setup
    Hp      = -60.0e3*meD.h[1]                                                  # softening modulus
    # constitutive model param.
    cmParam = (Kc = K, Del = Del, Hp = Hp,)
    @info "mesh & mp feature(s):" nel=Int64(meD.nel[end]) nno=meD.nno[end] nmp=mpD.nmp
    # plot & time stepping parameters
    tw,tC,it,ctr,toc,flag,ηmax,ηtot = 0.0,1.0/1.0,0,0,0.0,0,0,0    
    # action
    @info "launch bsmpm calculation cycle..."
    prog  = ProgressUnknown("working hard:", spinner=true,showspeed=true)
    while tw<=t
        # plot/save
        if tw >= ctr*tC
            ctr = __plotStuff(mpD,varPlot,ctr)
        end
        # set clock on/off
        tic = time_ns()
        # adaptative Δt & linear increase in gravity
        Δt,g  = get_Δt(mpD.v,meD.h,yd),get_g(tw,tg,meD.nD)
        # bsmpm cycle
        ϕ∂ϕ!(mpD,meD)
        mapsto!(mpD,meD,g,Δt,"p->N")                  
        solve!(meD,Δt)
        mapsto!(mpD,meD,g,Δt,"p<-N")
        ηmax = elastoplast!(mpD,meD,cmParam,cmType,isΔFbar,tw>te)
        if tw>te && flag == 0
            plot_coh(mpD.x,mpD.coh,mpD.phi,ϕ0)
            flag+=1
        end
        # update sim time
        tw,it,toc,ηtot = tw+Δt,it+1,((time_ns()-tic)),max(ηmax,ηtot)
        # update progress bas
        next!(prog;showvalues = [("[nel,np]",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),("iteration(s)",it),("ηmax,ηtot",(ηmax,ηtot)),("(✗) t/T",round(tw/t,digits=2))])
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = [("[nel,np]",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),("iteration(s)",it),("ηmax,ηtot",(ηmax,ηtot)),("(✓) t/T",1.0)])
    savefig(path_plot*varPlot*"_"*cmType*"_vollock_"*string(isΔFbar)*".png")
    @info "Figs saved in" path_plot
    println("[=> done! exiting...")
    return nothing
end