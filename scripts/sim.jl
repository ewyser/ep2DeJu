# julia -i -O3 -t auto --check-bounds=no --project=.
# include("./scripts/sim.jl")

# include dependencies
include("../src/superInclude.jl")
# arithmetic precision (double=Float64 or single=Float32)
typeD = Float64  
# relative path for figs & data
path_plot = "./out/"
if isdir(path_plot)==false
    mkdir(path_plot)    
end

@views function main()
    # ---------------------------------------------------------------------------
    nel = 80
    # ---------------------------------------------------------------------------
    # non-dimensional constant 
    # ---------------------------------------------------------------------------
    ν       = 0.3                                                               # poisson's ratio                                                        
    ni      = 2                                                                 # number of material point along 1d
    nstr    = 4                                                                 # number of stresses
    # ---------------------------------------------------------------------------
    # independant physical constant
    # ---------------------------------------------------------------------------
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]
    E       = 1.0e6                                                             # young's modulus             [Pa]
    K,G,Del = D(E,ν)                                                            # elastic matrix    
    ρ0      = 2700.0                                                            # density                     [kg/m^3]
    yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed          [m/s]
    c0      = 20.0e3                                                            # cohesion                    [Pa]
    ϕ0      = 20.0*pi/180                                                       # friction angle              [Rad]
    ψ0      = 0.0                                                               # dilantancy angle
    H       = -60.0e3                                                           # softening modulus           [Pa]
    cr      =   4.0e3                                                           # residual cohesion           [Pa]
    ϕr      = 7.5*pi/180                                                        # residual friction angle     [Rad]
    t       = 15.0                                                              # simulation time             [s]
    te      = 10.0                                                              # elastic loading             [s]
    tg      = te/1.5                                                            # gravity increase 
    # ---------------------------------------------------------------------------
    # mesh & mp setup
    # ---------------------------------------------------------------------------
    lx      = 64.1584                                                           # domain length along the x-direction
    lz      = 12.80                                                             # domain length along the z-direction
    meD,bc  = meshSetup(nel,lx,lz,typeD)                                        # mesh geometry setup
    mpD     = pointSetup(meD,ni,lz,c0,cr,ϕ0,ϕr,ρ0,nstr,typeD)                   # material point geometry setup
    Hp      = H*meD.h[1]                                                        # softening modulus
    # ---------------------------------------------------------------------------
    # display parameters & runtime
    # ---------------------------------------------------------------------------                                                            
    Δt      = 0.5*meD.h[1]/yd                                                   # unconditionally stable timestep
    nit     = ceil(t/Δt)                                                        # maximum number of interation
    nf      = max(2,ceil(round(1/Δt)/25))                                       # number of frame interval
    # runtime parameters
    it      = 1                                                                 # initialize iteration
    tw      = 0.0                                                               # initialize time while statement
    
    #char    = save2txt(meD,mpD,bc)
    #p       = [g;ρ0;ψ0;ν;E;Kc;Gc;cr;Hp;t;te;tg]
    #writedlm("/Users/manuwyser/Dropbox/PhD_Thesis/git_local/work_mpm/C_code_2D/scripts/setting_Exp2b/phys.txt" ,p)

    

    bcx = ones(meD.nno[3],1)
    bcx[bc.x] .= 0
    bc.x = bcx
    bcz = ones(meD.nno[3],1)
    bcz[bc.z] .= 0
    bc.z = bcz



    
    # plot & time stepping parameters
    tw,it,ctr,toc,flag = 0.0,0,0,0.0,0    
    tC = 1.0/1.0

    @info "** ϵp2-3De v1.0: finite strain formulation **"
    @info "mesh & mp feature(s):" nel=Int64(meD.nel[end]) nno=meD.nno[end] nmp=mpD.nmp
    println("[=> action!")
    prog  = ProgressUnknown("working hard:", spinner=true,showspeed=true)
    while tw<=t
        # plot/save
        if tw >= ctr*tC
            plot_Δϵp(mpD.xp,mpD.epII)
            #plot_v(mpD.xp,mpD.v)
            ctr+=1
        end
        # set clock in/off
        tic = time_ns()
        # adaptative Δt
        Δt  = get_Δt(mpD.vp,meD.h,yd)
        g   = get_g(tw,tg)
        topol!(mpD.p2e,mpD.p2n,meD.e2n,meD.xn,meD.zn,mpD.xp,meD.h,meD.nel,mpD.nmp,meD.nn)
        ϕ∂ϕ!(mpD.B,mpD.ϕ∂ϕ,mpD.xp,meD.xn,meD.zn,mpD.p2n,meD.h,meD.xB,meD.nn,mpD.nmp)
        accum!(meD.mn,meD.pn,meD.fen,meD.fin,mpD.σ,mpD.τ,mpD.J,mpD.vp,mpD.v,mpD.mp,mpD.ϕ∂ϕ,mpD.B,mpD.p2n,g,mpD.nmp,meD.nn)
        solve!(meD.fn,meD.an,meD.pn,meD.vn,meD.mn,meD.fen,meD.fin,bc.x,bc.z,meD.nno,Δt)
        flip!(mpD.vp,mpD.xp,mpD.ϕ∂ϕ,meD.an,meD.vn,mpD.p2n,mpD.nmp,Δt) 
        DMBC!(mpD.up,meD.pn,meD.un,meD.mn,mpD.ϕ∂ϕ,mpD.vp,mpD.mp,mpD.p2n,bc.x,bc.z,mpD.nmp,meD.nn,meD.nno,Δt)   # need to be improved
        deform!(mpD.ΔJ,mpD.J,mpD.Jbar,meD.Jn,mpD.mp,mpD.v,meD.Vn,mpD.ΔF,mpD.ΔFbar,mpD.F,mpD.Fbar,meD.un,mpD.ϕ∂ϕ,mpD.p2n,mpD.nmp)
        elast!(mpD.τ,mpD.ϵ,mpD.ΔJ,mpD.J,mpD.v,mpD.v0,mpD.l,mpD.l0,mpD.ΔF,mpD.ΔFbar,mpD.F,mpD.nmp,Del) # need to be improved
        if tw>te
            #plast!(mpD.τ,mpD.ϵ,mpD.epII,mpD.coh,mpD.phi,mpD.nmp,Del,Hp,cr)
            CPAplast!(mpD.τ,mpD.ϵ,mpD.epII,mpD.coh,mpD.phi,mpD.nmp,Del,Hp,cr)
            if flag==0 
                plot_coh(mpD.xp,mpD.coh,mpD.phi,ϕ0)
                flag+=1
            end
        end
        # update sim time
        tw,it,toc = tw+Δt,it+1,((time_ns()-tic))
        # update progress bas
        next!(prog;showvalues = [("[nel,np]",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),("iteration(s)",it),("(✗) t/T",round(tw/t,digits=2))])
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = [("[nel,np]",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),("iteration(s)",it),("(✓) t/T",1.0)])
    savefig(path_plot*"plot.png")
    @info "Figs saved in" path_plot
    println("[=> done! exiting...")
end
main()



# https://techytok.com/lesson-parallel-computing/