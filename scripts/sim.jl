typeD = Float64  # Precision (double=Float64 or single=Float32)

path_plot = "./out/"
if isdir(path_plot)==false
    mkdir(path_plot)    
end

include("../src/superInclude.jl")



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

    tw = 0.0
    it = 0
    itps = 0.0
    nout = 50
    wct  = 0.0
    flag = 0

    @printf("\no---------------------------------------------o")
    @printf("\n|             ** ϵp2-3De v1.0 **              |")
    @printf("\n|      -- finite strain formulation --        |")
    @printf("\no---------------------------------------------o")
    @printf("\n nel = [%d,%d,%d]",meD.nel[1],meD.nel[2],meD.nel[3])
    @printf("\n nno = [%d,%d,%d]",meD.nno[1],meD.nno[2],meD.nno[3])
    @printf("\n nmp = %d",mpD.nmp)
    @printf("\no---------------------------------------------o\n") 
    #prog = Progress(ceil(Int,t/Δt),dt=0.5,barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),barlen=16,showspeed=false)
    # action
    println("[=> action!")
    prog  = ProgressUnknown("working hard:", spinner=true,showspeed=true)
    while tw<t
        t0  = Base.time()
        Δt  = get_Δt(mpD.vp,meD.h,yd)
        g   = get_g(tw,tg)
        topol!(mpD.p2e,mpD.p2n,meD.e2n,meD.xn,meD.zn,mpD.xp,meD.h,meD.nel,mpD.nmp,meD.nn)
        ϕ∂ϕ!(mpD.B,mpD.ϕ,mpD.∂ϕx,mpD.∂ϕz,mpD.xp,meD.xn,meD.zn,mpD.p2n,meD.h,meD.xB,meD.nn,mpD.nmp)
        accum!(meD.mn,meD.pn,meD.fen,meD.fin,mpD.σ,mpD.τ,mpD.J,mpD.vp,mpD.v,mpD.mp,mpD.ϕ,mpD.B,mpD.p2n,g,mpD.nmp,meD.nn)
        solve!(meD.fn,meD.an,meD.pn,meD.vn,meD.mn,meD.fen,meD.fin,bc.x,bc.z,meD.nno,Δt)
        flip!(mpD.vp,mpD.xp,mpD.ϕ,meD.an,meD.vn,mpD.p2n,mpD.nmp,Δt) 
        DMBC!(mpD.up,meD.pn,meD.un,meD.mn,mpD.ϕ,mpD.vp,mpD.mp,mpD.p2n,bc.x,bc.z,mpD.nmp,meD.nn,meD.nno,Δt)   # need to be improved
        elast!(mpD.τ,mpD.ϵ,mpD.J,mpD.v,mpD.v0,mpD.l,mpD.l0,mpD.F,meD.un,mpD.∂ϕx,mpD.∂ϕz,mpD.p2n,mpD.nmp,Del) # need to be improved
        if(tw>te)
            #plast!(mpD.τ,mpD.ϵ,mpD.epII,mpD.coh,mpD.phi,mpD.nmp,Del,Hp,cr)
            CPAplast!(mpD.τ,mpD.ϵ,mpD.epII,mpD.coh,mpD.phi,mpD.nmp,Del,Hp,cr)
            if(flag==0)
                plot_coh(mpD.xp,mpD.coh,mpD.phi,ϕ0)
                flag+=1
            end
        end
        tw += Δt
        it += 1
        Δτ  = (Base.time()-t0)
        itps= 1/Δτ
        if it>1
            wct+= Δτ
        end
        if(mod(it,nout)==0)
            plot_Δϵp(mpD.xp,mpD.epII)       
        end 
        next!(prog;showvalues = [("[nel,np]",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),("iteration(s)",it),("(✗) t/T",round(tw/t,digits=2))])
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = [("[nel,np]",(meD.nel[1]*meD.nel[2],mpD.nmp)),("iteration(s)",it),("(✓) t/T",1.0)])
    savefig(path_plot*"plot.png")
    @info "Figs saved in" path_plot
    println("[=> done! exiting...")
end
main()



# https://techytok.com/lesson-parallel-computing/