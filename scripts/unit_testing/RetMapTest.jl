# include dependencies
include("../../src/superInclude.jl")
# main program
@views function RetMapTest(L::Vector{Float64},nel::Int64,varPlot::String,cmType::String; kwargs...)
    ϕ∂ϕType,fwrkDeform,trsfrAp,isΔFbar,isGRF = getKwargs(kwargs)
    @info "** ϵp$(length(L))De v$(getVersion()): $(fwrkDeform) strain formulation **"
    @info "init..."
    # mesh setup
    meD     = meshSetup(nel,L,typeD)                                            # mesh geometry setup
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    E,ν     = 1.0e6,0.3                                                         # Young's mod. [Pa], Poisson's ratio [-]
    K,G,Del = D(E,ν,meD.nD)                                                     # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    ρ0      = 2700.0                                                            # density [kg/m^3]
    yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
    t,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mp setup
    mpD     = pointSetup(meD,L,c0,cr,ϕ0,ϕr,ρ0,isGRF,typeD)                      # material point geometry setup
    Hp      = -60.0e3*meD.h[1]                                                  # softening modulus
    # constitutive model param.
    cmParam = (E = E, ν = ν, Kc = K, Gc = G, Del = Del, Hp = Hp,)
    # plot & time stepping parameters
    tw,tC,it,ctr,ηmax,ηtot = 0.0,1.0,0,0,0,0    
    # action
    @info "mesh & mp feature(s):" nel=Tuple(meD.nel) nno=Tuple(meD.nno) nmp=mpD.nmp
    @info "launch $(ϕ∂ϕType) calculation cycle using $(nthreads()) thread(s)..."
    prog  = ProgressUnknown("working hard:", spinner=true,showspeed=true)

    for p ∈ 1:mpD.nmp
        y = mpD.x[p,2]
        P = (-ρ0*9.81*y)
        τ = cmParam.Gc
        mpD.σ[:,p] .= (τ.*[0.0,0.0,1.0]+P.*[1.0,1.0,0.0])
        mpD.σ[1,p] *= rand()
    end

    Pnmp = []
    Qnmp = Pnmp
    for p ∈ 1:mpD.nmp
        χ         = 3.0/2.0
        P,ξ,ξn,J2 = getCamCParam(mpD.σ[:,p],3)
        q         = sqrt(χ)*ξn
        push!(Pnmp,P)
        push!(Qnmp,q)
    end

    gr() # We will continue onward using the GR backend
    tit = ""
    plot(Pnmp,Qnmp, show=true, markershape=:circle,markersize=1.0, color = :blue, seriestype = :scatter, title = tit,xlabel=L"p/p_c",ylabel=L"q/p_c",)

    #η = camCRetMap!(mpD,cmParam,fwrkDeform) # Borja (1990); De Souza Neto (2008)
    #println(η)
    @info "Figs saved in" path_plot
    return msg("(✓) Done! exiting...")
end
# include("./scripts/unit_testing/RetMapTest.jl")
# e.g., L = [64.1584,12.80] or L = [64.1584,5.0,12.80]                                                                                        
# RetMapTest([64.1584,12.80],40,"P","DP";shpfun=:bsmpm,fwrk=:finite,trsf=:mUSL,vollock=true)
