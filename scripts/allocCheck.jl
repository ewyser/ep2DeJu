# include dependencies
include("../src/superInclude.jl")
using BenchmarkTools

@warn "validation/test"
@views function allocCheck(nel::Int64,varPlot::String,cmType::String; kwargs...)
    ϕ∂ϕType,fwrkDeform,isΔFbar = getKwargs(kwargs)
    @info "** ϵp2-3De v$(getVersion()): $(fwrkDeform) strain formulation **"
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    K,G,Del = D(1.0e6,0.3)                                                      # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    ρ0      = 2700.0                                                            # density [kg/m^3]
    yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
    t,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    L       = [64.1584,12.80]                                                   # domain geometry
    meD     = meshSetup(nel,L,typeD)                                            # mesh geometry setup
    mpD     = pointSetup(meD,L,c0,cr,ϕ0,ϕr,ρ0,typeD)                            # material point geometry setup
    Hp      = -60.0e3*meD.h[1]                                                  # softening modulus
    # constitutive model param.
    cmParam = (Kc = K, Gc = G, Del = Del, Hp = Hp,)
    @info "mesh & mp feature(s):" nel=Int64(meD.nel[end]) nno=meD.nno[end] nmp=mpD.nmp
    # plot & time stepping parameters
    tw,tC,it,ctr,toc,flag,ηmax,ηtot = 0.0,1.0/1.0,0,0,0.0,0,0,0    
    # action
    
    @info "launch ϕ∂ϕ!()"
    @btime ϕ∂ϕ!($mpD,$meD,$ϕ∂ϕType) 
    @info "launch mapsto!(p->n)"
    @btime mapsto!($mpD,$meD,vec([0.0,0.0,9.81]),0.1,"p->n")   
    @info "launch solve!()"
    @btime solve!($meD,0.1)
    @info "launch mapsto!(p<-n)"
    @btime mapsto!($mpD,$meD,vec([0.0,0.0,9.81]),0.1,"p<-n")
    @info "launch elastoplast!(), elastic"
    @btime ηmax = elastoplast!($mpD,$meD,$cmParam,$cmType,$isΔFbar,$fwrkDeform,false)
    @info "launch elastoplast!(), elastoplastic"
    @btime ηmax = elastoplast!($mpD,$meD,$cmParam,$cmType,$isΔFbar,$fwrkDeform,true)

    return msg("(✓) Done! exiting...")
end
allocCheck(40,"P","MC")