# include dependencies
include("../../src/misc/rxiv/camCGolchinRetMap.jl")
# independant physical constant
g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
E,ν     = 1.0e6,0.3                                                         # Young's mod. [Pa], Poisson's ratio [-]
K,G,Del = D(E,ν,3)                                                     # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
ρ0      = 2700.0                                                            # density [kg/m^3]
yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
cmParam = (E = E, ν = ν, Kc = K, Gc = G, Del = Del,)
# camC param
χ     = 3.0/2.0
pc0   = -cmParam.Kc/3.0
pc,pt = pc0,-0.0*pc0
ϕcs   = 20.0*π/180.0
M     = 6.0*sin(ϕcs)/(3.0-sin(ϕcs))
ζ,γ   = 1.0,1.0
α,β   = 0.0,0.0

p1 = camCplotYieldFun(pc0,pt,γ,M,α,β)
display(plot(p1;layout=(1,1),size=(500,250)))
sleep(2.5)
savefig(path_plot*"pqSpace_camCYieldFun.png")

# include("./scripts/unit_testing/camCRetMapTest.jl")