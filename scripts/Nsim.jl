# Initialisation
using Printf, LinearAlgebra, DelimitedFiles
using Plots, LaTeXStrings
using Base.Threads
typeD = Float64  # Precision (double=Float64 or single=Float32)
include("../src/superInclude.jl")
include("sim.jl")

@views function mainNSim()
    savepath = string(pwd(),"/results/")
    nsim = 5
    for sim in 1:nsim
        mpD = main()
        savefig(string(savepath,"sim_",sim,"epII.png"))
        writedlm(string(savepath,"xp_sim",sim,".txt"),vec(mpD.xp[:,1]))
        writedlm(string(savepath,"zp_sim",sim,".txt"),vec(mpD.xp[:,2]))
        writedlm(string(savepath,"up_sim",sim,".txt"),vec(mpD.up[:,1]))
        writedlm(string(savepath,"wp_sim",sim,".txt"),vec(mpD.up[:,2]))
        writedlm(string(savepath,"epII_sim",sim,".txt"),vec(mpD.epII))
        writedlm(string(savepath,"coh0_sim",sim,".txt"),vec(mpD.coh./1e3))
    end
end
mainNSim()



# https://techytok.com/lesson-parallel-computing/