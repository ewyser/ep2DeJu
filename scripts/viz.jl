# Initialisation
using Printf, LinearAlgebra, DelimitedFiles
using Plots, LaTeXStrings
using Statistics
using Base.Threads
typeD = Float64  # Precision (double=Float64 or single=Float32)

@views function main()
    savepath = string(pwd(),"/results/")
    nsim = 20
    for sim in 1:nsim
        x = readdlm(string(savepath,"xp_sim",sim,".txt"),'\n', Float64)
        z = readdlm(string(savepath,"zp_sim",sim,".txt"),'\n', Float64)
        e = readdlm(string(savepath,"epII_sim",sim,".txt"),'\n', Float64)
        if sim == 1
            xp   = x
            zp   = z
            epII = e
        else
            xp   = hcat(xp,x)
            zp   = hcat(zp,z)
            epII = hcat(epII,e)
        end
        sleep(0.1)
        @printf("\rImporting & averaging outputs... [%d/%d]",sim,nsim)                      
    end
    xp=mean(xp,dims=2)
    zp=mean(zp,dims=2)
    epII=mean(epII,dims=2)

    
    gr(size=(2*250,2*125),legend=true,markersize=2.5)
    scatter(xp,zp,zcolor=epII,
            markershape=:circle,
            label="",
            show=true,
            aspect_ratio=1,
            c=:viridis,
            clims=(0.0,2.0),
            ylim=(-10.0,20.0),
            )
    savefig(string(savepath,"epII_average.png"))
end
main()



# https://techytok.com/lesson-parallel-computing/