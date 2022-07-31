# Initialisation
using Printf, LinearAlgebra, DelimitedFiles
using Plots
using Base.Threads
using LaTeXStrings
typeD = Float64  # Precision (double=Float64 or single=Float32)
include("../src/BSpline.jl")

@views function main()

default(titlefont=14, guidefont=14, tickfont=10, legendfont=14)

    @printf("\no---------------------------------------------o");
    @printf("\n|              ** ep23De v1.0 **              |");
    @printf("\no---------------------------------------------o");
    @printf("\n cubic B-Spline check");
    @printf("\no---------------------------------------------o\n"); 



    #

    xn = LinRange(-2, 12, 10)
    xn = xn[:]
    xp = LinRange(xn[3], xn[end-2], 4*80)
    xp = xp[:]
    dx = abs(xn[1]-xn[2])
    println(dx)

    xB = vec(hcat(xn[3],xn[end-2]))
    a = zeros(Float64,length(xp),length(xn),5)
    for mp in 1:length(xp)
        for nn in 1:length(xn)
            # compute basis functions
            ξ      = (xp[mp] - xn[nn])/dx 
            type   = whichType(xn[nn],xB,dx)
                a[mp,nn,end] = 1.0*type 
            ϕx,dϕx = ϕ∇ϕ(ξ,type,dx)
            η      = (xp[mp] - xn[nn])/dx
            type   = whichType(xn[nn],xB,dx)
            ϕz,dϕz = ϕ∇ϕ(η,type,dx)
            # convolution of basis function
            ϕ   =  ϕx*  ϕz                                        
            ∂ϕx = dϕx*  ϕz                                        
            ∂ϕz =  ϕx* dϕz
            a[mp,nn,1] = ϕx
            a[mp,nn,2] = ∂ϕx
            a[mp,nn,3] = ∂ϕz
            a[mp,nn,4] = xp[mp] 
        end
    end
    p = findall(x->x>0.0, vec(a[:,:,1]))
    x = vec(a[:,:,4])
    x = x[p]
    y = vec(a[:,:,1])
    y = y[p]
    c = vec(a[:,:,end])
    c = c[p]

    scatter(x,y,zcolor=c,markershape=:circle,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0)
    scatter!(xn,zeros(size(xn)),color="red",markersize=5,xlabel=L"x",ylabel=L"\phi_n(x_p)",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-2,2),colorbar_title=L"type",levels=5)
    savefig("Nsplineplot.png")

    x = vec(a[:,:,4])
    x = x[p]
    y = vec(a[:,:,2])
    y = y[p]
    c = vec(a[:,:,end])
    c = c[p]
    scatter(x,y,zcolor=c,markershape=:circle,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0)
    scatter!(xn,zeros(size(xn)),color="red",markersize=5,xlabel=L"x",ylabel=L"\partial_x\phi_n(x_p)",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-2,2),colorbar_title=L"type",levels=5)
    savefig("dNsplineplot.png")

    #=
    mn  = reshape(meD.mn,meD.nno[3],meD.nno[1],meD.nno[2])
    xn  = reshape(meD.xn,meD.nno[3],meD.nno[1],meD.nno[2])
    zn  = reshape(meD.zn,meD.nno[3],meD.nno[1],meD.nno[2])
    idy = convert(Int64,ceil(meD.nno[2]/2))
    mn  = mn[:,:,idy]
    xn  = xn[:,:,idy]
    zn  = zn[:,:,idy]
    #display(heatmap(1:meD.nno[1], 1:meD.nno[3], mn, xlabel="lx", ylabel="ly", title="lumped mass matrix",aspect_ratio=1))
    p = -(mpD.σ[1,:]+mpD.σ[2,:]+mpD.σ[3,:])/3.0/1e3
    scatter!(mpD.xp[:,1],mpD.xp[:,3],marker_z=p,markershape=:circle,label="",show=true,aspect_ratio=1)
    savefig("plot.png")
    =#
end
main()



# https://techytok.com/lesson-parallel-computing/