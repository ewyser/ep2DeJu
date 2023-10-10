#include("./scripts/splineCheck.jl")
# Initialisation
using Printf, LinearAlgebra, DelimitedFiles
using Plots
using Base.Threads
using LaTeXStrings
typeD = Float64  # Precision (double=Float64 or single=Float32)
include("../src/fun_fs/shpfun.jl")

path_plot = "./out/"
if isdir(path_plot)==false mkdir(path_plot) end

default(
    fontfamily="Computer Modern",
    linewidth=2,
    framestyle=:box,
    label=nothing,
    grid=false
    )

@views function ϕ∂ϕCheck(ϕ∂ϕType)
    @info "** ϵp2-3De v1.0: "*ϕ∂ϕType*" check **"
    default(titlefont=14, 
            guidefont=14,  
            tickfont=10, 
            legendfont=14
            )
    @info "init. arbitrary nodal & mp's coordinates"
    xn = LinRange(-2, 12, 10)
    xn = xn[:]
    xp = LinRange(xn[3], xn[end-2], 4*80)
    xp = xp[:]
    dx = abs(xn[1]-xn[2])
    xB = vec(hcat(xn[3],xn[end-2]))
    a  = zeros(Float64,length(xp),length(xn),5)
    @info "shape functions calculation(s)..."
    if ϕ∂ϕType == "bsmpm"
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
                a[mp,nn,1] = ϕx*  ϕz                                        
                a[mp,nn,2] = dϕx*  ϕz                                        
                a[mp,nn,3] =  ϕx* dϕz
                a[mp,nn,4] = xp[mp] 
            end
        end
    elseif ϕ∂ϕType == "gimpm"
        for mp in 1:length(xp)
            for nn in 1:length(xn)
                # compute basis functions
                ξ      = xp[mp] - xn[nn]
                η      = (xp[mp] - xn[nn])
                ϕx,dϕx = NdN(ξ,dx,0.5)
                ϕz,dϕz = NdN(η,dx,0.5)
                # convolution of basis function
                a[mp,nn,1] = ϕx*  ϕz                                        
                a[mp,nn,2] = dϕx*  ϕz                                        
                a[mp,nn,3] =  ϕx* dϕz
                a[mp,nn,4] = xp[mp] 
            end
        end
    end
    p = findall(x->x>0.0, vec(a[:,:,1]))
    x = vec(a[:,:,4])
    x = x[p]
    y = vec(a[:,:,1])
    y = y[p]
    c = vec(a[:,:,end])
    c = c[p]

    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    scatter(x,y,zcolor=c,markershape=:circle,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0)
    scatter!(xn,zeros(size(xn)),color="red",markersize=5,xlabel=L"$x$",ylabel=L"\phi_n(x_p)",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-2,2),colorbar_title="type",levels=5)
    savefig(path_plot*"check_"*ϕ∂ϕType*"_"*"Nsplineplot.png")

    x = vec(a[:,:,4])
    x = x[p]
    y = vec(a[:,:,2])
    y = y[p]
    c = vec(a[:,:,end])
    c = c[p]

    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    scatter(x,y,zcolor=c,markershape=:circle,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0)
    scatter!(xn,zeros(size(xn)),color="red",markersize=5,xlabel=L"$x$",ylabel=L"\partial_x\phi_n(x_p)",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-2,2),colorbar_title="type",levels=5)
    savefig(path_plot*"check_"*ϕ∂ϕType*"_"*"dNsplineplot.png")
    
    @info "Figs saved in" path_plot
    return println("[=> done! exiting...")
end
ϕ∂ϕCheck("bsmpm")



# https://techytok.com/lesson-parallel-computing/