#include("./scripts/unit_testing/shpfunTest.jl")
# Initialisation
using Printf, LinearAlgebra, DelimitedFiles
using Plots
using Base.Threads
using LaTeXStrings
include("../../src/superInclude.jl")

@views function ϕ∂ϕCheck(ϕ∂ϕType)
    @info "** ϵp2-3De v1.0: $(ϕ∂ϕType) check **"
    @info "init. arbitrary nodal & mp's coordinates"
    xn = LinRange(-2, 12, 10)
    xn = xn[:]
    xp = LinRange(xn[3], xn[end-2], 4*80)
    xp = xp[:]
    dx = abs(xn[1]-xn[2])
    xB = vec(hcat(xn[3],xn[end-2]))
    a  = zeros(Float64,length(xp),length(xn),5)
    PoU = zeros(Float64,length(xp))
    @info "shape functions calculation(s)..."
    if ϕ∂ϕType == "bsmpm"
        for mp ∈ eachindex(xp)
            for nn ∈ eachindex(xn)
                # compute basis functions
                ξ      = (xp[mp] - xn[nn])/dx 
                type   = whichType(xn[nn],xB,dx)
                    a[mp,nn,end] = 1.0*type 
                ϕx,dϕx = ϕ∇ϕ(ξ,type,dx)
                if type != 0
                    PoU[mp]+= ϕx
                end


                η      = (xp[mp] - xn[nn])/dx
                type   = whichType(xn[nn],xB,dx)
                ϕz,dϕz = ϕ∇ϕ(η,type,dx)

                # convolution of basis function
                a[mp,nn,1] = ϕx
                a[mp,nn,2] = dϕx
                a[mp,nn,3] = dϕz
                a[mp,nn,4] = xp[mp] 
            end
        end
    elseif ϕ∂ϕType == "gimpm"
        for mp ∈ eachindex(xp)
            for nn ∈ eachindex(xn)
                # compute basis functions
                ξ      = xp[mp] - xn[nn]
                η      = (xp[mp] - xn[nn])
                ϕx,dϕx = NdN(ξ,dx,0.5)
                ϕz,dϕz = NdN(η,dx,0.5)
                # convolution of basis function
                a[mp,nn,1] = ϕx
                a[mp,nn,2] = dϕx
                a[mp,nn,3] = dϕz
                a[mp,nn,4] = xp[mp] 
            end
        end
    end
    tol = 1e-3
    p = findall(x->x>tol, vec(a[:,:,1]))
    x = vec(a[:,:,4])
    x = x[p]
    y = vec(a[:,:,1])
    y = y[p]
    c = vec(a[:,:,end])
    c = c[p]


    CM    = zeros(RGB{Float64}, 4)
    CM[1] = RGB{Float64}(1,0,0)  # black
    CM[2] = RGB{Float64}(0,1,0) # yellow
    CM[3] = RGB{Float64}(0,0,1) # yellow
    CM[4] = RGB{Float64}(0,0.5,0) # red

    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    scatter(x,y,zcolor=c,markershape=:circle,label="",show=true,aspect_ratio=1,cmap=cgrad(CM,4;categorical=true),markerstrokecolor=:auto,markerstrokewidth=0)
    scatter!(xn,zeros(size(xn)),color="red",markersize=5,xlabel=L"$x$ [m]",ylabel=L"\phi_n(x_p)",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-2,2),colorbar_title="type",levels=5,title=L"$\log_{10}(\phi_n(x_p)) > $"*string(log10(tol)))
    sleep(2.5)
    savefig(path_plot*"check_$(ϕ∂ϕType)_ϕ.png")
    

    y = vec(a[:,:,2])
    y = y[p]

    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    scatter(x,y,zcolor=c,markershape=:circle,label="",show=true,aspect_ratio=1,cmap=cgrad(CM,4;categorical=true),markerstrokecolor=:auto,markerstrokewidth=0)
    scatter!(xn,zeros(size(xn)),color="red",markersize=5,xlabel=L"$x$ [m]",ylabel=L"\partial_x\phi_n(x_p)",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-2,2),colorbar_title="type",levels=5,title=L"$\log_{10}(\phi_n(x_p)) > $"*string(log10(tol)))
    sleep(2.5)
    savefig(path_plot*"check_$(ϕ∂ϕType)_∂ϕ.png")
    
    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    scatter(xp,PoU,zcolor=c,markershape=:circle,label="",show=true,aspect_ratio=1,cmap=cgrad(CM,4;categorical=true),markerstrokecolor=:auto,markerstrokewidth=0)
    scatter!(xn,zeros(size(xn)),color="red",markersize=5,xlabel=L"$x$ [m]",ylabel=L"$\sum_n\phi_n(x_p)$",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-2,2),colorbar_title="type",levels=5,title="partition of unity")
    sleep(2.5)
    savefig(path_plot*"check_$(ϕ∂ϕType)_PoU.png")

    @info "Figs saved in" path_plot
    return println("[=> done! exiting...")
end
ϕ∂ϕCheck("bsmpm")



# https://techytok.com/lesson-parallel-computing/