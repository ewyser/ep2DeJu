# include("./scripts/unit_testing/shpfunTest.jl")
# Initialisation
using LinearAlgebra, Plots, Test, Pkg, LaTeXStrings, Base.Threads
const path_plot = "./docs/out/"
include("../../src/fun_fs/shpfun.jl")
include("../../src/misc/plot.jl")
include("../../src/misc/utilities.jl")

@views function ϕ∂ϕCheck(ϕ∂ϕType)
    xn = LinRange(-2, 12, 10)
    xn = xn[:]
    xp = LinRange(xn[3], xn[end-2], 4*80)
    xp = xp[:]
    dx = abs(xn[1]-xn[2])
    xB = vec(hcat(xn[3],xn[end-2]))
    a  = zeros(Float64,length(xp),length(xn),5)
    PoU = zeros(Float64,length(xp))
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
                PoU[mp]+= ϕx
                a[mp,nn,1] = ϕx
                a[mp,nn,2] = dϕx
                a[mp,nn,3] = dϕz
                a[mp,nn,4] = xp[mp] 
            end
        end
    elseif ϕ∂ϕType == "smpm"
        for mp ∈ eachindex(xp)
            for nn ∈ eachindex(xn)
                # compute basis functions
                ξ      = xp[mp] - xn[nn]
                η      = (xp[mp] - xn[nn])
                ϕx,dϕx = NdN(ξ,dx,0.0)
                ϕz,dϕz = NdN(η,dx,0.0)
                # convolution of basis function
                PoU[mp]+= ϕx
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
    S = vec(a[:,:,1])
    S = S[p]
    c = vec(a[:,:,end])
    c = c[p]

    dS = vec(a[:,:,2])
    dS = dS[p]


    CM    = zeros(RGB{Float64}, 4)
    CM[1] = RGB{Float64}(1,0,0)  # black
    CM[2] = RGB{Float64}(0,1,0) # yellow
    CM[3] = RGB{Float64}(0,0,1) # yellow
    CM[4] = RGB{Float64}(0,0.5,0) # red

    Title = ["Shape functions","Derivatives","Partition of unity (PoU)"]
    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    if ϕ∂ϕType == "bsmpm"
        p1=scatter(x,S,zcolor=c,markershape=:circle,label="",show=true,aspect_ratio=1,cmap=cgrad(CM,4;categorical=true),markerstrokecolor=:auto,markerstrokewidth=0)
        p1=scatter!(xn,zeros(size(xn)),color="red",markersize=5,ylabel=L"\phi_n(x_p)",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-0.1,1.25),colorbar_title="type",levels=5,title=Title[1])
        p2=scatter(x,dS,zcolor=c,markershape=:circle,label="",show=true,aspect_ratio=1,cmap=cgrad(CM,4;categorical=true),markerstrokecolor=:auto,markerstrokewidth=0)
        p2=scatter!(xn,zeros(size(xn)),color="red",markersize=5,ylabel=L"\partial_x\phi_n(x_p)",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-1,1),colorbar_title="type",levels=5,title=Title[2])
        p3=scatter(xp,PoU,zcolor=c,markershape=:circle,label="",show=true,aspect_ratio=1,cmap=cgrad(CM,4;categorical=true),markerstrokecolor=:auto,markerstrokewidth=0)
        p3=scatter!(xn,zeros(size(xn)),color="red",markersize=5,xlabel=L"$x$ [m]",ylabel=L"$\sum_n\phi_n(x_p)$",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-0.1,1.5),colorbar_title="type",levels=5,title=Title[3])
    else
        p1=scatter(x,S,markershape=:circle,color="black",label="",show=true,aspect_ratio=1,markerstrokecolor=:auto,markerstrokewidth=0)
        p1=scatter!(xn,zeros(size(xn)),color="red",markersize=5,ylabel=L"\phi_n(x_p)",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-0.1,1.25),colorbar_title="type",levels=5,title=Title[1])
        p2=scatter(x,dS,markershape=:circle,color="black",label="",show=true,aspect_ratio=1,markerstrokecolor=:auto,markerstrokewidth=0)
        p2=scatter!(xn,zeros(size(xn)),color="red",markersize=5,ylabel=L"\partial_x\phi_n(x_p)",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-1,1),colorbar_title="type",levels=5,title=Title[2])
        p3=scatter(xp,PoU,markershape=:circle,color="black",label="",show=true,aspect_ratio=1,markerstrokecolor=:auto,markerstrokewidth=0)
        p3=scatter!(xn,zeros(size(xn)),color="red",markersize=5,xlabel=L"$x$ [m]",ylabel=L"$\sum_n\phi_n(x_p)$",markershape=:square,label="",show=true,aspect_ratio=1,c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-0.1,1.5),colorbar_title="type",levels=5,title=Title[3])
    end
    display(plot(p1,p2,p3; layout=(3,1), size=(550,500)))
    savefig(path_plot*"check_$(ϕ∂ϕType).png")

    return minimum(PoU),sum(PoU)/length(xp),maximum(PoU)
end

@info "** ϵp2De v$(getVersion()): partition of unity (PoU) testset **"
for shp in ["bsmpm","gimpm","smpm"]
    @testset "$(shp): PoU" begin
        Min,Mean,Max = ϕ∂ϕCheck(shp)
        @test Min  ≈ 1.0 atol=0.001
        @test Mean ≈ 1.0 atol=0.001
        @test Max  ≈ 1.0 atol=0.001
    end
end