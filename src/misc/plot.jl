# plot parameters
default(
    fontfamily  = "Computer Modern",
    titlefont   = 12, 
    guidefont   = 12,  
    tickfont    = 10, 
    legendfont  = 10,
    linewidth   = 2,
    framestyle  = :box,
    label       = nothing,
    grid        = false
    )
# plot routines
@views function plot_coh(xp,coh,phi,ϕ0)
    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    scatter(xp[:,1],xp[:,2],zcolor=coh./1e3,
    markershape=:circle,
    label="",
    show=true,
    aspect_ratio=1,
    c=:vik,
    clims=(10.0,30.0),
    markerstrokecolor=:auto,
    markerstrokewidth=0,
    ylim=(-10,20),
    )
    savefig(path_plot*"coh0.png")
    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    scatter(xp[:,1],xp[:,2],zcolor=phi,
    markershape=:circle,
    label="",
    show=true,
    aspect_ratio=1,
    c=:vik,
    clims=(ϕ0-ϕ0/5,ϕ0+ϕ0/5),
    markerstrokecolor=:auto,
    markerstrokewidth=0,
    ylim=(-10,20),
    )
    savefig(path_plot*"phi0.png")
end
@views function plotStuff(mpD,t,type,ctr)
    temp = L"$t = $"*string(round(t,digits=1))*" [s]"
    if type == "P"
        d   = -(mpD.σ[1,:]+mpD.σ[2,:]+mpD.σ[3,:])/3/1e3
        lab = L"$p=-\left(\sigma_{xx,p}+\sigma_{yy,p}+\sigma_{zz,p}\right)/3$"
        tit = "pressure, "*temp
    elseif type == "epII"
        d = mpD.ϵpII
        lab = L"$\epsilon_{\mathrm{II}}^{\mathrm{acc}}$"
        tit = "plastic strain, "*temp
    elseif type == "du"
        d = sqrt.(mpD.u[:,1].^2+mpD.u[:,2].^2)
        lab = L"$\Delta u$"
        tit = "displacement, "*temp
    else
        err_msg = "$(type): plot option undefined"
        throw(error(err_msg))
    end
    # plot
    gr(legend=true,markersize=2.5,markershape=:circle,markerstrokewidth=0.75,)#markerstrokecolor=:match,)
    p1 = scatter(mpD.x[:,1],mpD.x[:,2],zcolor=d,
    xlabel = L"$x-$direction",
    ylabel = L"$z-$direction",
    label  = lab,
    aspect_ratio=1,
    c=:viridis,
    ylim=(-10.0,20.0),
    title=tit,
    )
    display(plot(p1;layout=(1,1),size=(500,250))) 
    return ctr+=1
end