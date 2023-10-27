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
    xlab,ylab = L"$x-$direction",L"$z-$direction"
    gr(size=(2*250,2*125),legend=true,markersize=2.5,markershape=:circle,markerstrokewidth=0.75,)#markerstrokecolor=:match,)
    temp = L"$t = $"*string(round(t,digits=1))*" [s]"
    if type == "P"
        p = -(mpD.σ[1,:]+mpD.σ[2,:]+mpD.σ[3,:])/3/1e3
        scatter(mpD.x[:,1],mpD.x[:,2],zcolor=p,
            xlabel = xlab,
            ylabel = ylab,
            label=L"$p=-\left(\sigma_{xx,p}+\sigma_{yy,p}+\sigma_{zz,p}\right)/3$",
            aspect_ratio=1,
            c=:viridis,
            ylim=(-10.0,20.0),
            title="pressure: "*temp,
            show=true,
            )  
    elseif type == "epII"
        scatter(mpD.x[:,1],mpD.x[:,2],zcolor=mpD.ϵpII,
            xlabel = xlab,
            ylabel = ylab,    
            label=L"$\epsilon_{\mathrm{II}}^{\mathrm{acc}}$",
            aspect_ratio=1,
            c=:viridis,
            clims=(0.0,2.0),
            ylim=(-10.0,20.0),
            title="plastic strain: "*temp,
            show=true,
            ) 
    elseif type == "du"
        scatter(mpD.x[:,1],mpD.x[:,2],zcolor=sqrt.(mpD.u[:,1].^2+mpD.u[:,2].^2),
            markershape=:circle,
            xlabel = xlab,
            ylabel = ylab,
            label=L"$\Delta u$",
            aspect_ratio=1,
            c=:viridis,
            ylim=(-10.0,20.0),
            title="displacement: "*temp,
            show=true,
            )
    else
        err_msg = "$(type): plot option undefined"
        throw(error(err_msg))
    end
    return ctr+=1
end