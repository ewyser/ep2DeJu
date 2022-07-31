mutable struct mesh
    nel::Array{Float64}
    nno::Array{Int64}
    nn ::Int64
    L  ::Array{Float64}
    h  ::Array{Float64}
    xn ::Array{Float64}
    zn ::Array{Float64}
    mn ::Array{Float64}
    fen::Array{Float64}
    fin::Array{Float64}
    fn ::Array{Float64}
    an ::Array{Float64}
    pn ::Array{Float64}
    vn ::Array{Float64}
    un ::Array{Float64}
    pel::Array{Float64}
    e2n::Array{Int64}
    xB ::Array{Float64}
end
mutable struct point
    # scalars & vectors
    nmp ::Int64
    l0  ::Array{Float64}
    l   ::Array{Float64}
    v0  ::Array{Float64}
    v   ::Array{Float64}
    mp  ::Array{Float64}
    xp  ::Array{Float64}
    up  ::Array{Float64}
    vp  ::Array{Float64}
    pp  ::Array{Float64}
    coh ::Array{Float64}
    cohr::Array{Float64}
    phi ::Array{Float64}
    epII::Array{Float64}
    J   ::Array{Float64}
    # tensors
    ΔF  ::Array{Float64}
    F   ::Array{Float64}
    b   ::Array{Float64}
    bT  ::Array{Float64}
    ϵ   ::Array{Float64}
    ω   ::Array{Float64}
    σ   ::Array{Float64}
    τ   ::Array{Float64}
    dev ::Array{Float64}
    ep  ::Array{Float64}
    # additional quantities
    ϕ   ::Array{Float64}
    ∂ϕx ::Array{Float64}
    ∂ϕz ::Array{Float64}
    B   ::Array{Float64}
    # connectivity
    p2e ::Array{Int64}
    p2n ::Array{Int64}
end
mutable struct boundary
    x  ::Array{Int64}
    z  ::Array{Int64}
end


struct test
    X  ::Array{Float64}
end