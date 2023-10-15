Base.@kwdef struct mesh
    nD   ::Int64
    nel  ::Int64
    nno  ::Int64
    nn   ::Int64
    L    ::Vector{Float64}
    h    ::Vector{Float64}
    x    ::Matrix{Float64}
    # nodal quantities
    m    ::Vector{Float64} 
    fext ::Matrix{Float64}
    fint ::Matrix{Float64}
    D    ::Matrix{Float64}
    f    ::Matrix{Float64}
    a    ::Matrix{Float64}
    p    ::Matrix{Float64}
    v    ::Matrix{Float64}
    u    ::Matrix{Float64}
    pel  ::Matrix{Float64}
    ΔJ   ::Matrix{Float64}
    # mesh-to-node topology
    e2n  ::Matrix{Int64}
    xB   ::Vector{Float64} 
    bc   ::Matrix{Float64}
end
Base.@kwdef struct points
    nmp  ::Int64
    x    ::Matrix{Float64}
    u    ::Matrix{Float64}
    v    ::Matrix{Float64}
    p    ::Matrix{Float64}
    l0   ::Matrix{Float64}
    l    ::Matrix{Float64}
    V0   ::Vector{Float64}
    V    ::Vector{Float64}
    m    ::Vector{Float64}
    coh  ::Vector{Float64}
    cohr ::Vector{Float64}
    phi  ::Vector{Float64}
    ϵpII ::Vector{Float64}
    ϵpV  ::Vector{Float64}
    ΔJ   ::Vector{Float64}
    J    ::Vector{Float64}
    # tensor in matrix notation
    ΔF   ::Array{Float64}
    ΔFbar::Array{Float64}
    F    ::Array{Float64}
    ϵ    ::Array{Float64}
    b    ::Array{Float64}
    bT   ::Array{Float64}
    # tensor in voigt notation
    ω    ::Vector{Float64}
    σR   ::Matrix{Float64}
    σ    ::Matrix{Float64}
    τ    ::Matrix{Float64}
    dev  ::Matrix{Float64}
    ep   ::Matrix{Float64}
    # additional quantities
    ϕ∂ϕ  ::Array{Float64}
    B    ::Array{Float64}
    # connectivity
    p2e  ::Vector{Int64}
    p2n  ::Matrix{Int64}
end

Base.@kwdef struct Point
    x ::Matrix{Float64}
    #test
    a ::Matrix{Float64}
    b ::Matrix{Float64}
    c ::Matrix{Float64}
    d ::Matrix{Float64}
end