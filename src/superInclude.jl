# activate project
using Pkg
if splitpath(Base.active_project())[end-1]!="ep2DeJu"
    Pkg.activate(".")
end
# include dependencies & function call(s)
try 
    using LinearAlgebra, Plots, LaTeXStrings, Random, Base.Threads,ProgressMeter
catch
    @error "$(Base.active_project()) needs instantiation"
    @warn "By-default instantiation launched...may take a while"
    Pkg.instantiate()
    using LinearAlgebra, Plots, LaTeXStrings, Base.Threads,ProgressMeter
end
# arithmetic precision & relative path for figs & data
const typeD     = Float64  
const path_plot = "./docs/out/"
if isdir(path_plot)==false mkdir(path_plot) end
const path_test = "./docs/test/"
if isdir(path_test)==false mkdir(path_test) end
# include doc for: help?> ϵp2De()
include("./misc/doc.jl")
# include init
include("./misc/types.jl")
include("./misc/utilities.jl")
include("./misc/setup.jl")
include("./misc/physics.jl")
include("./misc/plot.jl")
# include core functions
if @isdefined perf 
    if perf
        @info "performance mode on: perf = $(perf)"
        include("./fun_fs/shpfun.jl")
        include("./fun_fs/mapsto.jl")
        include("./fun_fs/solve.jl")
        include("./misc/rxiv/elastoplast.jl")
    else
        @info "performance mode off: perf = $(perf)"
        include("./fun_fs/shpfun.jl")
        include("./fun_fs/mapsto.jl")
        include("./fun_fs/solve.jl")
        include("./fun_fs/elastoplast.jl")
        include("./fun_fs/plast.jl")
    end
else
    @info "ϵp23De() init by-default mode"
    include("./fun_fs/shpfun.jl")
    include("./fun_fs/mapsto.jl")
    include("./fun_fs/solve.jl")
    include("./fun_fs/elastoplast.jl")
        include("./fun_fs/RetMap/J2RetMap.jl")
        include("./fun_fs/RetMap/MCRetMap.jl")
        include("./fun_fs/RetMap/DPRetMap.jl")
end
