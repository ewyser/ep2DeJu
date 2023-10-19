# include dependencies & function call(s) for svSolver.jl
using Pkg, LinearAlgebra, Plots, LaTeXStrings, Base.Threads,ProgressMeter
# arithmetic precision (double=Float64 or single=Float32)
const typeD = Float64  
# relative path for figs & data
const path_plot = "./docs/out/"
if isdir(path_plot)==false mkdir(path_plot) end
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
    @info "ϵp2De() init by-default mode"
    include("./fun_fs/shpfun.jl")
    include("./fun_fs/mapsto.jl")
    include("./fun_fs/solve.jl")
    include("./fun_fs/elastoplast.jl")
    include("./fun_fs/plast.jl")
end
