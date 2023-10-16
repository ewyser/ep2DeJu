# include dependencies & function call(s) for svSolver.jl
using Pkg, LinearAlgebra, Plots, LaTeXStrings, Base.Threads,ProgressMeter
# include doc for: help?> ϵp2De()
include("./misc/doc.jl")
# include init
include("./misc/types.jl")
include("./misc/setup.jl")
include("./misc/utilities.jl")
include("./misc/plot.jl")
# include core functions
if @isdefined perf 
    if perf
        @info "performance mode on: perf = $(perf)"
        include("./misc/rxiv/shpfun.jl")
        include("./misc/rxiv/mapsto.jl")
        include("./fun_fs/solve.jl")
        include("./fun_fs/elastoplast.jl")
        include("./fun_fs/plast.jl")
    else
        @info "performance mode off: perf = $(perf)"
        include("./fun_fs/shpfun.jl")
        include("./fun_fs/mapsto.jl")
        include("./fun_fs/solve.jl")
        include("./fun_fs/elastoplast.jl")
        include("./fun_fs/plast.jl")
    end
else
    @info "performance mode off: perf not existing"
    include("./fun_fs/shpfun.jl")
    include("./fun_fs/mapsto.jl")
    include("./fun_fs/solve.jl")
    include("./fun_fs/elastoplast.jl")
    include("./fun_fs/plast.jl")
end
