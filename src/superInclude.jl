# include dependencies & function call(s) for svSolver.jl
using LinearAlgebra, Plots, LaTeXStrings, Base.Threads,ProgressMeter
# include doc for: help?> ϵp2De()
include("./misc/doc.jl")
# include init
include("./misc/functionsT.jl")
include("./misc/RFG.jl")
# include core functions
include("./fun_fs/shpfun.jl")
include("./fun_fs/mapsto.jl")
include("./fun_fs/solve.jl")
include("./fun_fs/elastoplast.jl")
include("./fun_fs/plast.jl")

