# include dependencies & function call(s) for svSolver.jl
using Printf, LinearAlgebra, DelimitedFiles, Plots, LaTeXStrings, Base.Threads,ProgressMeter

include("./fun_fs/types.jl")
include("./fun_fs/topol.jl")
include("./fun_fs/BSpline.jl")
include("./fun_fs/accum.jl")
include("./fun_fs/solve.jl")
include("./fun_fs/flip.jl")
include("./fun_fs/DMBC.jl")
include("./fun_fs/elast.jl")
include("./fun_fs/plast.jl")

include("./fun/functionsT.jl")
include("./fun/RFS.jl")