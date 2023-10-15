# julia -i -O3 -t auto --check-bounds=no --project=.
# include("./scripts/test_struct_vs_namedTuple.jl")
# test(40,"P","MC")
# test(40,"P","MC";shpfun="bsmpm",fwrk="finite",vollock=true)

# include dependencies
include("../src/misc/types.jl")
using BenchmarkTools
# arithmetic precision (double=Float64 or single=Float32)
typeD = Float64  
# relative path for figs & data
path_plot = "./out/"
if isdir(path_plot)==false mkdir(path_plot) end

@views function add_namedTuple(mpD_namedTuple)
    for k in 1:1000
        for r in 1:size(mpD_namedTuple.x,1)
            mpD_namedTuple.x[r] = rand()
        end
    end
end
@views function add_struct(mpD_struct)
    for k in 1:1000
        for r in 1:size(mpD_struct.x,1)
            mpD_struct.x[r] = rand()
        end
    end
end

@views function test()
    @warn "comparison struct & NamedTuple"

    mpD_namedTuple = (
        x    = zeros(Float64,10000,10000),
    )
    mpD_struct = Point(
        x = zeros(Float64,10000,10000),
        a = zeros(Float64,1000,1000),
        b = zeros(Float64,1000,1000),
        c = zeros(Float64,1000,1000),
        d = zeros(Float64,1000,1000),
    )

    @time add_namedTuple(mpD_namedTuple)
    @time add_struct(mpD_struct)

    return mpD_struct,mpD_namedTuple
end
test()