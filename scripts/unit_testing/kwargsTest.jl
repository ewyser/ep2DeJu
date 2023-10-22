#include("./scripts/unit_testing/kwargsTest.jl")

# include dependencies
include("../../src/superInclude.jl")
# main program
@views function kwargsCheck(nel::Int64,varPlot::String,cmType::String; kwargs...)
    ϕ∂ϕType,fwrkDeform,isΔFbar = getKwargs(kwargs)
    println(ϕ∂ϕType)
    println(fwrkDeform)
    println(isΔFbar)
    return println("[=> done! exiting...")
end

println()
@warn "validation/test"

@info "(Baseline) test 0: by-default call"
kwargsCheck(40,"P","MC")
@info "Test 1: different kwargs"
kwargsCheck(40,"P","MC";shpfun=:gimpm,fwrk=:infinitesimal,vollock=false)
@info "Test 2: missing kwargs"
kwargsCheck(40,"P","MC";shpfun=:gimpm)
@info "Test 3: invalid symbol"
try
    kwargsCheck(40,"P","MC";shpfun=:gifadsfmpm,fwrk=:infinitesimal)
catch err1
    @error err1
end

