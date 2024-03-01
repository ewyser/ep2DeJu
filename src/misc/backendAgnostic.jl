# https://juliagpu.github.io/KernelAbstractions.jl/stable/
# include("./src/misc/backendAgnostic.jl")
using KernelAbstractions, CUDA, LinearAlgebra


@kernel inbounds=true function kernel!(A,nx,ny)
    i,j = @index(Global, NTuple)
    if 1<i<size(A,1) && 1<j<size(A,2)
        A[i,j] = 2.0 * A[i,j]
    end
end
@kernel inbounds=true function test!(sig,np)
    k = @index(Global)
    if k<=np
        λ,n = eigen(sig[:,:,k],sortby=nothing)
    end
end

function configDevice(dev)
    blockSize     = KernelAbstractions.isgpu(dev) ? 256 : 1024
    global mul_k! = kernel!(dev, blockSize)
    return @info "kernel launch parameters & alias(es) done"
end

# initialize
A     = ones(1024,1024)#A     = CuArray(A)
nx,ny = size(A)

# allocate memory, get backend & kernel launch parameter
AT = CUDA.functional() ? CUDA.CuArray : Array
A_D= AT(A)

dev = get_backend(A_D)
configDevice(dev)


# agnostic backend kernel call
for k in 1:10000
    mul_k!(A_D,nx,ny;ndrange=size(A_D))
    KernelAbstractions.synchronize(dev)
end

# print result
#println(A)



sig = ones(3,3,10)
sig_D = sig
global eigN! = test!(CPU(), KernelAbstractions.isgpu(dev) ? 256 : 1024)
np = 10
eigN!(sig_D,np, ndrange=np)
KernelAbstractions.synchronize(dev)