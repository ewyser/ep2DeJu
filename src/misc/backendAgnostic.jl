# https://juliagpu.github.io/KernelAbstractions.jl/stable/
#include("./src/misc/backendAgnostic.jl")
using KernelAbstractions


@kernel function mul2_kernel!(A,nx,ny)
    i,j = @index(Global, NTuple)
    if 1<i<size(A,1) && 1<j<size(A,2)
        A[i,j] = 2 * A[i,j]
    end
end

# allocate array on device memory 
A     = ones(5, 5)#A     = CuArray(A)
nx,ny = size(A)

# get backend where memory is allocated
dev       = get_backend(A)
blockSize = KernelAbstractions.isgpu(get_backend(a)) ? 256 : 1024

# agnostic backend kernel call
mul2_kernel!(dev, blockSize)(A,nx,ny, ndrange=size(A))
synchronize(dev)

# print result
println(A)