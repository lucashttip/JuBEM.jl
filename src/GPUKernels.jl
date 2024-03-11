using CUDA

function gpu_nonsingint(y)

    index = threadIdx().x
    stride = blockDim().x

    for i = index:stride:length(y)
        
    end

end