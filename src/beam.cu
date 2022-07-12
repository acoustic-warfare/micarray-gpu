#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

extern "C" void cudaBeamWrapper(int *res, const int *first, const int *last, int n_bytes);

__global__ void beamKernel(int *res, const int *a, const int *b, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size)
    {
        res[i] = a[i] * b[i];
    }
}

// Cuda Wrapper for `beamKernel` used by C or Cython code
void cudaBeamWrapper(int *res, const int *first, const int *last, int n_bytes)
{
    // Setup buffers for GPU
    int *dev_res = nullptr;
    int *dev_first = nullptr;
    int *dev_last = nullptr;

    // Allocate memory on GPU for three vectors
    cudaMalloc((void **)&dev_res, n_bytes * sizeof(int));
    cudaMalloc((void **)&dev_first, n_bytes * sizeof(int));
    cudaMalloc((void **)&dev_last, n_bytes * sizeof(int));

    // Copy allocated host memory to device
    cudaMemcpy(dev_first, first, n_bytes * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_last, last, n_bytes * sizeof(int), cudaMemcpyHostToDevice);

    // Compute the result using one thread per element in vector
    // 2 is number of computational blocks and (n_bytes + 1) / 2 is a number of threads in a block
    beamKernel<<<2, (n_bytes + 1) / 2>>>(dev_res, dev_first, dev_last, n_bytes);

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaDeviceSynchronize();

    // Copy output vector from GPU buffer to host memory.
    cudaMemcpy(res, dev_res, n_bytes * sizeof(int), cudaMemcpyDeviceToHost);

    // Release allocated memory
    cudaFree(dev_res);
    cudaFree(dev_first);
    cudaFree(dev_last);

    cudaDeviceReset();
}



// int main(int argc, char **argv)
// {
//     const int arraySize = BYTES;
//     int res[arraySize] = {0};
//     int first[arraySize];
//     int last[arraySize];

//     // Inititate random values
//     int i;
//     for (i = 0; i < BYTES; i++)
//     {
//         first[i] = rand();
//     }
//     for (i = 0; i < BYTES; i++)
//     {
//         last[i] = rand();
//     }
    
//     cudaBeamWrapper(res, first, last, arraySize);
//     int loop;
//     for (loop = 0; loop < BYTES; loop++)
//         printf("%d ", res[loop]);
//     printf("\n");
//     cudaDeviceReset();

//     return 0;
// }