#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

__global__ void mandelKernel(float lowerX, float lowerY, float stepX, float stepY, int maxIterations, int resX, int resY, size_t pitch, int* img) {
    // To avoid error caused by the floating number, use the following pseudo code
    //
    // float x = lowerX + thisX * stepX;
    // float y = lowerY + thisY * stepY;

    int thisX = blockIdx.x * blockDim.x + threadIdx.x;
    int thisY = blockIdx.y * blockDim.y + threadIdx.y;

    if (thisX >= resX || thisY >= resY)
        return;

    float x = lowerX + thisX * stepX;
    float y = lowerY + thisY * stepY;

    // int index = thisY * pitch / sizeof(int) + thisX; // Use pitch for proper indexing
    int* index = (int*)((char*)img + thisY * pitch);

    float z_re = x, z_im = y;
    int i;
    for (i = 0; i < maxIterations; ++i) {
        if (z_re * z_re + z_im * z_im > 4.f)
            break;

        float new_re = z_re * z_re - z_im * z_im;
        float new_im = 2.f * z_re * z_im;
        z_re = x + new_re;
        z_im = y + new_im;
    }

    // img[index] = i;
    index[thisX] = i;
}

// Host front-end function that allocates the memory and launches the GPU kernel
void hostFE (float upperX, float upperY, float lowerX, float lowerY, int* img, int resX, int resY, int maxIterations)
{
    float stepX = (upperX - lowerX) / resX;
    float stepY = (upperY - lowerY) / resY;
    size_t pitch = 0;

    int *h_img;
    cudaHostAlloc((void **)&h_img, resX * resY * sizeof(int), cudaHostAllocDefault);

    int *d_img;
    cudaMallocPitch((void **)&d_img, &pitch, resX * sizeof(int), resY);

    dim3 threadsPerBlock(25, 25);
    dim3 numBlocks(resX / threadsPerBlock.x, resY / threadsPerBlock.y);

    mandelKernel<<<numBlocks, threadsPerBlock>>>(lowerX, lowerY, stepX, stepY, maxIterations, resX, resY, pitch, d_img);

    cudaDeviceSynchronize();

    cudaMemcpy2D(h_img, resX * sizeof(int), d_img, pitch, resX * sizeof(int), resY, cudaMemcpyDeviceToHost);
    memcpy(img, h_img, resX * resY * sizeof(int));

    cudaFree(d_img);
    cudaFreeHost(h_img);
}