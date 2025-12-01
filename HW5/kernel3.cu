#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

__global__ void mandelKernel(float lowerX, float lowerY, float stepX, float stepY, int maxIterations, int resX, int resY, size_t pitch, int* img, int pixels) {
    // To avoid error caused by the floating number, use the following pseudo code
    //
    // float x = lowerX + thisX * stepX;
    // float y = lowerY + thisY * stepY;

    int thisX = ( blockIdx.x * blockDim.x + threadIdx.x ) * pixels;
    int thisY = blockIdx.y * blockDim.y + threadIdx.y;

    if (thisY >= resY) return;
    int* index = (int*)((char*)img + thisY * pitch);

    for (int i = 0; i < pixels && (thisX + i) < resX; i++)
    {
        float x = lowerX + (thisX + i) * stepX;
        float y = lowerY + thisY * stepY;

        float z_re = x, z_im = y;
        int j;
        for (j = 0; j < maxIterations; ++j) {
            if (z_re * z_re + z_im * z_im > 4.f) break;

            float new_re = (z_re * z_re) - (z_im * z_im);
            float new_im = 2.f * z_re * z_im;
            z_re = x + new_re;
            z_im = y + new_im;
        }

        
        index[thisX + i] = j;
    }
}

// Host front-end function that allocates the memory and launches the GPU kernel
void hostFE(float upperX, float upperY, float lowerX, float lowerY, int* img, int resX, int resY, int maxIterations) {
    float stepX = (upperX - lowerX) / resX;
    float stepY = (upperY - lowerY) / resY;
    int pixels = 4;

    int* h_img;
    cudaHostAlloc((void**)&h_img, resX * resY * sizeof(int), cudaHostAllocDefault);

    int* d_img;
    size_t pitch = 0;
    cudaMallocPitch((void**)&d_img, &pitch, resX * sizeof(int), resY);

    dim3 blockSize(16, 16);
    dim3 gridSize(resX / ( blockSize.x * pixels ), resY / blockSize.y);

    mandelKernel<<<gridSize, blockSize>>>(lowerX, lowerY, stepX, stepY, maxIterations, resX, resY, pitch, d_img, pixels);

    // cudaMemcpy2D(h_img, resX * sizeof(int), d_img, pitch, resX * sizeof(int), resY, cudaMemcpyDeviceToHost);
    // memcpy(img, h_img, resX * resY + sizeof(int));
    cudaMemcpy2D(h_img, resX * sizeof(int), d_img, pitch, resX * sizeof(int), resY, cudaMemcpyDeviceToHost);
    memcpy(img, h_img, resX * resY * sizeof(int));
    
    cudaFree(d_img);
    cudaFreeHost(h_img);
}