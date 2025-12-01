#include <stdio.h>
#include <stdlib.h>
#include "hostFE.h"
#include "helper.h"

void hostFE(int filterWidth, float *filter, int imageHeight, int imageWidth,
            float *inputImage, float *outputImage, cl_device_id *device,
            cl_context *context, cl_program *program)
{
    cl_int status;
    int filterSize = filterWidth * filterWidth;
    
    cl_command_queue command_queue = clCreateCommandQueue(*context, *device, 0, &status);
    
    cl_mem inputImageBuffer = clCreateBuffer(*context, CL_MEM_READ_ONLY, imageHeight * imageWidth * sizeof(float), NULL, &status);
    cl_mem outputImageBuffer = clCreateBuffer(*context, CL_MEM_WRITE_ONLY, imageHeight * imageWidth * sizeof(float), NULL, &status);
    cl_mem filterBuffer = clCreateBuffer(*context, CL_MEM_READ_ONLY, filterSize * sizeof(float), NULL, &status);
    
    status = clEnqueueWriteBuffer(command_queue, inputImageBuffer, CL_TRUE, 0, imageHeight * imageWidth * sizeof(float), inputImage, 0, NULL, NULL);
    status = clEnqueueWriteBuffer(command_queue, filterBuffer, CL_TRUE, 0, filterSize * sizeof(float), filter, 0, NULL, NULL);

    cl_kernel kernel = clCreateKernel(*program, "convolution", &status);
    status = clSetKernelArg(kernel, 0, sizeof(cl_mem), &inputImageBuffer);
    status = clSetKernelArg(kernel, 1, sizeof(cl_mem), &outputImageBuffer);
    status = clSetKernelArg(kernel, 2, sizeof(cl_mem), &filterBuffer);
    status = clSetKernelArg(kernel, 3, sizeof(int), &filterWidth);
    status = clSetKernelArg(kernel, 4, sizeof(int), &imageWidth);
    status = clSetKernelArg(kernel, 5, sizeof(int), &imageHeight);

    size_t global_work_size[2] = {imageWidth, imageHeight};
    
    status = clEnqueueNDRangeKernel(command_queue, kernel, 2, NULL, global_work_size, NULL, 0, NULL, NULL);
    
    status = clFinish(command_queue);
    
    status = clEnqueueReadBuffer(command_queue, outputImageBuffer, CL_TRUE, 0, imageHeight * imageWidth * sizeof(float), outputImage, 0, NULL, NULL);
    
    clReleaseMemObject(inputImageBuffer);
    clReleaseMemObject(outputImageBuffer);
    clReleaseMemObject(filterBuffer);
    clReleaseCommandQueue(command_queue);
    clReleaseKernel(kernel);
}
