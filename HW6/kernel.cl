__kernel void convolution(
    __global float* inputImage, 
    __global float* outputImage, 
    __global float* filter, 
    const int filterWidth, 
    const int imageWidth, 
    const int imageHeight) 
{
    int x = get_global_id(0);
    int y = get_global_id(1);

    int filterHalf = filterWidth / 2;
    float sum = 0.0;
    for (int i = -filterHalf; i <= filterHalf; i++) {
        for (int j = -filterHalf; j <= filterHalf; j++) {
            // Get the input image pixel
            int currentX = x + j;
            int currentY = y + i;

            if (currentX >= 0 && currentX < imageWidth && currentY >= 0 && currentY < imageHeight) {
                // Get the corresponding filter weight
                float filterValue = filter[(i + filterHalf) * filterWidth + (j + filterHalf)];
                // Accumulate the weighted input pixel
                sum += inputImage[currentY * imageWidth + currentX] * filterValue;
            }
        }
    }

    outputImage[y * imageWidth + x] = sum;
}
