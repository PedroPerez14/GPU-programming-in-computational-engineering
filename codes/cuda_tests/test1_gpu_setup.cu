#include <iostream>
#include <cuda_runtime.h>

__global__ void hello() {
    printf("\n************");
    printf("\nHi from GPU!");
    printf("\n************\n\n");
}

int main() {

    // Step 1: Get the number of GPUs available
    int numGPUs = 0;
    //******* complete code here *******/
    cudaGetDeviceCount(&numGPUs);
    std::cout << "GPUs number: " << numGPUs << std::endl;
    //**********************************/

    if (numGPUs == 0) {
        std::cerr << "No GPUs found!" << std::endl;
        return -1;
    }



    // Step 2: List properties of each GPU
    for (int i = 0; i < numGPUs; i++) {
        //******* complete code here *******/
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        //**********************************/

        std::cout << "GPU NAME: " << prop.name << std::endl;
    } 


    // Step 3: Select a GPU (e.g., GPU 0)
    //******* complete code here *******/
    int selectedDeviceID = 0;
    cudaSetDevice(selectedDeviceID);
    //**********************************/
    
    

    // Step 4: Verify the selected GPU
    //******* complete code here *******/
    cudaDeviceProp selectedProp;
    cudaGetDeviceProperties(&selectedProp, selectedDeviceID);
    std::cout << "CUDA COMPUTE CAPABILITY: " << selectedProp.major << "." << selectedProp.minor << std::endl;
    std::cout << "GPU SMs: " << selectedProp.multiProcessorCount << std::endl;
    std::cout << "MAX THREADS PER BLOCK: " << selectedProp.maxThreadsPerBlock << std::endl;
    std::cout << "MAX BLOCKS PER GRID: " << selectedProp.maxGridSize[0] << std::endl;
    //**********************************/
    



    // Step 5: Get total and free memory
    size_t freeMemory, totalMemory;
    //******* complete code here *******/
    cudaMemGetInfo(&freeMemory, &totalMemory);
    std::cout << "FREE MEMORY: " << freeMemory << std::endl;
    std::cout << "TOTAL MEMORY: " << freeMemory << std::endl;
    //**********************************/




    // Launch kernel
    hello <<<1,1>>> ();

    // Check for launch errors
    //******* complete code here *******/
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        std::cout << "ERROR: " << cudaGetErrorString(err) << std::endl;
    }
    //**********************************/



    // Synchronize and check for execution errors
    //******* complete code here *******/
    cudaError_t syncErr = cudaDeviceSynchronize();
    if(syncErr != cudaSuccess){
        std::cout << "Sync error: " << cudaGetErrorString(syncErr) << std::endl;
        return -1;
    }
    //**********************************/



    // Flush stdout to ensure printf output is displayed
    fflush(stdout);

    return 0;
}
