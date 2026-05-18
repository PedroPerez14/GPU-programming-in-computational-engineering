#include <iostream>
#include <cuda_runtime.h>

// Kernel function to add elements of two arrays
__global__ void vectorAdd(int N, double *A, double *B, double *C) {
    
    int ithread = threadIdx.x;
    int iblock = blockIdx.x; 
    int idx = ( blockIdx.x * blockDim.x ) + threadIdx.x;

    if (idx < N) {     
        C[idx] = A[idx] + B[idx];
    
        printf("<<gpu>> Thread %d in block %d - d_C[%d] = %.0lf \n", ithread, iblock, idx, C[idx]);
    
    }

}


__device__ void addDoubleData(double *A, double *B, double *C) {
     
    (*C) = (*A) + (*B);

}


__global__ void deviceVectorAdd(int N, double *A, double *B, double *C) {
    
    int ithread = threadIdx.x;
    int iblock = blockIdx.x; 
    int idx = ( blockIdx.x * blockDim.x ) + threadIdx.x;

    if (idx < N) { 
        addDoubleData( &(A[idx]), &(B[idx]), &(C[idx]) );    
        
        printf("<<gpu>> Thread %d in block %d - d_C[%d] = %.0lf \n", ithread, iblock, idx, C[idx]);  
    }

}


int main() {

    cudaError_t err; // CUDA error
    double exec_time=0.0; //timer 
    clock_t stime1, stime2;
    int threadsPerBlock;
    int blocksPerGrid;



    // Size of the vectors
    int N = 1000; // 1 million elements

    // Allocate memory on the host
    size_t size = N*sizeof(double);
    double *h_A = (double*) malloc(size);
    double *h_B = (double*) malloc(size);
    double *h_C = (double*) malloc(size);

    // Initialize host vectors
    for (int i = 0; i < N; i++) {
        h_A[i] = (double)i;
        h_B[i] = 0.5;
    }



    // Non-coalescent array
    const int stride = 8; //Separation between memory positions
    double *h_A_nc = (double*) malloc(size*stride);
    double *h_B_nc = (double*) malloc(size*stride);

    // Initialize host vectors
    for (int i = 0; i < N*stride; i++) {
        h_A_nc[i] = (double)(i/stride);
        
        //std::cout << "  h_A_nc[" << i << "] = " << h_A_nc[i] << std::endl;

        h_B_nc[i] = 0.5;
    }    

  



   // Select a GPU
    int selectedGPU = 0; // Change this to select your asigned GPU
    cudaSetDevice(selectedGPU);  

    // Get the device propoerties
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, selectedGPU);    

    // Get total and free memory
    size_t freeMemory, totalMemory;
    cudaMemGetInfo(&freeMemory, &totalMemory);

    std::cout << "  Total Global Memory: " << totalMemory / (1024 * 1024) << " MB" << std::endl;
    std::cout << "  Free Memory: " << freeMemory / (1024 * 1024) << " MB" << std::endl;





    // Allocate memory on the device
    double *d_A, *d_B, *d_C;
    cudaMalloc((void**) &d_A, size);
    cudaMalloc((void**) &d_B, size);
    cudaMalloc((void**) &d_C, size);

    // Non-coalescent arrays allocation on the device
    double *d_A_nc, *d_B_nc;
    cudaMalloc((void**) &d_A_nc, size*stride);
    cudaMalloc((void**) &d_B_nc, size*stride);

    cudaMemGetInfo(&freeMemory, &totalMemory);
    std::cout << "  Updated free Memory: " << freeMemory / (1024 * 1024) << " MB" << std::endl;




    // Copy data from host to device
    cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);

    // Non-coalescent arrays copy from host to device
    cudaMemcpy(d_A_nc, h_A_nc, size*stride, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B_nc, h_B_nc, size*stride, cudaMemcpyHostToDevice);    

    getchar(); //Pause point



    //Launch parallel kernel
    //Start IO time .....................................
    stime1=clock();

    //All threads in a single block
    std::cout << "  Number of tasks: " << N << std::endl;
    std::cout << "  Maximun threads per block: " << prop.maxThreadsPerBlock << std::endl;
    getchar();

    vectorAdd <<<1, N>>> (N, d_A, d_B, d_C);


    //Launch the optimized kernel
    /*
    threadsPerBlock = 256;
    blocksPerGrid = N/threadsPerBlock+1;

    vectorAdd <<<blocksPerGrid, threadsPerBlock>>> (N, d_A, d_B, d_C);
    */

    //Launch device funtion
    /*
    threadsPerBlock = 256;
    blocksPerGrid = N/threadsPerBlock+1;

    deviceVectorAdd <<<blocksPerGrid, threadsPerBlock>>> (N, d_A, d_B, d_C); 
    */


    // Sincronizar y comprobar errores
    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        std::cerr << "Error en cudaDeviceSynchronize: " << cudaGetErrorString(err) << std::endl;
    }

    stime2=clock();
    exec_time += double(stime2-stime1)/CLOCKS_PER_SEC;    
    std::cerr << "  Execution time: " << exec_time << std::endl; 
    //End IO time .....................................  

    getchar(); //Pause point





    // Copy the result from device to host
    cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);

    //Verify the result
    for (int i = 0; i < N; i++) {
        std::cout << "  h_C[" << i << "] = " << h_C[i] << std::endl;
    }




    // Free device memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    // Check the updated free memory
    cudaDeviceSynchronize();
    cudaMemGetInfo(&freeMemory, &totalMemory);
    std::cout << "  Final free Memory: " << freeMemory / (1024 * 1024) << " MB" << std::endl;


    // Free host memory
    free(h_A);
    free(h_B);
    free(h_C);

    return 0;
}
