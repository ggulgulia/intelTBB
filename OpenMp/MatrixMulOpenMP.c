#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h>

void serialMatrixMultiply(double *A, double *B, double *C, size_t rows, size_t columns, size_t k)
{
    for(size_t row = 0; row< rows; ++row)
    {
        for(size_t col = 0; col< columns; ++col)
        {
            C[row*columns + col] = A[row*k + 0]*B[0*columns + col];
            for(size_t i = 1; i < k; ++i)
            {
                C[row*columns + col] += A[row*k + i]*B[i*columns +col];
            }
        }
    }
}

void parallelMatrixMultiply(double *A, double *B, double *C, size_t rows, size_t columns, size_t k)
{
    size_t row, col;
    #pragma omp parallel num_threads(4)
    {
        #pragma omp for schedule(static, 32) private(row, col)
        for(row = 0; row< rows; ++row)
        {
           // #pragma omp for schedule(static, 32) private(col)
            for(col = 0; col< columns; ++col)
            {
                C[row*columns + col] = A[row*k + 0]*B[0*columns + col];
                for(size_t i = 0; i < k; ++i)
                {
                    C[row*columns + col] += A[row*k + i]*B[i*columns +col];
                }
            }
        }
    }
}

int main()
{
    printf("Application to test the parallel speedup from OpenMP\n");
    
    const size_t max_size = 2048;
    const size_t ref_size = 256; //for incrementing the size of matrix on each sample
    int samples = max_size/ref_size ;  //number of samples to be run for serial and tbb 
    
    size_t max_sizeSq = max_size*max_size;
    
    double *A;
    double *B;
    double *C;
    
    A = (double*) calloc(max_sizeSq, sizeof(double));
    B = (double*) calloc(max_sizeSq, sizeof(double));
    C = (double*) calloc(max_sizeSq, sizeof(double));
    
    struct timespec begin1, end1, begin2, end2;
    double diff1, diff2;
    
    int rows_t[samples];
    int columns_t[samples];
    int k_t[samples];
    
    for(int i=0; i<samples; ++i)
    {
        rows_t[i] = (i+1)*ref_size;
        columns_t[i] = (i+1)*ref_size;
        k_t[i] = (i+1)*ref_size;
    }
    
    for(int i = 0;  i<samples; ++i)
    {
        printf("\n===============================================================\n");
        printf("Begin of test %i\n", i+1);
        printf("Reference Size of square matrix: %i\n\n", rows_t[i]);
        
        clock_gettime(CLOCK_MONOTONIC, &begin1);
        serialMatrixMultiply(A, B, C, rows_t[i], columns_t[i], k_t[i]);
        clock_gettime(CLOCK_MONOTONIC, &end1);
        
        diff1 = ((double)end1.tv_sec + 1.0e-9*end1.tv_nsec) -
                ((double)begin1.tv_sec + 1.0e-9*begin1.tv_nsec);  
                
        double totalFlops = (2*k_t[i] -1)*rows_t[i]*columns_t[i]; 
        double gflops_per_sec = totalFlops/(diff1*1000000000);
            
        printf("time for serial multiply: %f\n", diff1);
        printf("total glops needed by serial multiply: %f\n\n", gflops_per_sec);
        
        
        clock_gettime(CLOCK_MONOTONIC, &begin2);
        parallelMatrixMultiply(A, B, C, rows_t[i], columns_t[i], k_t[i]);
        clock_gettime(CLOCK_MONOTONIC, &end2);
        diff2 = ((double)end2.tv_sec + 1.0e-9*end2.tv_nsec) -
                ((double)begin2.tv_sec + 1.0e-9*begin2.tv_nsec);
                
        double gflops_per_sec1 = totalFlops/(diff2*1000000000);
        
        printf("time for parallel multiply: %f\n", diff2); 
        printf("total glops needed by parallel multiply: %f\n\n", gflops_per_sec1);
        
        printf("\nEnd of test %i\n", i+1);
       printf("===============================================================\n\n");
        
    } 

    
    free(A);
    free(B);
    free(C);
    
    return 0;
}
    
