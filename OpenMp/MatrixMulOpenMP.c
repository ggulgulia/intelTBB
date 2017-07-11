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
            for(size_t i = 0; i < k; ++i)
            {
                C[row*columns + col] += A[row*k + i]*B[i*columns +col];
            }
        }
    }
}

int main()
{
    printf("Application to test the parallel speedup from OpenMP\n");
    
    const size_t max_size = 1024;
    const size_t ref_size = 256; //for incrementing the size of matrix on each sample
    int samples = max_size/ref_size ;  //number of samples to be run for serial and tbb 
    
    size_t max_sizeSq = max_size*max_size;
    
    double *A;
    double *B;
    double *C;
    
    A = (double*) calloc(max_sizeSq, sizeof(double));
    B = (double*) calloc(max_sizeSq, sizeof(double));
    C = (double*) calloc(max_sizeSq, sizeof(double));
    
    clock_t start1, end1, start2, end2;
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
        printf("Reference Size of square matrix: %i\n", rows_t[i]);
        start1 = clock();
        
        serialMatrixMultiply(A, B, C, rows_t[i], columns_t[i], k_t[i]);
        end1 = clock();
        diff1 = ((double) end1-start1)/CLOCKS_PER_SEC; 
        
        printf("time for serial multiply: %f\n", diff1); 
        
        printf("\nEnd of test %i\n", i+1);
       printf("===============================================================\n\n");
        
    } 

    
    free(A);
    free(B);
    free(C);
    
    return 0;
}
    
