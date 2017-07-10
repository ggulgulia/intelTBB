#include <iostream>
#include "tbb/tbb.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range2d.h"
#include <chrono>

using namespace tbb;

const size_t m_rows = 2048;
const size_t m_columns = 2048;
const size_t m_k = 2048;

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

class MatrixMultiply2D
{
    private:
        double *m_A;
        double *m_B;
        double *m_C;

    public:
        void operator()(const blocked_range2d<size_t>& rnge ) const
        {
            double *usrA = m_A;
            double *usrB = m_B;
            double *usrC = m_C;

            for(size_t row=rnge.rows().begin(); row<rnge.rows().end(); ++row)
            {
                for(size_t col=rnge.cols().begin(); col <rnge.cols().end(); ++col)
                {
                    double sum = 0.0;
                
                     for(size_t k=0; k<m_k; ++k)
                     {
                        usrC[row*(rnge.cols().end()-rnge.cols().begin()) + col]  += usrA[row*m_k + k]*usrB[k*(rnge.cols().end()-rnge.cols().begin()) + col];     
                    }
                }
            }
        }

        MatrixMultiply2D(double *A, double *B, double *C):
        m_A(A), m_B(B), m_C(C)
        {       }
};

// this line spawns error, if commented the code compiles fine!
void parallelMatrixMultiply(double *A, double *B, double *C, size_t rows, size_t columns, size_t k)
{
    parallel_for( blocked_range2d<size_t>(0, rows, 32, 0, columns, 32), MatrixMultiply2D(A,B,C) );
}


int main(int argc, char **argv)
{
    double *A = new double[m_rows*m_k];
    double *B = new double[m_k*m_columns];
    double *C = new double[m_rows*m_columns];

    std::chrono::time_point<std::chrono::system_clock> start, start1;
    std::chrono::time_point<std::chrono::system_clock> end, end1;

/***************************************************************************/
/*********************** CALLING SERIAL MATRIX *****************************/
    start = std::chrono::system_clock::now();

    serialMatrixMultiply(A, B, C, m_rows, m_columns, m_k);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedTime = end - start;
    double totalFlops = (2*m_k -1)*m_rows*m_columns;
    
    double flops_per_sec = totalFlops/elapsedTime.count();
    double gflops_per_sec = flops_per_sec/1000000000;
    std::cout << "elaplsed time " << elapsedTime.count() << "\n";
    std:: cout << " total glops needed by serial multiply:" << gflops_per_sec << "\n";
    
/************************ END OF SERIAL MATRIX *********************************/
/*******************************************************************************/


/***************************************************************************/
/********************* CALLING TBB PARALLEL MATRIX *************************/
    start1 = std::chrono::system_clock::now();

    parallelMatrixMultiply(A, B, C, m_rows, m_columns, m_k);

    end1 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedTime1 = end1 - start1;
    double totalFlops1 = (2*m_k -1)*m_rows*m_columns;
    
    double flops_per_sec1 = totalFlops1/elapsedTime1.count();
    double gflops_per_sec1 = flops_per_sec1/1000000000;
    std::cout << "elaplsed time " << elapsedTime1.count() << "\n";
    std:: cout << " total glops needed by parallel multiply:" << gflops_per_sec1 << "\n"; 
    
/********************* END OF PARALLEL MATRIX *********************************/
/******************************************************************************/

    double speedup = elapsedTime.count()/elapsedTime1.count();
    std::cout << "Application speed up : " << speedup << "\n";

    delete[] A;
    delete[] B;
    delete[] C;
   
    return 0;
}
