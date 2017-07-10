#include "tbb/tbb.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include <iostream>
#include <chrono>

using namespace tbb;

const size_t size = 10000000000;

double serialReduce(double  *array, size_t n)
{
    double sum = 0.0;
    for(size_t i = 0; i<n; i++)
    {
        sum += array[i];
    }
    
    return sum;
}

class myReduce
{
    private:
        
        double *m_array;
        
    public:
       double m_sum;
       void operator() (const blocked_range<size_t> &range)
       {
            double *array = m_array;
            double sum = m_sum;
            size_t end = range.end();
            for(size_t i=range.begin(); i<end; ++i)
            {
                sum += array[i];
            }

            m_sum = sum;
        }
        
        myReduce(myReduce &object, split) : m_array(object.m_array), m_sum(0.0) 
        {   }
        
        myReduce( double *array_t): m_array(array_t), m_sum(0.0)
        {   }
        
        void join(const myReduce &y) 
        {
            m_sum += y.m_sum;
        }
        
        ~myReduce()
        {
            delete[] m_array;
        }
};

double parallelReduce(double *array, size_t n)
{
    myReduce reduce(array);
    parallel_reduce(blocked_range<size_t> (0, 25000000), reduce);
    
    return reduce.m_sum;
}

int main(int argc, char **argv)
{
    
    std::cout << "\n Application to compare serial reduce and parallel reduce using TBB\n\n";
    double *dummyArray = new double[size];
    
    for(size_t i = 0; i<size; ++i)
    {
        dummyArray[i]  = 1.0;
    }
    
    double serialResult, parallelResult;
    
    
    

/*********************** CALLING SERIAL REDUCE *****************************/
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    
   // serialResult = serialReduce(dummyArray, size);
    
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedTime = end - start;
    std::cout << "elapsed time for serial reduce = " << elapsedTime.count() << "\n";

/************************ END OF SERIAL REDUCE *********************************/


/*********************** CALLING PARALLEL REDUCE *****************************/
    std::chrono::time_point<std::chrono::system_clock> start1, end1;
    start1 = std::chrono::system_clock::now();
    
    parallelResult = parallelReduce(dummyArray, size);
    
    end1 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedTime1 = end1 - start1;
    std::cout << "elapsed time for parallel reduce = " << elapsedTime1.count() << "\n";

/************************ END OF PARALLEL REDUCE *********************************/
    
    
    return 0;
}
