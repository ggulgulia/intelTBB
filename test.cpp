 #include <iostream>
#include "tbb/tbb.h"
#include "tbb/parallel_for.h"

using namespace tbb;

class MySeries
{
	public:
		double mySum;
		MySeries(MySeries &x, split): mySum(0) {	}
		MySeries() : mySum(0) {	}
		void operator()(const blocked_range<size_t>& r)
		{
			std::cout << "Calling blocked range operator\n";
			double sum = mySum;
			size_t end = r.end();
			for(long i = r.begin(); i!=end; ++i)
			{
				double denom = i*2.0+1.0;
				if(i%2==1)
				{
					denom = -denom;
				}
				sum += 1.0/denom;
			}
			mySum = sum;
		}
		
		void join(MySeries& y)
		{
			mySum += y.mySum;
		}
};

int main(int argc, char **argv)
{
	MySeries x;
	parallel_reduce(blocked_range<size_t>(0,1000000),x);
	std::cout << 4*x.mySum << std::endl;
	return 0;
}
			
