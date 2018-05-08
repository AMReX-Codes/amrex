#include <AMReX.H>

#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

#include "StructOfArrays.H"

using namespace amrex;

struct computeKey
{
    __host__ __device__
    computeKey() {}
    
    template <typename Tuple>
    __host__ __device__
    int operator()(Tuple tup) const
    {
        return thrust::get<3>(tup);
    }
};

int main(int argc, char* argv[]) {
    
    StructOfArrays<2, 2> soa;
    
    for (int i = 100; i >= 0; --i) {
        soa.GetRealData(0).push_back(i);
        soa.GetRealData(1).push_back(i);
        
        soa.GetIntData(0).push_back(i);
        soa.GetIntData(1).push_back(i);
    }
    
    thrust::device_vector<int> keys;
    keys.resize(soa.size());
    
    thrust::transform(thrust::device, soa.begin(), soa.end(), keys.begin(), computeKey());
    thrust::sort_by_key(thrust::device, keys.begin(), keys.end(), soa.begin());
    
    StructOfArrays<2, 2> soa2;
    soa2.GetRealData(0).push_back(101);
    soa2.GetRealData(1).push_back(101);
    
    soa2.GetIntData(0).push_back(101);
    soa2.GetIntData(1).push_back(101);
    
    soa.insert(soa.end(), soa2.begin(), soa2.end());

    for (int i = 0; i < static_cast<int>(soa.size()); ++i) {
        std::cout << soa.GetRealData(0)[i] << " " << soa.GetRealData(0)[i] << " " << soa.GetIntData(0)[i] << " " << soa.GetIntData(1)[i] << std::endl;
    }   
}
