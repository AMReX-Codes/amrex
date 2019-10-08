#ifndef AMREX_CUDA_CONTAINERS_H_
#define AMREX_CUDA_CONTAINERS_H_

#include <numeric>

#include <AMReX_Vector.H>
#include <AMReX_PODVector.H>
#include <AMReX_CudaAllocators.H>

namespace amrex {

namespace Cuda {

#ifdef AMREX_USE_CUDA

    template <class T>
    using DeviceVector = PODVector<T, ArenaAllocator<T> >;
	
    template <class T>
    using HostVector = PODVector<T>;
    
    template <class T>
    using ManagedDeviceVector = PODVector<T, ManagedArenaAllocator<T> >;
    
    template <class T>
    using ManagedVector = PODVector<T, ManagedArenaAllocator<T> >;

    template <class T>
    using PinnedDeviceVector = PODVector<T, PinnedArenaAllocator<T> >;

    template <class T>
    using PolymorphicDeviceVector = PODVector<T, PolymorphicAllocator<T> >;
    
#else
    //! When Cuda is off, all these containers revert to amrex::Vector. 
    template <class T>
    using DeviceVector = PODVector<T>;

    template <class T>
    using HostVector = PODVector<T>;

    template <class T>
    using ManagedVector = PODVector<T>;

    template <class T>
    using ManagedDeviceVector = PODVector<T>;

    template <class T>
    using PinnedDeviceVector = PODVector<T>;

    template <class T>
    using PolymorphicDeviceVector = PODVector<T>;
#endif

    template<class InIter, class OutIter>
    OutIter thrust_copy (InIter begin, InIter end, OutIter result)
    {
#ifdef AMREX_USE_CUDA
        return thrust::copy(begin, end, result);
#else
        return std::copy(begin, end, result);
#endif
    }

    template<class InIter, class OutIter>
    OutIter inclusive_scan (InIter begin, InIter end, OutIter result)
    {
#ifdef AMREX_USE_CUDA
        return thrust::inclusive_scan(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
                                      begin, end, result);
#else
        return std::partial_sum(begin, end, result);
#endif
    }

    template<class InIter, class OutIter>
    OutIter exclusive_scan(InIter begin, InIter end, OutIter result)
    {
#ifdef AMREX_USE_CUDA
        return thrust::exclusive_scan(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
                                      begin, end, result);
#else
        if (begin == end) return result;
        
        typename std::iterator_traits<InIter>::value_type sum = *begin;
        *result++ = sum - *begin;
        
        while (++begin != end) {
            sum = std::move(sum) + *begin;
            *result++ = sum - *begin;
        }
        return ++result;
#endif
    }
}


}


#endif
