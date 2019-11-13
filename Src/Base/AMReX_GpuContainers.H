#ifndef AMREX_GPU_CONTAINERS_H_
#define AMREX_GPU_CONTAINERS_H_

#include <numeric>

#include <AMReX_Vector.H>
#include <AMReX_PODVector.H>
#include <AMReX_GpuAllocators.H>

namespace amrex {

namespace Gpu {

#ifdef AMREX_USE_GPU

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

    static struct HostToDevice {} hostToDevice;
    static struct DeviceToHost {} deviceToHost;
    static struct DeviceToDevice {} deviceToDevice;
    
    template<class InIter, class OutIter>
    void copy (HostToDevice, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;
        auto size = std::distance(begin, end);
        htod_memcpy(&(*result), &(*begin), size*sizeof(value_type));
    }

    template<class InIter, class OutIter>
    void copy (DeviceToHost, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;
        auto size = std::distance(begin, end);
        dtoh_memcpy(&(*result), &(*begin), size*sizeof(value_type));
    }

    template<class InIter, class OutIter>
    void copy (DeviceToDevice, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;
        auto size = std::distance(begin, end);
        dtod_memcpy(&(*result), &(*begin), size*sizeof(value_type));
    }

    template<class InIter, class OutIter>
    void copyAsync (HostToDevice, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;
        auto size = std::distance(begin, end);
        htod_memcpy_async(&(*result), &(*begin), size*sizeof(value_type));
    }

    template<class InIter, class OutIter>
    void copyAsync (DeviceToHost, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;
        auto size = std::distance(begin, end);
        dtoh_memcpy_async(&(*result), &(*begin), size*sizeof(value_type));
    }

    template<class InIter, class OutIter>
    void copyAsync (DeviceToDevice, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;
        auto size = std::distance(begin, end);
        dtod_memcpy_async(&(*result), &(*begin), size*sizeof(value_type));
    }
    
    template<class InIter, class OutIter>
    OutIter inclusive_scan (InIter begin, InIter end, OutIter result)
    {
#ifdef AMREX_USE_GPU
        return thrust::inclusive_scan(Gpu::The_ThrustCachedPolicy(),
                                      begin, end, result);
#else
        return std::partial_sum(begin, end, result);
#endif
    }

    template<class InIter, class OutIter>
    OutIter exclusive_scan(InIter begin, InIter end, OutIter result)
    {
#ifdef AMREX_USE_GPU
        return thrust::exclusive_scan(Gpu::The_ThrustCachedPolicy(),
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
