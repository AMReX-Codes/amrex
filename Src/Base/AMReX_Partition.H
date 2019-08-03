#ifndef AMREX_PARTITION_H_
#define AMREX_PARTITION_H_

#include <AMReX_Gpu.H>
#include <AMReX_Scan.H>
#include <algorithm>

namespace amrex {

#ifdef AMREX_USE_GPU
// The return value is the total number of trues.
template <typename T, typename F>
int Partition (Gpu::DeviceVector<T>& v, F && f)
{
    int n = v.size();
    Gpu::DeviceVector<T> v2(n);
    T const* AMREX_RESTRICT pv = v.dataPtr();
    T * AMREX_RESTRICT pv2 = v2.dataPtr();
    int tot = Scan::PrefixSum<int> (n,
         [=] AMREX_GPU_DEVICE (int i) -> int
         {
             return f(pv[i]);
         },
         [=] AMREX_GPU_DEVICE (int i, int const& s)
         {
             // We store true elements from the beginning and false
             // elements reversely from the end.  If all elements
             // before pv[i] are true, the exclusive sum so far would
             // be i.  But the actual value is s.
             if (f(pv[i])) {
                 // For true element, s spots from the beginning have
                 // been taken.
                 pv2[s] = pv[i];
             } else {
                 // There are i-s elements before this element that
                 // are false.  From the end, i-s spots have been
                 // taken.
                 pv2[n-1-(i-s)] = pv[i];
             }
         },
         Scan::Type::exclusive);
    v.swap(v2);
    return tot;
}

template <typename T, typename F>
int StablePartition (Gpu::DeviceVector<T>& v, F && f)
{
    int n = Partition(v, std::forward<F>(f));
    int n2 = static_cast<int>(v.size()) - n;
    if (n2 > 1) {
        int npairs = n2/2;
        T* p = v.dataPtr() + n;
        amrex::ParallelFor(npairs, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            std::swap(p[i], p[n2-1-i]);
        });
        Gpu::synchronize();
    }
    return n;
}

#else

// The return value is the total number of trues.
template <typename T, typename F>
int Partition (Gpu::DeviceVector<T>& v, F && f)
{
    auto it = std::partition(v.begin(), v.end(), f);
    return static_cast<int>(std::distance(b.begin(), it));
}

template <typename T, typename F>
int StablePartition (Gpu::DeviceVector<T>& v, F && f)
{
    auto it = std::stable_partition(v.begin(), v.end(), f);
    return static_cast<int>(std::distance(b.begin(), it));
}

#endif

}

#endif
