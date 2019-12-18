#ifndef AMREX_PARTITION_H_
#define AMREX_PARTITION_H_

#include <AMReX_Gpu.H>
#include <AMReX_Scan.H>

#include <algorithm>

namespace amrex {

#ifdef AMREX_USE_GPU

namespace detail
{
    template <typename T, typename F>    
    int amrex_partition_helper (T const* AMREX_RESTRICT pv, T* AMREX_RESTRICT pv2, int n, F && f)
    {
        return Scan::PrefixSum<int> (n,
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
    }

    template <typename T>
    void amrex_stable_partition_helper (T* p, int n2)
    {
        if (n2 > 1) {
            int npairs = n2/2;
            amrex::ParallelFor(npairs, [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                std::swap(p[i], p[n2-1-i]);
            });
            Gpu::synchronize();
        }
    }
}

template <typename T, typename F>
int Partition (T* data, int beg, int end, F && f)    
{
    int n = end - beg;
    Gpu::DeviceVector<T> v2(n);
    int tot = detail::amrex_partition_helper(data + beg, v2.dataPtr(), n, std::forward<F>(f));
    Gpu::copy(Gpu::deviceToDevice, v2.begin(), v2.end(), data + beg);
    return tot;
}
    
template <typename T, typename F>
int Partition (T* data, int n, F && f)
{
    return Partition(data, 0, n, std::forward<F>(f));
}

template <typename T, typename F>
int Partition (Gpu::DeviceVector<T>& v, F && f)
{
    int n = v.size();
    Gpu::DeviceVector<T> v2(n);
    int tot = detail::amrex_partition_helper(v.dataPtr(), v2.dataPtr(), n, std::forward<F>(f));
    v.swap(v2);
    return tot;
}

template <typename T, typename F>
int StablePartition (T* data, int beg, int end, F && f)
{
    int n = Partition(data, beg, end, std::forward<F>(f));
    int n2 = end - beg;
    detail::amrex_stable_partition_helper(data + beg + n, n2);
    return n;
}
    
template <typename T, typename F>
int StablePartition (T* data, int n, F && f)
{
    return StablePartition(data, 0, n, std::forward<F>(f));
}

template <typename T, typename F>
int StablePartition (Gpu::DeviceVector<T>& v, F && f)
{
    int n = Partition(v, std::forward<F>(f));
    int n2 = static_cast<int>(v.size()) - n;
    detail::amrex_stable_partition_helper(v.dataPtr() + n, n2);    
    return n;
}

#else

template <typename T, typename F>
int Partition (T* data, int beg, int end, F && f)    
{
    auto it = std::partition(data + beg, data + end, f);
    return static_cast<int>(std::distance(data + beg, it));    
}

template <typename T, typename F>
int Partition (T* data, int n, F && f)
{
    return Partition(data, 0, n, std::forward<F>(f));
}
        
template <typename T, typename F>
int Partition (Gpu::DeviceVector<T>& v, F && f)
{
    auto it = std::partition(v.begin(), v.end(), f);
    return static_cast<int>(std::distance(b.begin(), it));
}

template <typename T, typename F>
int StablePartition (T* data, int beg, int end, F && f)    
{
    auto it = std::stable_partition(data + beg, data + end, f);
    return static_cast<int>(std::distance(data + beg, it));    
}

template <typename T, typename F>
int StablePartition (T* data, int n, F && f)
{
    return StablePartition(data, 0, n, std::forward<F>(f));
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
