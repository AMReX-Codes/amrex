#ifndef AMREX_PARTITION_H_
#define AMREX_PARTITION_H_

#include <AMReX_Gpu.H>
#include <AMReX_Scan.H>
#include <AMReX_Algorithm.H>

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
                amrex::Swap(p[i], p[n2-1-i]);
            });
            Gpu::synchronize();
        }
    }
}

/**
 * \brief A GPU-capable partition function for contiguous data.
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is not stable, if you want that behavior use 
 * amrex::StablePartition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param data pointer to the data to be partitoned
 * \param beg index at which to start
 * \param end index at which to stop (exclusive)
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */    
template <typename T, typename F>
int Partition (T* data, int beg, int end, F && f)    
{
    int n = end - beg;
    Gpu::DeviceVector<T> v2(n);
    int tot = detail::amrex_partition_helper(data + beg, v2.dataPtr(), n, std::forward<F>(f));
    Gpu::copy(Gpu::deviceToDevice, v2.begin(), v2.end(), data + beg);
    return tot;
}

/**
 * \brief A GPU-capable partition function for contiguous data.
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is not stable, if you want that behavior use 
 * amrex::StablePartition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param data pointer to the data to be partitoned
 * \param N the number of elements in the array
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */        
template <typename T, typename F>
int Partition (T* data, int n, F && f)
{
    return Partition(data, 0, n, std::forward<F>(f));
}

/**
 * \brief A GPU-capable partition function for contiguous data.
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is not stable, if you want that behavior use 
 * amrex::StablePartition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param v a Gpu::DeviceVector with the data to be partitioned.
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */
template <typename T, typename F>
int Partition (Gpu::DeviceVector<T>& v, F && f)
{
    int n = v.size();
    Gpu::DeviceVector<T> v2(n);
    int tot = detail::amrex_partition_helper(v.dataPtr(), v2.dataPtr(), n, std::forward<F>(f));
    v.swap(v2);
    return tot;
}

/**
 * \brief A GPU-capable partition function for contiguous data.
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is stable, meaning that, within each side of the resulting
 * array, order is maintained  - if element i was before element j in the input,
 * then it will also be before j in the output. If you don't care about this property,
 * use amrex::Partition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param data pointer to the data to be partitoned
 * \param beg index at which to start
 * \param end index at which to stop (exclusive)
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */
template <typename T, typename F>
int StablePartition (T* data, int beg, int end, F && f)
{
    int n = Partition(data, beg, end, std::forward<F>(f));
    int n2 = end - beg - n;
    detail::amrex_stable_partition_helper(data + beg + n, n2);
    return n;
}

/**
 * \brief A GPU-capable partition function for contiguous data.
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is stable, meaning that, within each side of the resulting
 * array, order is maintained  - if element i was before element j in the input,
 * then it will also be before j in the output. If you don't care about this property,
 * use amrex::Partition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param data pointer to the data to be partitoned
 * \param N the number of elements in the array
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */            
template <typename T, typename F>
int StablePartition (T* data, int n, F && f)
{
    return StablePartition(data, 0, n, std::forward<F>(f));
}

/**
 * \brief A GPU-capable partition function for contiguous data.
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is stable, meaning that, within each side of the resulting
 * array, order is maintained  - if element i was before element j in the input,
 * then it will also be before j in the output. If you don't care about this property,
 * use amrex::Partition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param v a Gpu::DeviceVector with the data to be partitioned.
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */    
template <typename T, typename F>
int StablePartition (Gpu::DeviceVector<T>& v, F && f)
{
    int n = Partition(v, std::forward<F>(f));
    int n2 = static_cast<int>(v.size()) - n;
    detail::amrex_stable_partition_helper(v.dataPtr() + n, n2);    
    return n;
}

#else

/**
 * \brief A wrapper around std::partition.
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is not stable, if you want that behavior use 
 * amrex::StablePartition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param data pointer to the data to be partitoned
 * \param beg index at which to start
 * \param end index at which to stop (exclusive)
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */
template <typename T, typename F>
int Partition (T* data, int beg, int end, F && f)    
{
    auto it = std::partition(data + beg, data + end, f);
    return static_cast<int>(std::distance(data + beg, it));    
}

/**
 * \brief A wrapper around std::partition.
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is not stable, if you want that behavior use 
 * amrex::StablePartition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param data pointer to the data to be partitoned
 * \param N the number of elements in the array
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */            
template <typename T, typename F>
int Partition (T* data, int n, F && f)
{
    return Partition(data, 0, n, std::forward<F>(f));
}

/**
 * \brief A wrapper around std::partition.
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is not stable, if you want that behavior use 
 * amrex::StablePartition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param v a Gpu::DeviceVector with the data to be partitioned.
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */    
template <typename T, typename F>
int Partition (Gpu::DeviceVector<T>& v, F && f)
{
    auto it = std::partition(v.begin(), v.end(), f);
    return static_cast<int>(std::distance(v.begin(), it));
}

/**
 * \brief A wrapper around std::stable_partition
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is stable, meaning that, within each side of the resulting
 * array, order is maintained  - if element i was before element j in the input,
 * then it will also be before j in the output. If you don't care about this property,
 * use amrex::Partition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param data pointer to the data to be partitoned
 * \param beg index at which to start
 * \param end index at which to stop (exclusive)
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */    
template <typename T, typename F>
int StablePartition (T* data, int beg, int end, F && f)    
{
    auto it = std::stable_partition(data + beg, data + end, f);
    return static_cast<int>(std::distance(data + beg, it));    
}

/**
 * \brief A wrapper around std::stable_partition
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is stable, meaning that, within each side of the resulting
 * array, order is maintained  - if element i was before element j in the input,
 * then it will also be before j in the output. If you don't care about this property,
 * use amrex::Partition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param data pointer to the data to be partitoned
 * \param N the number of elements in the array
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */                
template <typename T, typename F>
int StablePartition (T* data, int n, F && f)
{
    return StablePartition(data, 0, n, std::forward<F>(f));
}

/**
 * \brief A wrapper around std::stable_partition
 *
 * After calling this, all the items for which the predicate is true
 * will be before the items for which the predicate is false in the
 * input array. 
 *
 * This version is stable, meaning that, within each side of the resulting
 * array, order is maintained  - if element i was before element j in the input,
 * then it will also be before j in the output. If you don't care about this property,
 * use amrex::Partition instead.
 * 
 * \tparam T type of the data to be partitioned.
 * \tparam F type of the predicate function.
 *
 * \param v a Gpu::DeviceVector with the data to be partitioned.
 * \param f predicate function that returns 1 or 0 for each input
 *
 * Returns the index of the first element for which f is 0.
 */        
template <typename T, typename F>
int StablePartition (Gpu::DeviceVector<T>& v, F && f)
{
    auto it = std::stable_partition(v.begin(), v.end(), f);
    return static_cast<int>(std::distance(v.begin(), it));
}

#endif

}

#endif
