#ifndef AMREX_ALG_VECTOR_H_
#define AMREX_ALG_VECTOR_H_
#include <AMReX_Config.H>

#include <AMReX_AlgPartition.H>
#include <AMReX_FabArray.H>
#include <AMReX_INT.H>
#include <AMReX_LayoutData.H>
#include <AMReX_TableData.H>

#include <fstream>
#include <string>
#include <type_traits>
#include <utility>

namespace amrex {

template <typename T, typename Allocator = DefaultAllocator<T> >
class AlgVector
{
public:
    using value_type = T;
    using allocator_type = Allocator;

    using Vec = PODVector<T,Allocator>;

    AlgVector () = default;
    explicit AlgVector (Long global_size);
    explicit AlgVector (AlgPartition partition);

    AlgVector (AlgVector<T, Allocator> const&) = delete;
    AlgVector& operator= (AlgVector<T, Allocator> const&) = delete;

    AlgVector (AlgVector<T, Allocator> &&) noexcept = default;
    AlgVector& operator= (AlgVector<T, Allocator> &&) noexcept = default;

    ~AlgVector () = default;

    void define (Long global_size);
    void define (AlgPartition partition);

    [[nodiscard]] bool empty () const { return m_partition.empty(); }

    [[nodiscard]] AlgPartition const& partition () const { return m_partition; }

    [[nodiscard]] Long numLocalRows () const { return m_end - m_begin; }
    [[ nodiscard]] Long numGlobalRows () const { return m_partition.numGlobalRows(); }

    //! Inclusive global index begin.
    [[nodiscard]] Long globalBegin () const { return m_begin; }
    //! Exclusive global index end.
    [[nodiscard]] Long globalEnd () const { return m_end; }

    //! Return raw data pointer. Unlike view(), the starting index is
    //! zero. Use the local row index when accessing the data.
    [[nodiscard]] T const* data () const { return m_data.data(); }
    //! Return raw data pointer. Unlike view(), the starting index is
    //! zero. Use the local row index when accessing the data.
    [[nodiscard]] T      * data ()       { return m_data.data(); }

    //! Return a Table1D view of the data. Unlike data(), indexing starts at
    //! the global row index. Use the global row index when accessing the
    //! data.
    [[nodiscard]] AMREX_FORCE_INLINE
    Table1D<T const, Long> view () const {
        return Table1D<T const, Long>{m_data.data(), m_begin, m_end};
    }

    //! Return a Table1D view of the data. Unlike data(), indexing starts at
    //! the global row index. Use the global row index when accessing the
    //! data.
    [[nodiscard]] AMREX_FORCE_INLINE
    Table1D<T const, Long> const_view () const {
        return Table1D<T const, Long>{m_data.data(), m_begin, m_end};
    }

    //! Return a Table1D view of the data. Unlike data(), indexing starts at
    //! the global row index. Use the global row index when accessing the
    //! data.
    [[nodiscard]] AMREX_FORCE_INLINE
    Table1D<T, Long> view () {
        return Table1D<T, Long>{m_data.data(), m_begin, m_end};
    }

    //! Set every element to `val`. For GPU builds, this function
    //! synchronizes with the host before returning.
    void setVal (T val);

    //! Set every element to `val`. For GPU builds, this function
    //! is asynchronous with respect to the host.
    void setValAsync (T val);

    //! Copy values from `rhs` into this vector. For GPU builds, this
    //! function synchronizes with the hostr before returning.
    void copy (AlgVector<T,Allocator> const& rhs);

    //! Copy values from `rhs` into this vector. For GPU guilds, this
    //! function is asynchronous with respect to the host.
    void copyAsync (AlgVector<T,Allocator> const& rhs);

    //! Add `rhs` to this vector in-place. For GPU builds, this function
    //! synchronizes with the hostr before returning.
    void plus (AlgVector<T,Allocator> const& rhs);

    //! Add `rhs` to this vector in-place. For GPU guilds, this function is
    //! asynchronous with respect to the host.
    void plusAsync (AlgVector<T,Allocator> const& rhs);

    //! Scale this vector by `scale_factor` in-place. For GPU builds, this
    //! function synchronizes with the hostr before returning.
    void scale (T scale_factor);

    //! Scale this vector by `scale_factor` in-place. For GPU guilds, this
    //! function is asynchronous with respect to the host.
    void scaleAsync (T scale_factor);

    /**
     * \brief Return the sum of elements.
     *
     * By default, this returns the global sum over all MPI ranks. If
     * `local` is true, it returns the sum over the locally stored entries
     * only.
     *
     * \param local If true, compute the sum over local entries only.
     * \return The sum of elements.
     */
    [[nodiscard]] T sum (bool local = false) const;

    /**
     * \brief Return the infinity norm.
     *
     * By default, this returns the global norm over all MPI ranks. If
     * `local` is true, it returns the norm over the locally stored entries
     * only.
     *
     * \param local If true, compute the norm over local entries only.
     * \return The infinity norm.
     */
    [[nodiscard]] T norminf (bool local = false) const;

    /**
     * \brief Return the 1-norm.
     *
     * By default, this returns the global norm over all MPI ranks. If
     * `local` is true, it returns the norm over the locally stored entries
     * only.
     *
     * \param local If true, compute the norm over local entries only.
     * \return The 1-norm.
     */
    [[nodiscard]] T norm1 (bool local = false) const;

    /**
     * \brief Return the 2-norm.
     *
     * By default, this returns the global norm over all MPI ranks. If
     * `local` is true, it returns the norm over the locally stored entries
     * only.
     *
     * \param local If true, compute the norm over local entries only.
     * \return The 2-norm.
     */
    [[nodiscard]] T norm2 (bool local = false) const;

    /**
     * \brief Return the 1-norm and 2-norm as a pair.
     *
     * By default, this returns the global norms over all MPI ranks. If
     * `local` is true, it returns the norms over the locally stored entries
     * only.
     *
     * \param local If true, compute the norms over local entries only.
     * \return The 1-norm and 2-norm as a pair.
     */
    [[nodiscard]] std::pair<T,T> norm1and2 (bool local = false) const;

    template <typename FAB,
              std::enable_if_t<amrex::IsBaseFab<FAB>::value &&
                               std::is_same_v<T,typename FAB::value_type>, int> = 0>
    void copyFrom (FabArray<FAB> const& fa);

    template <typename FAB,
              std::enable_if_t<amrex::IsBaseFab<FAB>::value &&
                               std::is_same_v<T,typename FAB::value_type>,int> = 0>
    void copyTo (FabArray<FAB> & fa) const;

    /*
     * \brief Print the vector to files.
     *
     * Each process writes its local portion of the vector to a separate
     * file. This function is provided for debugging purpose. Using it in
     * large-scale runs may overwhelm the file system.
     *
     * File format:
     * - The first line contains two integers describing the local vector:
     *   the global row begin index and the global row end index.
     * - This is followed by one line per element. Each line contains the
     *   global row index and vector value.
     *
     * \param file Base file name. The full name on process `i` is `{file}.{i}`.
     */
    void printToFile (std::string const& file) const;

private:
    AlgPartition m_partition;
    Long m_begin = 0;
    Long m_end = 0;
    Vec m_data;
};

template <typename T, typename Allocator>
AlgVector<T,Allocator>::AlgVector (Long global_size)
    : m_partition(global_size),
      m_begin(m_partition[ParallelDescriptor::MyProc()]),
      m_end(m_partition[ParallelDescriptor::MyProc()+1]),
      m_data(m_end-m_begin)
{}

template <typename T, typename Allocator>
AlgVector<T,Allocator>::AlgVector (AlgPartition partition)
    : m_partition(std::move(partition)),
      m_begin(m_partition[ParallelDescriptor::MyProc()]),
      m_end(m_partition[ParallelDescriptor::MyProc()+1]),
      m_data(m_end-m_begin)
{}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::define (Long global_size)
{
    m_partition.define(global_size);
    m_begin = m_partition[ParallelDescriptor::MyProc()];
    m_end = m_partition[ParallelDescriptor::MyProc()+1];
    m_data.resize(m_end-m_begin);
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::define (AlgPartition partition)
{
    m_partition = std::move(partition);
    m_begin = m_partition[ParallelDescriptor::MyProc()];
    m_end = m_partition[ParallelDescriptor::MyProc()+1];
    m_data.resize(m_end-m_begin);
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::setVal (T val)
{
    setValAsync(val);
    Gpu::streamSynchronize();
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::setValAsync (T val)
{
    Long n = m_data.size();
    T* p = m_data.data();
    ParallelForOMP(n, [=] AMREX_GPU_DEVICE (Long i) noexcept { p[i] = val; });
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::copy (AlgVector<T,Allocator> const& rhs)
{
    copyAsync(rhs);
    Gpu::streamSynchronize();
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::copyAsync (AlgVector<T,Allocator> const& rhs)
{
    Long n = m_data.size();
    AMREX_ASSERT(m_data.size() == rhs.m_data.size());
    T* dst = m_data.data();
    T const* src = rhs.data();
    Gpu::dtod_memcpy_async(dst, src, n*sizeof(T));
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::plus (AlgVector<T,Allocator> const& rhs)
{
    plusAsync(rhs);
    Gpu::streamSynchronize();
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::plusAsync (AlgVector<T,Allocator> const& rhs)
{
    Long n = m_data.size();
    AMREX_ASSERT(m_data.size() == rhs.m_data.size());
    T* dst = m_data.data();
    T const* src = rhs.data();
    ParallelForOMP(n, [=] AMREX_GPU_DEVICE (Long i) noexcept { dst[i] += src[i]; });
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::scale (T scale_factor)
{
    scaleAsync(scale_factor);
    Gpu::streamSynchronize();
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::scaleAsync (T scale_factor)
{
    Long n = m_data.size();
    T* p = m_data.data();
    ParallelForOMP(n, [=] AMREX_GPU_DEVICE (Long i) noexcept { p[i] *= scale_factor; });
}

template <typename T, typename Allocator>
T AlgVector<T,Allocator>::sum (bool local) const
{
    Long n = m_data.size();
    T const* p = m_data.data();
    T r = Reduce::Sum<T>(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
                         {
                             return p[i];
                         });
    if (!local) {
        ParallelAllReduce::Sum(r, ParallelContext::CommunicatorSub());
    }
    return r;
}

template <typename T, typename Allocator>
T AlgVector<T,Allocator>::norminf (bool local) const
{
    Long n = m_data.size();
    T const* p = m_data.data();
    T r = Reduce::Max<T>(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
                         {
                             return amrex::Math::abs(p[i]);
                         });
    if (!local) {
        ParallelAllReduce::Max(r, ParallelContext::CommunicatorSub());
    }
    return r;
}

template <typename T, typename Allocator>
T AlgVector<T,Allocator>::norm1 (bool local) const
{
    Long n = m_data.size();
    T const* p = m_data.data();
    T r = Reduce::Sum<T>(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
                         {
                             return amrex::Math::abs(p[i]);
                         });
    if (!local) {
        ParallelAllReduce::Sum(r, ParallelContext::CommunicatorSub());
    }
    return r;
}

template <typename T, typename Allocator>
T AlgVector<T,Allocator>::norm2 (bool local) const
{
    Long n = m_data.size();
    T const* p = m_data.data();
    T r = Reduce::Sum<T>(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
                         {
                             return p[i]*p[i];
                         });
    if (!local) {
        ParallelAllReduce::Sum(r, ParallelContext::CommunicatorSub());
    }
    return std::sqrt(r);
}

template <typename T, typename Allocator>
std::pair<T,T> AlgVector<T,Allocator>::norm1and2 (bool local) const
{
    Long n = m_data.size();
    T const* p = m_data.data();
    ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
    ReduceData<T,T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (Long i) -> ReduceTuple
                   {
                       return {amrex::Math::abs(p[i]), p[i]*p[i]};
                   });
    auto hv = reduce_data.value(reduce_op);
    if (local) {
        return std::make_pair(amrex::get<0>(hv), std::sqrt(amrex::get<1>(hv)));
    } else {
        std::array<T,2> a{amrex::get<0>(hv), amrex::get<1>(hv)};
        ParallelAllReduce::Sum(a.data(), 2, ParallelContext::CommunicatorSub());
        return std::make_pair(a[0], std::sqrt(a[1]));
    }
}

template <typename T, typename Allocator>
template <typename FAB, std::enable_if_t<amrex::IsBaseFab<FAB>::value &&
                                         std::is_same_v<T,typename FAB::value_type>, int> >
void AlgVector<T,Allocator>::copyFrom (FabArray<FAB> const& fa)
{
    AMREX_ASSERT(fa.is_cell_centered());

    LayoutData<T*> dptrs(fa.boxArray(), fa.DistributionMap());
    T* p = m_data.data();
    for (MFIter mfi(fa); mfi.isValid(); ++mfi) {
        dptrs[mfi] = p;
        p += mfi.validbox().numPts();
    }

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
#pragma omp parallel
#endif
    for (MFIter mfi(fa); mfi.isValid(); ++mfi) {
        fa[mfi].template copyToMem<RunOn::Device>(mfi.validbox(), 0, 1, dptrs[mfi]);
    }
}

template <typename T, typename Allocator>
template <typename FAB, std::enable_if_t<amrex::IsBaseFab<FAB>::value &&
                                         std::is_same_v<T,typename FAB::value_type>, int> >
void AlgVector<T,Allocator>::copyTo (FabArray<FAB> & fa) const
{
    AMREX_ASSERT(fa.is_cell_centered());

    LayoutData<T const*> dptrs(fa.boxArray(), fa.DistributionMap());
    T const* p = m_data.data();
    for (MFIter mfi(fa); mfi.isValid(); ++mfi) {
        dptrs[mfi] = p;
        p += mfi.validbox().numPts();
    }

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
#pragma omp parallel
#endif
    for (MFIter mfi(fa); mfi.isValid(); ++mfi) {
        fa[mfi].template copyFromMem<RunOn::Device>(mfi.validbox(), 0, 1, dptrs[mfi]);
    }
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::printToFile (std::string const& file) const
{
    std::ofstream ofs(file+"."+std::to_string(ParallelDescriptor::MyProc()));
    ofs << m_begin << " " << m_end << "\n";
#ifdef AMREX_USE_GPU
    Gpu::PinnedVector<T> hv(m_data.size());
    Gpu::dtoh_memcpy_async(hv.data(), m_data.data(), m_data.size()*sizeof(T));
    Gpu::streamSynchronize();
    T const* p = hv.data();
#else
    T const* p = m_data.data();
#endif
    for (Long i = 0, N = m_data.size(); i < N; ++i) {
        ofs << i+m_begin << " " << p[i] << "\n";
    }
}

}

#endif
