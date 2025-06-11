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

    [[nodiscard]] T const* data () const { return m_data.data(); }
    [[nodiscard]] T      * data ()       { return m_data.data(); }

    [[nodiscard]] AMREX_FORCE_INLINE
    Table1D<T const, Long> view () const {
        return Table1D<T const, Long>{m_data.data(), m_begin, m_end};
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    Table1D<T const, Long> const_view () const {
        return Table1D<T const, Long>{m_data.data(), m_begin, m_end};
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    Table1D<T, Long> view () {
        return Table1D<T, Long>{m_data.data(), m_begin, m_end};
    }

    void setVal (T val);
    void setValAsync (T val);

    void copy (AlgVector<T> const& rhs);
    void copyAsync (AlgVector<T> const& rhs);

    void plus (AlgVector<T> const& rhs);
    void plusAsync (AlgVector<T> const& rhs);

    void scale (T scale_factor);
    void scaleAsync (T scale_factor);

    [[nodiscard]] T sum (bool local = false) const;

    [[nodiscard]] T norminf (bool local = false) const;
    [[nodiscard]] T norm2 (bool local = false) const;

    template <typename FAB,
              std::enable_if_t<amrex::IsBaseFab<FAB>::value &&
                               std::is_same_v<T,typename FAB::value_type>, int> = 0>
    void copyFrom (FabArray<FAB> const& fa);

    template <typename FAB,
              std::enable_if_t<amrex::IsBaseFab<FAB>::value &&
                               std::is_same_v<T,typename FAB::value_type>,int> = 0>
    void copyTo (FabArray<FAB> & fa) const;

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
{
    static_assert(std::is_floating_point<T>::value, "AlgVector is for floating point type only");
}

template <typename T, typename Allocator>
AlgVector<T,Allocator>::AlgVector (AlgPartition partition)
    : m_partition(std::move(partition)),
      m_begin(m_partition[ParallelDescriptor::MyProc()]),
      m_end(m_partition[ParallelDescriptor::MyProc()+1]),
      m_data(m_end-m_begin)
{
    static_assert(std::is_floating_point<T>::value, "AlgVector is for floating point type only");
}

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
    ParallelFor(n, [=] AMREX_GPU_DEVICE (Long i) noexcept { p[i] = val; });
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::copy (AlgVector<T> const& rhs)
{
    copyAsync(rhs);
    Gpu::streamSynchronize();
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::copyAsync (AlgVector<T> const& rhs)
{
    Long n = m_data.size();
    AMREX_ASSERT(m_data.size() == rhs.m_data.size());
    T* dst = m_data.data();
    T const* src = rhs.data();
#ifdef AMREX_USE_GPU
    Gpu::dtod_memcpy_async(dst, src, n*sizeof(T));
#else
    std::memcpy(dst, src, n*sizeof(T));
#endif
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::plus (AlgVector<T> const& rhs)
{
    plusAsync(rhs);
    Gpu::streamSynchronize();
}

template <typename T, typename Allocator>
void AlgVector<T,Allocator>::plusAsync (AlgVector<T> const& rhs)
{
    Long n = m_data.size();
    AMREX_ASSERT(m_data.size() == rhs.m_data.size());
    T* dst = m_data.data();
    T const* src = rhs.data();
    ParallelFor(n, [=] AMREX_GPU_DEVICE (Long i) noexcept { dst[i] += src[i]; });
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
    ParallelFor(n, [=] AMREX_GPU_DEVICE (Long i) noexcept { p[i] *= scale_factor; });
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
    T const* p = m_data;
#endif
    for (Long i = 0, N = m_data.size(); i < N; ++i) {
        ofs << i+m_begin << " " << p[i] << "\n";
    }
}

template <class V, class Enable = void> struct IsAlgVector : std::false_type {};
//
template <class V>
struct IsAlgVector<V, std::enable_if_t<std::is_same_v<AlgVector<typename V::value_type,
                                                                typename V::allocator_type>,
                                                      V> > >
    : std::true_type {};

template <typename V1, typename F>
std::enable_if_t<IsAlgVector<std::decay_t<V1> >::value>
ForEach (V1 & x, F const& f)
{
    Long n = x.numLocalRows();
    auto* px = x.data();
    ParallelFor(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        f(px[i]);
    });
}

template <typename V1, typename V2, typename F>
std::enable_if_t<IsAlgVector<std::decay_t<V1> >::value &&
                 IsAlgVector<std::decay_t<V2> >::value>
ForEach (V1 & x, V2 & y, F const& f)
{
    AMREX_ASSERT(x.numLocalRows() == y.numLocalRows());
    Long n = x.numLocalRows();
    auto* AMREX_RESTRICT px = x.data();
    auto* AMREX_RESTRICT py = y.data();
    ParallelFor(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        f(px[i], py[i]);
    });
}

template <typename V1, typename V2, typename V3, typename F>
std::enable_if_t<IsAlgVector<std::decay_t<V1> >::value &&
                 IsAlgVector<std::decay_t<V2> >::value &&
                 IsAlgVector<std::decay_t<V3> >::value>
ForEach (V1 & x, V2 & y, V3 & z, F const& f)
{
    AMREX_ASSERT(x.numLocalRows() == y.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == z.numLocalRows());
    Long n = x.numLocalRows();
    auto* AMREX_RESTRICT px = x.data();
    auto* AMREX_RESTRICT py = y.data();
    auto* AMREX_RESTRICT pz = z.data();
    ParallelFor(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        f(px[i], py[i], pz[i]);
    });
}

template <typename V1, typename V2, typename V3, typename V4, typename F>
std::enable_if_t<IsAlgVector<std::decay_t<V1> >::value &&
                 IsAlgVector<std::decay_t<V2> >::value &&
                 IsAlgVector<std::decay_t<V3> >::value &&
                 IsAlgVector<std::decay_t<V4> >::value>
ForEach (V1 & x, V2 & y, V3 & z, V4 & a, F const& f)
{
    AMREX_ASSERT(x.numLocalRows() == y.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == z.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == a.numLocalRows());
    Long n = x.numLocalRows();
    auto* AMREX_RESTRICT px = x.data();
    auto* AMREX_RESTRICT py = y.data();
    auto* AMREX_RESTRICT pz = z.data();
    auto* AMREX_RESTRICT pa = a.data();
    ParallelFor(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        f(px[i], py[i], pz[i], pa[i]);
    });
}

template <typename V1, typename V2, typename V3, typename V4, typename V5, typename F>
std::enable_if_t<IsAlgVector<std::decay_t<V1> >::value &&
                 IsAlgVector<std::decay_t<V2> >::value &&
                 IsAlgVector<std::decay_t<V3> >::value &&
                 IsAlgVector<std::decay_t<V4> >::value &&
                 IsAlgVector<std::decay_t<V5> >::value>
ForEach (V1 & x, V2 & y, V3 & z, V4 & a, V5 & b, F const& f)
{
    AMREX_ASSERT(x.numLocalRows() == y.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == z.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == a.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == b.numLocalRows());
    Long n = x.numLocalRows();
    auto* AMREX_RESTRICT px = x.data();
    auto* AMREX_RESTRICT py = y.data();
    auto* AMREX_RESTRICT pz = z.data();
    auto* AMREX_RESTRICT pa = a.data();
    auto* AMREX_RESTRICT pb = b.data();
    ParallelFor(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        f(px[i], py[i], pz[i], pa[i], pb[i]);
    });
}

template <typename T>
T Dot (AlgVector<T> const& x, AlgVector<T> const& y, bool local = false)
{
    AMREX_ASSERT(x.numLocalRows() == y.numLocalRows());
    Long n = x.numLocalRows();
    auto const* px = x.data();
    auto const* py = y.data();
    T r = Reduce::Sum<T>(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
                         {
                             return px[i] * py[i];
                         });
    if (!local) {
        ParallelAllReduce::Sum(r, ParallelContext::CommunicatorSub());
    }
    return r;
}

template <typename T>
void Axpy (AlgVector<T>& y, T a, AlgVector<T> const& x, bool async = false)
{
    ForEach(y, x, [=] AMREX_GPU_DEVICE (T& yi, T const& xi) { yi += a*xi; });
    if (!async) { Gpu::streamSynchronize(); }
}

template <typename T>
void LinComb (AlgVector<T>& y, T a, AlgVector<T> const& xa, T b, AlgVector<T> const& xb, bool async = false)
{
    ForEach(y, xa, xb, [=] AMREX_GPU_DEVICE (T& yi, T const& xai, T const& xbi) {
                           yi = a*xai + b*xbi;
                       });
    if (!async) { Gpu::streamSynchronize(); }
}

}

#endif
