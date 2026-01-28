#ifndef AMREX_REDUCE_H_
#define AMREX_REDUCE_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_OpenMP.H>
#include <AMReX_MFIter.H>
#include <AMReX_ReduceOpType.H>
#include <AMReX_ValLocPair.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <functional>
#include <limits>

namespace amrex {

/// \cond DOXYGEN_IGNORE
namespace Reduce::detail {

#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
    template <std::size_t I, typename T, typename P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void for_each_parallel (T& d, T const& s, Gpu::Handler const& h)
    {
        P().parallel_update(amrex::get<I>(d), amrex::get<I>(s), h);
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void for_each_parallel (T& d, T const& s, Gpu::Handler const& h)
    {
        P().parallel_update(amrex::get<I>(d), amrex::get<I>(s), h);
        for_each_parallel<I+1,T,P1,Ps...>(d, s, h);
    }
#else
    template <std::size_t I, typename T, typename P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void for_each_parallel (T& d, T const& s)
    {
        P().parallel_update(amrex::get<I>(d), amrex::get<I>(s));
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void for_each_parallel (T& d, T const& s)
    {
        P().parallel_update(amrex::get<I>(d), amrex::get<I>(s));
        for_each_parallel<I+1,T,P1,Ps...>(d, s);
    }
#endif
#endif

    template <std::size_t I, typename T, typename P>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void for_each_local (T& d, T const& s)
    {
        P().local_update(amrex::get<I>(d), amrex::get<I>(s));
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void for_each_local (T& d, T const& s)
    {
        P().local_update(amrex::get<I>(d), amrex::get<I>(s));
        for_each_local<I+1,T,P1,Ps...>(d, s);
    }

    template <std::size_t I, typename T, typename P>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr void for_each_init (T& t)
    {
        P().init(amrex::get<I>(t));
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    constexpr void for_each_init (T& t)
    {
        P().init(amrex::get<I>(t));
        for_each_init<I+1,T,P1,Ps...>(t);
    }
}
/// \endcond

/// \ingroup amrex_collectives
struct ReduceOpSum
{

#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s, Gpu::Handler const& h) const noexcept {
        T r = Gpu::blockReduceSum(s,h);
        if (h.threadIdx() == 0) { d += r; }
    }
#else
    template <typename T, int MT=AMREX_GPU_MAX_THREADS>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s) const noexcept {
        T r = Gpu::blockReduceSum<MT>(s);
        if (threadIdx.x == 0) { d += r; }
    }
#endif
#endif

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void local_update (T& d, T const& s) const noexcept { d += s; }

    template <typename T>
    constexpr void init (T& t) const noexcept { t = 0; }
};

/// \ingroup amrex_collectives
struct ReduceOpMin
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s, Gpu::Handler const& h) const noexcept {
        T r = Gpu::blockReduceMin(s,h);
        if (h.threadIdx() == 0) { d = amrex::min(d,r); }
    }
#else
    template <typename T, int MT=AMREX_GPU_MAX_THREADS>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s) const noexcept {
        T r = Gpu::blockReduceMin<MT>(s);
        if (threadIdx.x == 0) { d = amrex::min(d,r); }
    }
#endif
#endif

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void local_update (T& d, T const& s) const noexcept { d = amrex::min(d,s); }

    template <typename T>
    constexpr std::enable_if_t<std::numeric_limits<T>::is_specialized>
    init (T& t) const noexcept { t = std::numeric_limits<T>::max(); }

    template <typename T>
    constexpr std::enable_if_t<!std::numeric_limits<T>::is_specialized>
    init (T& t) const noexcept { t = T::max(); }
};

/// \ingroup amrex_collectives
struct ReduceOpMax
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s, Gpu::Handler const& h) const noexcept {
        T r = Gpu::blockReduceMax(s,h);
        if (h.threadIdx() == 0) { d = amrex::max(d,r); }
    }
#else
    template <typename T, int MT=AMREX_GPU_MAX_THREADS>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s) const noexcept {
        T r = Gpu::blockReduceMax<MT>(s);
        if (threadIdx.x == 0) { d = amrex::max(d,r); }
    }
#endif
#endif

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void local_update (T& d, T const& s) const noexcept { d = amrex::max(d,s); }

    template <typename T>
    constexpr std::enable_if_t<std::numeric_limits<T>::is_specialized>
    init (T& t) const noexcept { t = std::numeric_limits<T>::lowest(); }

    template <typename T>
    constexpr std::enable_if_t<!std::numeric_limits<T>::is_specialized>
    init (T& t) const noexcept { t = T::lowest(); }
};

/// \ingroup amrex_collectives
struct ReduceOpLogicalAnd
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral_v<T>>
    parallel_update (T& d, T s, Gpu::Handler const& h) const noexcept {
        T r = Gpu::blockReduceLogicalAnd(s,h);
        if (h.threadIdx() == 0) { d = d && r; }
    }
#else
    template <typename T, int MT=AMREX_GPU_MAX_THREADS>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral_v<T>>
    parallel_update (T& d, T s) const noexcept {
        T r = Gpu::blockReduceLogicalAnd<MT>(s);
        if (threadIdx.x == 0) { d = d && r; }
    }
#endif
#endif

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral_v<T>>
    local_update (T& d, T s) const noexcept { d = d && s; }

    template <typename T>
    constexpr std::enable_if_t<std::is_integral_v<T>>
    init (T& t) const noexcept { t = true; }
};

/// \ingroup amrex_collectives
struct ReduceOpLogicalOr
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral_v<T>>
    parallel_update (T& d, T s, Gpu::Handler const& h) const noexcept {
        T r = Gpu::blockReduceLogicalOr(s,h);
        if (h.threadIdx() == 0) { d = d || r; }
    }
#else
    template <typename T, int MT=AMREX_GPU_MAX_THREADS>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral_v<T>>
    parallel_update (T& d, T s) const noexcept {
        T r = Gpu::blockReduceLogicalOr<MT>(s);
        if (threadIdx.x == 0) { d = d || r; }
    }
#endif
#endif

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral_v<T>>
    local_update (T& d, T s) const noexcept { d = d || s; }

    template <typename T>
    constexpr std::enable_if_t<std::is_integral_v<T>>
    init (T& t) const noexcept { t = false; }
};

template <typename... Ps> class ReduceOps;

#ifdef AMREX_USE_GPU

/// \ingroup amrex_collectives
template <typename... Ts>
class ReduceData
{
public:
    using Type = GpuTuple<Ts...>;

    template <typename... Ps>
    explicit ReduceData (ReduceOps<Ps...>& reduce_op)
        : m_max_blocks(Gpu::Device::maxBlocksPerLaunch()),
          m_host_tuple((Type*)(The_Pinned_Arena()->alloc(sizeof(Type)))),
          m_device_tuple((Type*)(The_Arena()->alloc((AMREX_GPU_MAX_STREAMS)
                                                    * m_max_blocks * sizeof(Type)))),
          m_fn_value([&reduce_op,this] () -> Type { return this->value(reduce_op); })
    {
        reduce_op.resetResultReadiness();
        static_assert(std::is_trivially_copyable<Type>(),
                      "ReduceData::Type must be trivially copyable");
        static_assert(std::is_trivially_destructible<Type>(),
                      "ReduceData::Type must be trivially destructible");

        new (m_host_tuple) Type();
        m_nblocks.fill(0);
    }

    ~ReduceData () {
        The_Pinned_Arena()->free(m_host_tuple);
        The_Arena()->free(m_device_tuple);
    }

    ReduceData (ReduceData<Ts...> const&) = delete;
    ReduceData (ReduceData<Ts...> &&) = delete;
    void operator= (ReduceData<Ts...> const&) = delete;
    void operator= (ReduceData<Ts...> &&) = delete;

    Type value ()
    {
        return m_fn_value();
    }

    template <typename... Ps>
    Type value (ReduceOps<Ps...> & reduce_op)
    {
        return reduce_op.value(*this);
    }

    Type* devicePtr () { return m_device_tuple; }
    Type* devicePtr (gpuStream_t const& s) {
        return m_device_tuple+(Gpu::Device::streamIndex(s))*m_max_blocks;
    }

    Type* hostPtr () { return m_host_tuple; }

    GpuArray<int,AMREX_GPU_MAX_STREAMS>& nBlocks () { return m_nblocks; }
    int& nBlocks (gpuStream_t const& s) { return m_nblocks[Gpu::Device::streamIndex(s)]; }

    int maxBlocks () const { return m_max_blocks; }

    int maxStreamIndex () const { return m_max_stream_index; }
    void updateMaxStreamIndex (gpuStream_t const& s) {
        m_max_stream_index = std::max(m_max_stream_index,Gpu::Device::streamIndex(s));
    }

private:
    int m_max_blocks;
    int m_max_stream_index = 0;
    Type* m_host_tuple = nullptr;
    Type* m_device_tuple = nullptr;
    GpuArray<int,AMREX_GPU_MAX_STREAMS> m_nblocks;
    std::function<Type()> m_fn_value;
};

/// \cond DOXYGEN_IGNORE
namespace Reduce::detail {

    // call_f_intvect_box

    template <typename F, int dim>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_box (F const& f, IntVectND<dim> iv, IndexTypeND<dim>) noexcept ->
        decltype(amrex::detail::call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv))
    {
        return amrex::detail::call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv);
    }

    template <typename F, int dim>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_box (F const& f, IntVectND<dim> iv, IndexTypeND<dim> t) noexcept ->
        decltype(f(BoxND<dim>(iv, iv, t)))
    {
        return f(BoxND<dim>(iv, iv, t));
    }

    // call_f_intvect_n

    template <typename F, typename T, int dim>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_n (F const& f, IntVectND<dim> iv, T n) noexcept ->
        decltype(amrex::detail::call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n))
    {
        return amrex::detail::call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n);
    }

    // mf_call_f

    /// \cond DOXYGEN_IGNORE
    struct iterate_box {};
    struct iterate_box_comp {};
    /// \endcond

    template <typename I, typename F, typename T, typename... Ps,
              std::enable_if_t<std::is_same_v<iterate_box,I>,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void mf_call_f (F const& f, int ibox, int i, int j, int k, int, T& r) noexcept
    {
        auto const& pr = f(ibox,i,j,k);
        Reduce::detail::for_each_local<0, T, Ps...>(r, pr);
    }

    template <typename I, typename F, typename T, typename... Ps,
              std::enable_if_t<std::is_same_v<iterate_box_comp,I>,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void mf_call_f (F const& f, int ibox, int i, int j, int k, int ncomp, T& r) noexcept
    {
        for (int n = 0; n < ncomp; ++n) {
            auto const& pr = f(ibox,i,j,k,n);
            Reduce::detail::for_each_local<0, T, Ps...>(r, pr);
        }
    }
}
/// \endcond

/// \ingroup amrex_collectives
template <typename... Ps>
class ReduceOps
{
public:

    /// \cond DOXYGEN_IGNORE
    // This is public for CUDA
    template <typename I, typename MF, typename D, typename F>
    void eval_mf (I, MF const& mf, IntVect const& nghost, int ncomp, D& reduce_data, F const& f)
    {
        using ReduceTuple = typename D::Type;
        const int nboxes = mf.local_size();
        if (nboxes > 0) {
            auto const& parforinfo = mf.getParForInfo(nghost);
            auto nblocks_per_box = parforinfo.getNBlocksPerBox(AMREX_GPU_MAX_THREADS);
            AMREX_ASSERT(Long(nblocks_per_box)*Long(nboxes) < Long(std::numeric_limits<int>::max()));
            const int nblocks = nblocks_per_box * nboxes;
            const BoxIndexer* dp_boxes = parforinfo.getBoxes();

            auto const& stream = Gpu::gpuStream();
            auto pdst = reduce_data.devicePtr(stream);
            int nblocks_ec = std::min(nblocks, reduce_data.maxBlocks());
            AMREX_ASSERT(Long(nblocks_ec)*2 <= Long(std::numeric_limits<int>::max()));
            reduce_data.nBlocks(stream) = nblocks_ec;
            reduce_data.updateMaxStreamIndex(stream);

#ifdef AMREX_USE_SYCL
            // device reduce needs local(i.e., shared) memory
            constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
            amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, shared_mem_bytes, stream,
                          [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
            {
                Dim1 blockIdx {gh.blockIdx()};
                Dim1 threadIdx{gh.threadIdx()};
#else
            amrex::launch_global<AMREX_GPU_MAX_THREADS>
                <<<nblocks_ec, AMREX_GPU_MAX_THREADS, 0, stream>>>
                ([=] AMREX_GPU_DEVICE () noexcept
            {
#endif
                ReduceTuple r;
                Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(r);
                ReduceTuple& dst = pdst[blockIdx.x];
                if (threadIdx.x == 0) {
                    dst = r;
                }
                for (int iblock = blockIdx.x; iblock < nblocks; iblock += nblocks_ec) {
                    int ibox = iblock / nblocks_per_box;
                    auto icell = std::uint64_t(iblock-ibox*nblocks_per_box)*AMREX_GPU_MAX_THREADS + threadIdx.x;

                    BoxIndexer const& indexer = dp_boxes[ibox];
                    if (icell < indexer.numPts()) {
                        auto [i, j, k] = indexer(icell);
                        Reduce::detail::mf_call_f<I, F, ReduceTuple, Ps...>
                            (f, ibox, i, j, k, ncomp, r);
                    }
                }
#ifdef AMREX_USE_SYCL
                Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r, gh);
#else
                Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r);
#endif
            });
        }
    }
    /// \endcond

    template <typename MF, typename D, typename F>
    std::enable_if_t<IsFabArray<MF>::value
#ifndef AMREX_USE_CUDA
                     && IsCallable<F, int, int, int, int>::value
#endif
                    >
    eval (MF const& mf, IntVect const& nghost, D& reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        const int nboxes = mf.local_size();
        if (nboxes == 0) {
            return;
        } else if (!mf.isFusingCandidate()) {
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& b = amrex::grow(mfi.validbox(), nghost);
                const int li = mfi.LocalIndex();
                this->eval(b, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
                {
                    return f(li, i, j, k);
                });
            }
        } else {
            eval_mf(Reduce::detail::iterate_box{},
                    mf, nghost, 0, reduce_data, std::forward<F>(f));
        }
    }

    template <typename MF, typename D, typename F>
    std::enable_if_t<IsFabArray<MF>::value
#ifndef AMREX_USE_CUDA
                     && IsCallable<F, int, int, int, int, int>::value
#endif
                     >
    eval (MF const& mf, IntVect const& nghost, int ncomp, D& reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;

        const int nboxes = mf.local_size();

        if (nboxes == 0) {
            return;
        } else if (!mf.isFusingCandidate()) {
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                Box const& b = amrex::grow(mfi.validbox(), nghost);
                const int li = mfi.LocalIndex();
                this->eval(b, ncomp, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept -> ReduceTuple
                {
                    return f(li, i, j, k, n);
                });
            }
        } else {
            eval_mf(Reduce::detail::iterate_box_comp{},
                    mf, nghost, ncomp, reduce_data, std::forward<F>(f));
        }
    }

    template <typename D, typename F, int dim>
    void eval (BoxND<dim> const& box, D & reduce_data, F const& f)
    {
        using ReduceTuple = typename D::Type;
        auto const& stream = Gpu::gpuStream();
        auto dp = reduce_data.devicePtr(stream);
        int& nblocks = reduce_data.nBlocks(stream);
        const BoxIndexerND<dim> indexer(box);
        IndexTypeND<dim> ixtype = box.ixType();
        constexpr int nitems_per_thread = 4;
        Long nblocks_ec = (box.numPts() + nitems_per_thread*AMREX_GPU_MAX_THREADS-1)
            / (nitems_per_thread*AMREX_GPU_MAX_THREADS);
        nblocks_ec = std::min<Long>(nblocks_ec, reduce_data.maxBlocks());
        reduce_data.updateMaxStreamIndex(stream);
#ifdef AMREX_USE_SYCL
        // device reduce needs local(i.e., shared) memory
        constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
        amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, shared_mem_bytes, stream,
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
        {
            Dim1 blockIdx {gh.blockIdx()};
            Dim1 threadIdx{gh.threadIdx()};
            Dim1 gridDim  {gh.gridDim()};
#else
        amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, 0, stream,
        [=] AMREX_GPU_DEVICE () noexcept
        {
#endif
            ReduceTuple r;
            Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(r);
            ReduceTuple& dst = *(dp+blockIdx.x);
            if (threadIdx.x == 0 && static_cast<int>(blockIdx.x) >= nblocks) {
                dst = r;
            }
            for (std::uint64_t icell = std::uint64_t(AMREX_GPU_MAX_THREADS)*blockIdx.x+threadIdx.x,
                 stride = std::uint64_t(AMREX_GPU_MAX_THREADS)*gridDim.x;
                 icell < indexer.numPts();
                 icell += stride)
            {
                auto iv = indexer.intVect(icell);
                auto pr = Reduce::detail::call_f_intvect_box(f, iv, ixtype);
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
            }
#ifdef AMREX_USE_SYCL
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r, gh);
#else
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r);
#endif
        });
        nblocks = std::max(nblocks, static_cast<int>(nblocks_ec));
    }

    template <typename N, typename D, typename F, int dim,
              typename M=std::enable_if_t<std::is_integral_v<N>> >
    void eval (BoxND<dim> const& box, N ncomp, D & reduce_data, F const& f)
    {
        using ReduceTuple = typename D::Type;
        auto const& stream = Gpu::gpuStream();
        auto dp = reduce_data.devicePtr(stream);
        int& nblocks = reduce_data.nBlocks(stream);
        const BoxIndexerND<dim> indexer(box);
        constexpr int nitems_per_thread = 4;
        Long nblocks_ec = (box.numPts() + nitems_per_thread*AMREX_GPU_MAX_THREADS-1)
            / (nitems_per_thread*AMREX_GPU_MAX_THREADS);
        nblocks_ec = std::min<Long>(nblocks_ec, reduce_data.maxBlocks());
        reduce_data.updateMaxStreamIndex(stream);
#ifdef AMREX_USE_SYCL
        // device reduce needs local(i.e., shared) memory
        constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
        amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, shared_mem_bytes, stream,
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
        {
            Dim1 blockIdx {gh.blockIdx()};
            Dim1 threadIdx{gh.threadIdx()};
            Dim1 gridDim  {gh.gridDim()};
#else
        amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, 0, stream,
        [=] AMREX_GPU_DEVICE () noexcept
        {
#endif
            ReduceTuple r;
            Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(r);
            ReduceTuple& dst = *(dp+blockIdx.x);
            if (threadIdx.x == 0 && static_cast<int>(blockIdx.x) >= nblocks) {
                dst = r;
            }
            for (std::uint64_t icell = std::uint64_t(AMREX_GPU_MAX_THREADS)*blockIdx.x+threadIdx.x,
                 stride = std::uint64_t(AMREX_GPU_MAX_THREADS)*gridDim.x;
                 icell < indexer.numPts();
                 icell += stride)
            {
                auto iv = indexer.intVect(icell);
                for (N n = 0; n < ncomp; ++n) {
                    auto pr = Reduce::detail::call_f_intvect_n(f, iv, n);
                    Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
                }
            }
#ifdef AMREX_USE_SYCL
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r, gh);
#else
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r);
#endif
        });
        nblocks = std::max(nblocks, static_cast<int>(nblocks_ec));
    }

    template <typename N, typename D, typename F,
              typename M=std::enable_if_t<std::is_integral_v<N>> >
    void eval (N n, D & reduce_data, F const& f)
    {
        if (n <= 0) { return; }
        using ReduceTuple = typename D::Type;
        auto const& stream = Gpu::gpuStream();
        auto dp = reduce_data.devicePtr(stream);
        int& nblocks = reduce_data.nBlocks(stream);
        constexpr int nitems_per_thread = 4;
        int nblocks_ec = (n + nitems_per_thread*AMREX_GPU_MAX_THREADS-1)
            / (nitems_per_thread*AMREX_GPU_MAX_THREADS);
        nblocks_ec = std::min(nblocks_ec, reduce_data.maxBlocks());
        reduce_data.updateMaxStreamIndex(stream);
#ifdef AMREX_USE_SYCL
        // device reduce needs local(i.e., shared) memory
        constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
        amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, shared_mem_bytes, stream,
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
        {
            Dim1 blockIdx {gh.blockIdx()};
            Dim1 threadIdx{gh.threadIdx()};
            Dim1 gridDim  {gh.gridDim()};
#else
        amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, 0, stream,
        [=] AMREX_GPU_DEVICE () noexcept
        {
#endif
            ReduceTuple r;
            Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(r);
            ReduceTuple& dst = *(dp+blockIdx.x);
            if (threadIdx.x == 0 && static_cast<int>(blockIdx.x) >= nblocks) {
                dst = r;
            }
            for (N i = N(AMREX_GPU_MAX_THREADS)*blockIdx.x+threadIdx.x,
                 stride = N(AMREX_GPU_MAX_THREADS)*gridDim.x;
                 i < n;
                 i += stride)
            {
                auto pr = f(i);
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r,pr);
            }
#ifdef AMREX_USE_SYCL
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r, gh);
#else
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r);
#endif
        });
        nblocks = amrex::max(nblocks, nblocks_ec);
    }

    template <typename D>
    typename D::Type value (D & reduce_data)
    {
        auto hp = reduce_data.hostPtr();

        if (m_result_is_ready) {
            return *hp;
        }

        using ReduceTuple = typename D::Type;
        auto const& stream = Gpu::gpuStream();
        auto dp = reduce_data.devicePtr();
        auto const& nblocks = reduce_data.nBlocks();
#if defined(AMREX_USE_SYCL)
        if (reduce_data.maxStreamIndex() == 0 && nblocks[0] <= 4096) {
            const int N = nblocks[0];
            if (N == 0) {
                Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(*hp);
            } else {
                Gpu::PinnedVector<ReduceTuple> tmp(N);
                Gpu::dtoh_memcpy_async(tmp.data(), dp, sizeof(ReduceTuple)*N);
                Gpu::streamSynchronize();
                for (int i = 1; i < N; ++i) {
                    Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(tmp[0], tmp[i]);
                }
                *hp = tmp[0];
            }
        } else
#endif
        {
            int maxblocks = reduce_data.maxBlocks();
#ifdef AMREX_USE_SYCL
            // device reduce needs local(i.e., shared) memory
            constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
#ifndef AMREX_NO_SYCL_REDUCE_WORKAROUND
            // xxxxx SYCL todo: reduce bug workaround
            Gpu::DeviceVector<ReduceTuple> dtmp(1);
            auto presult = dtmp.data();
#else
            auto presult = hp;
#endif
            amrex::launch<AMREX_GPU_MAX_THREADS>(1, shared_mem_bytes, stream,
            [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
            {
                ReduceTuple r;
                Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(r);
                ReduceTuple dst = r;
                for (int istream = 0, nstreams = nblocks.size(); istream < nstreams; ++istream) {
                    auto dp_stream = dp+istream*maxblocks;
                    for (int i = gh.item->get_global_id(0), stride = gh.item->get_global_range(0);
                         i < nblocks[istream]; i += stride) {
                        Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, dp_stream[i]);
                    }
                }
                Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r, gh);
                if (gh.threadIdx() == 0) { *presult = dst; }
            });
#ifndef AMREX_NO_SYCL_REDUCE_WORKAROUND
            Gpu::dtoh_memcpy_async(hp, dtmp.data(), sizeof(ReduceTuple));
#endif
#else
            amrex::launch<AMREX_GPU_MAX_THREADS>(1, 0, stream,
            [=] AMREX_GPU_DEVICE () noexcept
            {
                ReduceTuple r;
                Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(r);
                ReduceTuple dst = r;
                for (int istream = 0, nstreams = nblocks.size(); istream < nstreams; ++istream) {
                    auto dp_stream = dp+istream*maxblocks;
                    for (int i = threadIdx.x; i < nblocks[istream]; i += AMREX_GPU_MAX_THREADS) {
                        Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, dp_stream[i]);
                    }
                }
                Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r);
                if (threadIdx.x == 0) { *hp = dst; }
            });
#endif
            Gpu::streamSynchronize();
        }

        m_result_is_ready = true;
        return *hp;
    }

private:
    bool m_result_is_ready = false;

public:
    void resetResultReadiness () { m_result_is_ready = false; }
};

namespace Reduce {

/// \ingroup amrex_collectives
template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
T Sum (N n, T const* v, T init_val = 0)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {v[i]}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return amrex::get<0>(hv) + init_val;
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral_v<N>> >
T Sum (N n, F const& f, T init_val = 0)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {f(i)}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return amrex::get<0>(hv) + init_val;
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
T Min (N n, T const* v, T init_val = std::numeric_limits<T>::max())
{
    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {v[i]}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return std::min(amrex::get<0>(hv),init_val);
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral_v<N>> >
T Min (N n, F const& f, T init_val = std::numeric_limits<T>::max())
{
    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {f(i)}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return std::min(amrex::get<0>(hv),init_val);
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
T Max (N n, T const* v, T init_val = std::numeric_limits<T>::lowest())
{
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {v[i]}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return std::max(amrex::get<0>(hv),init_val);
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral_v<N>> >
T Max (N n, F const& f, T init_val = std::numeric_limits<T>::lowest())
{
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {f(i)}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return std::max(amrex::get<0>(hv),init_val);
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
std::pair<T,T> MinMax (N n, T const* v)
{
    ReduceOps<ReduceOpMin,ReduceOpMax> reduce_op;
    ReduceData<T,T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple {
            return {v[i],v[i]};
        });
    auto hv = reduce_data.value(reduce_op);
    return std::make_pair(amrex::get<0>(hv), amrex::get<1>(hv));
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral_v<N>> >
std::pair<T,T> MinMax (N n, F const& f)
{
    ReduceOps<ReduceOpMin,ReduceOpMax> reduce_op;
    ReduceData<T,T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple {
            T tmp = f(i);
            return {tmp,tmp};
        });
    auto hv = reduce_data.value(reduce_op);
    return std::make_pair(amrex::get<0>(hv), amrex::get<1>(hv));
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename P, typename M=std::enable_if_t<std::is_integral_v<N>> >
bool AnyOf (N n, T const* v, P const& pred)
{
    Gpu::LaunchSafeGuard lsg(true);
    Gpu::DeviceScalar<int> ds(0);
    int* dp = ds.dataPtr();
    auto ec = Gpu::ExecutionConfig(n);
    ec.numBlocks.x = std::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch());

#ifdef AMREX_USE_SYCL
    const int num_ints = std::max(Gpu::Device::warp_size, int(ec.numThreads.x)/Gpu::Device::warp_size) + 1;
    const std::size_t shared_mem_bytes = num_ints*sizeof(int);
    amrex::launch<AMREX_GPU_MAX_THREADS>(ec.numBlocks.x, shared_mem_bytes, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept {
        int* has_any = &(static_cast<int*>(gh.sharedMemory())[num_ints-1]);
        if (gh.threadIdx() == 0) { *has_any = *dp; }
        gh.sharedBarrier();

        if (!(*has_any))
        {
            int r = false;
            for (N i = AMREX_GPU_MAX_THREADS*gh.blockIdx()+gh.threadIdx(), stride = AMREX_GPU_MAX_THREADS*gh.gridDim();
                 i < n && !r; i += stride)
            {
                r = pred(v[i]) ? 1 : 0;
            }

            r = Gpu::blockReduce<Gpu::Device::warp_size>
                (r, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::Plus<int> >(), 0, gh);
            if (gh.threadIdx() == 0 && r) { *dp = 1; }
        }
    });
#else
    amrex::launch<AMREX_GPU_MAX_THREADS>(ec.numBlocks.x, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        __shared__ int has_any;
        if (threadIdx.x == 0) { has_any = *dp; }
        __syncthreads();

        if (!has_any)
        {
            int r = false;
            for (N i = AMREX_GPU_MAX_THREADS*blockIdx.x+threadIdx.x, stride = AMREX_GPU_MAX_THREADS*gridDim.x;
                 i < n && !r; i += stride)
            {
                r = pred(v[i]) ? 1 : 0;
            }
            r = Gpu::blockReduce<Gpu::Device::warp_size>
                (r, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::Plus<int> >(), 0);
            if (threadIdx.x == 0 && r) *dp = 1;
        }
    });
#endif
    return ds.dataValue();
}

/// \ingroup amrex_collectives
template <typename P, int dim>
bool AnyOf (BoxND<dim> const& box, P const& pred)
{
    Gpu::LaunchSafeGuard lsg(true);
    Gpu::DeviceScalar<int> ds(0);
    int* dp = ds.dataPtr();
    const BoxIndexerND<dim> indexer(box);
    auto ec = Gpu::ExecutionConfig(box.numPts());
    ec.numBlocks.x = std::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch());

#ifdef AMREX_USE_SYCL
    const int num_ints = std::max(Gpu::Device::warp_size, int(ec.numThreads.x)/Gpu::Device::warp_size) + 1;
    const std::size_t shared_mem_bytes = num_ints*sizeof(int);
    amrex::launch<AMREX_GPU_MAX_THREADS>(ec.numBlocks.x, shared_mem_bytes, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept {
        int* has_any = &(static_cast<int*>(gh.sharedMemory())[num_ints-1]);
        if (gh.threadIdx() == 0) { *has_any = *dp; }
        gh.sharedBarrier();

        if (!(*has_any))
        {
            int r = false;
            for (std::uint64_t icell = std::uint64_t(AMREX_GPU_MAX_THREADS)*gh.blockIdx()+gh.threadIdx(),
                 stride = std::uint64_t(AMREX_GPU_MAX_THREADS)*gh.gridDim();
                 icell < indexer.numPts() && !r;
                 icell += stride)
            {
                auto iv = indexer.intVect(icell);
                r = amrex::detail::call_f_intvect(pred, iv) ? 1 : 0;
            }
            r = Gpu::blockReduce<Gpu::Device::warp_size>
                (r, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::Plus<int> >(), 0, gh);
            if (gh.threadIdx() == 0 && r) { *dp = 1; }
        }
    });
#else
    AMREX_LAUNCH_KERNEL(AMREX_GPU_MAX_THREADS, ec.numBlocks, ec.numThreads, 0,
                        Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        __shared__ int has_any;
        if (threadIdx.x == 0) { has_any = *dp; }
        __syncthreads();

        if (!has_any)
        {
            int r = false;
            for (std::uint64_t icell = std::uint64_t(AMREX_GPU_MAX_THREADS)*blockIdx.x+threadIdx.x,
                 stride = std::uint64_t(AMREX_GPU_MAX_THREADS)*gridDim.x;
                 icell < indexer.numPts() && !r;
                 icell += stride)
            {
                auto iv = indexer.intVect(icell);
                r = amrex::detail::call_f_intvect(pred, iv) ? 1 : 0;
            }
            r = Gpu::blockReduce<Gpu::Device::warp_size>
                (r, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::Plus<int> >(), 0);
            if (threadIdx.x == 0 && r) *dp = 1;
        }
    });
#endif
    return ds.dataValue();
}

}

#else

template <typename... Ts>
class ReduceData
{
public:
    using Type = GpuTuple<Ts...>;

    template <typename... Ps>
    explicit ReduceData (ReduceOps<Ps...>& reduce_op)
        : m_tuple(OpenMP::in_parallel() ? 1 : OpenMP::get_max_threads()),
          m_fn_value([&reduce_op,this] () -> Type { return this->value(reduce_op); })
    {
        reduce_op.resetResultReadiness();
        for (auto& t : m_tuple) {
            Reduce::detail::for_each_init<0, Type, Ps...>(t);
        }
    }

    ~ReduceData () = default;
    ReduceData (ReduceData<Ts...> const&) = delete;
    ReduceData (ReduceData<Ts...> &&) = delete;
    void operator= (ReduceData<Ts...> const&) = delete;
    void operator= (ReduceData<Ts...> &&) = delete;

    Type value () { return m_fn_value(); }

    template <typename... Ps>
    Type value (ReduceOps<Ps...>& reduce_op)
    {
        return reduce_op.value(*this);
    }

    Vector<Type>& reference () { return m_tuple; }

    Type& reference (int tid)
    {
        if (m_tuple.size() == 1) {
            // No OpenMP or already inside OpenMP parallel when reduce_data is constructed
            return m_tuple[0];
        } else {
            return m_tuple[tid];
        }
    }

private:
    Vector<Type> m_tuple;
    std::function<Type()> m_fn_value;
};

namespace Reduce::detail {

    // call_f_intvect

    template <typename F, int dim>
    AMREX_FORCE_INLINE
    auto call_f_intvect (F const& f, IntVectND<dim> iv) noexcept ->
        decltype(amrex::detail::call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv))
    {
        return amrex::detail::call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv);
    }

    // call_f_intvect_n

    template <typename F, typename T, int dim>
    AMREX_FORCE_INLINE
    auto call_f_intvect_n (F const& f, IntVectND<dim> iv, T n) noexcept ->
        decltype(amrex::detail::call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n))
    {
        return amrex::detail::call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n);
    }
}

template <typename... Ps>
class ReduceOps
{
private:

    // call_f_box

    template <typename D, typename F, int dim>
    AMREX_FORCE_INLINE
    static auto call_f_box (BoxND<dim> const& box, typename D::Type & r, F const& f)
        noexcept -> std::enable_if_t<std::is_same_v<std::decay_t<decltype(
            Reduce::detail::call_f_intvect(f, IntVectND<dim>(0))
        )>, typename D::Type>>
    {
        using ReduceTuple = typename D::Type;
        For(box,
            [&] (IntVectND<dim> iv) {
                auto pr = Reduce::detail::call_f_intvect(f, iv);
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
            });
    }

    template <typename D, typename F, int dim>
    AMREX_FORCE_INLINE
    static auto call_f_box (BoxND<dim> const& box, typename D::Type & r, F const& f)
        noexcept -> std::enable_if_t<std::is_same_v<std::decay_t<decltype(f(box))>,
                                                  typename D::Type>>
    {
        using ReduceTuple = typename D::Type;
        Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, f(box));
    }

public:

    template <typename MF, typename D, typename F>
    std::enable_if_t<IsFabArray<MF>::value && IsCallable<F, int, int, int, int>::value>
    eval (MF const& mf, IntVect const& nghost, D & reduce_data, F const& f)
    {
        using ReduceTuple = typename D::Type;
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        {
            ReduceTuple rr;
            Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(rr);
            for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
                Box const& b = mfi.growntilebox(nghost);
                const int li = mfi.LocalIndex();
                const auto lo = amrex::lbound(b);
                const auto hi = amrex::ubound(b);
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(rr, f(li,i,j,k));
                }}}
            }
            Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(
                reduce_data.reference(OpenMP::get_thread_num()), rr);
        }
    }

    template <typename MF, typename D, typename F>
    std::enable_if_t<IsFabArray<MF>::value && IsCallable<F, int, int, int, int, int>::value>
    eval (MF const& mf, IntVect const& nghost, int ncomp, D & reduce_data, F const& f)
    {
        using ReduceTuple = typename D::Type;
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        {
            ReduceTuple rr;
            Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(rr);
            for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
                Box const& b = mfi.growntilebox(nghost);
                const int li = mfi.LocalIndex();
                const auto lo = amrex::lbound(b);
                const auto hi = amrex::ubound(b);
                for (int n = 0; n < ncomp; ++n) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(rr, f(li,i,j,k,n));
                }}}}
            }
            Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(
                reduce_data.reference(OpenMP::get_thread_num()), rr);
        }
    }

    template <typename D, typename F, int dim>
    void eval (BoxND<dim> const& box, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple rr;
        Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(rr);
        call_f_box<D>(box, rr, std::forward<F>(f));
        Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(
             reduce_data.reference(OpenMP::get_thread_num()), rr);
    }

    template <typename N, typename D, typename F, int dim,
              typename M=std::enable_if_t<std::is_integral_v<N>> >
    void eval (BoxND<dim> const& box, N ncomp, D & reduce_data, F const& f)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple rr;
        Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(rr);
        For(box, ncomp,
            [&] (IntVectND<dim> iv, int n) {
                auto pr = Reduce::detail::call_f_intvect_n(f, iv, n);
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(rr, pr);
            });
        Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(
             reduce_data.reference(OpenMP::get_thread_num()), rr);
    }

    template <typename N, typename D, typename F,
              typename M=std::enable_if_t<std::is_integral_v<N>> >
    void eval (N n, D & reduce_data, F const& f)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple rr;
        Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(rr);
        for (N i = 0; i < n; ++i) {
            Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(rr, f(i));
        }
        Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(
             reduce_data.reference(OpenMP::get_thread_num()), rr);
    }

    template <typename D>
    typename D::Type value (D & reduce_data)
    {
        auto& rrv = reduce_data.reference();
        if (! m_result_is_ready) {
            using ReduceTuple = typename D::Type;
            if (rrv.size() > 1) {
                for (int i = 1, N = rrv.size(); i < N; ++i) {
                    Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(rrv[0], rrv[i]);
                }
            }
            m_result_is_ready = true;
        }
        return rrv[0];
    }

    bool m_result_is_ready = false;

    void resetResultReadiness () { m_result_is_ready = false; }
};

namespace Reduce {

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral_v<N>> >
T Sum (N n, F const& f, T init_val = 0)
{
    T r = init_val;
#ifdef AMREX_USE_OMP
#pragma omp parallel for reduction(+:r)
#endif
    for (N i = 0; i < n; ++i) {
        r += f(i);
    }
    return r;
}

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
T Sum (N n, T const* v, T init_val = 0)
{
    return Sum(n, [=] (N i) -> T { return v[i]; }, init_val);
}

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral_v<N>> >
T Min (N n, F const& f, T init_val = std::numeric_limits<T>::max())
{
    T r = init_val;
#ifdef AMREX_USE_OMP
#pragma omp parallel for reduction(min:r)
#endif
    for (N i = 0; i < n; ++i) {
        r = std::min(r,f(i));
    }
    return r;
}

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
T Min (N n, T const* v, T init_val = std::numeric_limits<T>::max())
{
    return Reduce::Min(n, [=] (N i) -> T { return v[i]; }, init_val);
}

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral_v<N>> >
T Max (N n, F const& f, T init_val = std::numeric_limits<T>::lowest())
{
    T r = init_val;
#ifdef AMREX_USE_OMP
#pragma omp parallel for reduction(max:r)
#endif
    for (N i = 0; i < n; ++i) {
        r = std::max(r,f(i));
    }
    return r;
}

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
T Max (N n, T const* v, T init_val = std::numeric_limits<T>::lowest())
{
    return Reduce::Max(n, [=] (N i) -> T { return v[i]; }, init_val);
}

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral_v<N>> >
std::pair<T,T> MinMax (N n, F const& f)
{
    T r_min = std::numeric_limits<T>::max();
    T r_max = std::numeric_limits<T>::lowest();
#ifdef AMREX_USE_OMP
#pragma omp parallel for reduction(min:r_min) reduction(max:r_max)
#endif
    for (N i = 0; i < n; ++i) {
        T tmp = f(i);
        r_min = std::min(r_min,tmp);
        r_max = std::max(r_max,tmp);
    }
    return std::make_pair(r_min,r_max);
}

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
std::pair<T,T> MinMax (N n, T const* v)
{
    return Reduce::MinMax<T>(n, [=] (N i) -> T { return v[i]; });
}

template <typename T, typename N, typename P, typename M=std::enable_if_t<std::is_integral_v<N>> >
bool AnyOf (N n, T const* v, P&& pred)
{
    return std::any_of(v, v+n, std::forward<P>(pred));
}

template <typename P, int dim>
bool AnyOf (BoxND<dim> const& box, P const& pred)
{
    for (auto iv : box.iterator()) { // NOLINT(readability-use-anyofallof)
        if (Reduce::detail::call_f_intvect(pred, iv)) { return true; }
    }
    return false;
}

}

#endif

/**
 * \brief Return a GpuTuple containing the identity element for each operation in ReduceOps.
 * For example 0, +inf and -inf for ReduceOpSum, ReduceOpMin and ReduceOpMax respectively.
 */
template <typename... Ts, typename... Ps>
AMREX_GPU_HOST_DEVICE
constexpr GpuTuple<Ts...>
IdentityTuple (GpuTuple<Ts...>, ReduceOps<Ps...>) noexcept
{
    GpuTuple<Ts...> r{};
    Reduce::detail::for_each_init<0, decltype(r), Ps...>(r);
    return r;
}

/**
 * \brief Return a GpuTuple containing the identity element for each ReduceOp in TypeList.
 * For example 0, +inf and -inf for ReduceOpSum, ReduceOpMin and ReduceOpMax respectively.
 */
template <typename... Ts, typename... Ps>
AMREX_GPU_HOST_DEVICE
constexpr GpuTuple<Ts...>
IdentityTuple (GpuTuple<Ts...>, TypeList<Ps...>) noexcept
{
    GpuTuple<Ts...> r{};
    Reduce::detail::for_each_init<0, decltype(r), Ps...>(r);
    return r;
}

namespace Reduce {

/// \ingroup amrex_collectives
template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
T Sum (N n, T* v, T init_val = 0)
{
    return Sum(n, (T const*)v, init_val);
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
T Min (N n, T* v, T init_val = std::numeric_limits<T>::max())
{
    return Min(n, (T const*)v, init_val);
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
T Max (N n, T* v, T init_val = std::numeric_limits<T>::lowest())
{
    return Max(n, (T const*)v, init_val);
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename M=std::enable_if_t<std::is_integral_v<N>> >
std::pair<T,T> MinMax (N n, T* v)
{
    return MinMax(n, (T const*)v);
}

/// \ingroup amrex_collectives
template <typename T, typename N, typename P, typename M=std::enable_if_t<std::is_integral_v<N>> >
bool AnyOf (N n, T* v, P const& pred)
{
    return AnyOf(n, (T const*)v, pred);
}

}

/// \cond DOXYGEN_IGNORE
namespace reducer_detail
{
    template <typename F, std::size_t...Ns>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect (std::index_sequence<Ns...>, F const& f, int iop, IntVectND<1> iv)
        noexcept -> decltype(f(0, 0, 0, 0))
    {
        return f(iop, iv[0], 0, 0);
    }

    template <typename F, std::size_t...Ns>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect (std::index_sequence<Ns...>, F const& f, int iop, IntVectND<2> iv)
        noexcept -> decltype(f(0, 0, 0, 0))
    {
        return f(iop, iv[0], iv[1], 0);
    }

    template <typename F, int dim, std::size_t...Ns>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect (std::index_sequence<Ns...>, F const& f, int iop, IntVectND<dim> iv)
        noexcept -> decltype(f(0, iv))
    {
        return f(iop, iv);
    }

    template <typename F, int dim, std::size_t...Ns>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect (std::index_sequence<Ns...>, F const& f, int iop, IntVectND<dim> iv)
        noexcept -> decltype(f(0, iv[Ns]...))
    {
        return f(iop, iv[Ns]...);
    }

    template <typename F, std::size_t...Ns>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_n (std::index_sequence<Ns...>, F const& f, int iop, IntVectND<1> iv, int n)
        noexcept -> decltype(f(0, 0, 0, 0, 0))
    {
        return f(iop, iv[0], 0, 0, n);
    }

    template <typename F, std::size_t...Ns>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_n (std::index_sequence<Ns...>, F const& f, int iop, IntVectND<2> iv, int n)
        noexcept -> decltype(f(0, 0, 0, 0, 0))
    {
        return f(iop, iv[0], iv[1], 0, n);
    }

    template <typename F, int dim, std::size_t...Ns>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_n (std::index_sequence<Ns...>, F const& f, int iop, IntVectND<dim> iv, int n)
        noexcept -> decltype(f(0, iv, n))
    {
        return f(iop, iv, n);
    }

    template <typename F, int dim, std::size_t...Ns>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_n (std::index_sequence<Ns...>, F const& f, int iop, IntVectND<dim> iv, int n)
        noexcept -> decltype(f(0, iv[Ns]..., n))
    {
        return f(iop, iv[Ns]..., n);
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T init_value (ReduceOpType op)
    {
        T r;
        if (ReduceOpType::max == op) {
            ReduceOpMax{}.init(r);
        } else if (ReduceOpType::min == op) {
            ReduceOpMin{}.init(r);
        } else {
            ReduceOpSum{}.init(r);
        }
        return r;
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void local_update (ReduceOpType op, T& dst, T const& src)
    {
        if (ReduceOpType::max == op) {
            ReduceOpMax{}.local_update(dst, src);
        } else if (ReduceOpType::min == op) {
            ReduceOpMin{}.local_update(dst, src);
        } else {
            ReduceOpSum{}.local_update(dst, src);
        }
    }

#ifdef AMREX_USE_GPU

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void parallel_update (ReduceOpType op, T& dst, T const& src
#ifdef AMREX_USE_SYCL
                          , Gpu::Handler const& gh
#endif
        )
    {
        if (ReduceOpType::max == op) {
            ReduceOpMax{}.parallel_update(dst, src
#ifdef AMREX_USE_SYCL
                                          , gh
#endif
                );
        } else if (ReduceOpType::min == op) {
            ReduceOpMin{}.parallel_update(dst, src
#ifdef AMREX_USE_SYCL
                                          , gh
#endif
                );
        } else {
            ReduceOpSum{}.parallel_update(dst, src
#ifdef AMREX_USE_SYCL
                                          , gh
#endif
                );
        }
    }

    template <typename OP, typename T, int dim, typename F>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void block_update (int iop, OP op, T* dst, std::uint64_t icell,
                       std::uint64_t stride, BoxIndexerND<dim> const& indexer,
                       F const& f
#ifdef AMREX_USE_SYCL
                       , Gpu::Handler const& gh
#endif
        )
    {
        T r;
        op.init(r);
        for ( ; icell < indexer.numPts(); icell += stride) {
            auto pr = call_f_intvect(std::make_index_sequence<dim>(),
                                     f, iop, indexer.intVect(icell));
            op.local_update(r, pr);
        }
        op.parallel_update(*dst, r
#ifdef AMREX_USE_SYCL
                          , gh
#endif
            );
    }

    template <typename I, typename T, int dim, typename F,
              std::enable_if_t<std::is_same_v<Reduce::detail::iterate_box,I>,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void block_update (int iop, ReduceOpType op, T* dst, std::uint64_t icell,
                       std::uint64_t stride, BoxIndexerND<dim> const& indexer,
                       int, F const& f
#ifdef AMREX_USE_SYCL
                       , Gpu::Handler const& gh
#endif
        )
    {
        if (ReduceOpType::max == op) {
            block_update(iop, ReduceOpMax{}, dst, icell, stride, indexer, f
#ifdef AMREX_USE_SYCL
                         , gh
#endif
                );
        } else if (ReduceOpType::min == op) {
            block_update(iop, ReduceOpMin{}, dst, icell, stride, indexer, f
#ifdef AMREX_USE_SYCL
                         , gh
#endif
                );
        } else {
            block_update(iop, ReduceOpSum{}, dst, icell, stride, indexer, f
#ifdef AMREX_USE_SYCL
                         , gh
#endif
                );
        }
    }

    template <typename OP, typename T, int dim, typename F>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void block_update (int iop, OP op, T* dst, std::uint64_t icell,
                       std::uint64_t stride, BoxIndexerND<dim> const& indexer,
                       int ncomp, F const& f
#ifdef AMREX_USE_SYCL
                       , Gpu::Handler const& gh
#endif
        )
    {
        T r;
        op.init(r);
        for ( ; icell < indexer.numPts(); icell += stride) {
            auto iv = indexer.intVect(icell);
            for (int n = 0; n < ncomp; ++n) {
                auto pr = call_f_intvect_n(std::make_index_sequence<dim>(),
                                           f, iop, iv, n);
                op.local_update(r, pr);
            }
        }
        op.parallel_update(*dst, r
#ifdef AMREX_USE_SYCL
                          , gh
#endif
            );
    }

    template <typename I, typename T, int dim, typename F,
              std::enable_if_t<std::is_same_v<Reduce::detail::iterate_box_comp,I>,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void block_update (int iop, ReduceOpType op, T* dst, std::uint64_t icell,
                       std::uint64_t stride, BoxIndexerND<dim> const& indexer,
                       int ncomp, F const& f
#ifdef AMREX_USE_SYCL
                       , Gpu::Handler const& gh
#endif
        )
    {
        if (ReduceOpType::max == op) {
            block_update(iop, ReduceOpMax{}, dst, icell, stride, indexer, ncomp, f
#ifdef AMREX_USE_SYCL
                         , gh
#endif
                );
        } else if (ReduceOpType::min == op) {
            block_update(iop, ReduceOpMin{}, dst, icell, stride, indexer, ncomp, f
#ifdef AMREX_USE_SYCL
                         , gh
#endif
                );
        } else {
            block_update(iop, ReduceOpSum{}, dst, icell, stride, indexer, ncomp, f
#ifdef AMREX_USE_SYCL
                         , gh
#endif
                );
        }
    }

    template <typename OP, typename T, typename N, typename F>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void block_update (int iop, OP op, T* dst, N i, N n, N stride, F const& f
#ifdef AMREX_USE_SYCL
                       , Gpu::Handler const& gh
#endif
        )
    {
        T r;
        op.init(r);
        for ( ; i < n; i += stride) {
            auto pr = f(iop, i);
            op.local_update(r, pr);
        }
        op.parallel_update(*dst, r
#ifdef AMREX_USE_SYCL
                          , gh
#endif
            );
    }

    template <typename T, typename N, typename F>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void block_update (int iop, ReduceOpType op, T* dst, N i, N n, N stride, F const& f
#ifdef AMREX_USE_SYCL
                       , Gpu::Handler const& gh
#endif
        )
    {
        if (ReduceOpType::max == op) {
            block_update(iop, ReduceOpMax{}, dst, i, n, stride, f
#ifdef AMREX_USE_SYCL
                         , gh
#endif
                );
        } else if (ReduceOpType::min == op) {
            block_update(iop, ReduceOpMin{}, dst, i, n, stride, f
#ifdef AMREX_USE_SYCL
                         , gh
#endif
                );
        } else {
            block_update(iop, ReduceOpSum{}, dst, i, n, stride, f
#ifdef AMREX_USE_SYCL
                         , gh
#endif
                );
        }
    }

    template <typename OP, typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T final_update (int iop, int nops, OP op, T const* dp, int maxblocks,
                    GpuArray<int,AMREX_GPU_MAX_STREAMS> const& nblocks
#ifdef AMREX_USE_SYCL
                    , Gpu::Handler const& gh
#endif
        )
    {
        T r;
        op.init(r);
        T dst = r;
        for (int istream = 0; istream < AMREX_GPU_MAX_STREAMS; ++istream) {
            auto dp_stream = dp + maxblocks*(istream*nops + iop);
#ifdef AMREX_USE_SYCL
            for (int i = gh.threadIdx(); i < nblocks[istream]; i += int(gh.blockDim())) {
#else
            for (int i = threadIdx.x; i < nblocks[istream]; i += int(blockDim.x)) {
#endif
                op.local_update(r, dp_stream[i]);
            }
        }
        op.parallel_update(dst, r
#ifdef AMREX_USE_SYCL
                           , gh
#endif
            );
        return dst;
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T final_update (int iop, int nops, ReduceOpType op, T const* dp, int maxblocks,
                    GpuArray<int,AMREX_GPU_MAX_STREAMS> const& nblocks
#ifdef AMREX_USE_SYCL
                    , Gpu::Handler const& gh
#endif
        )
    {
        if (ReduceOpType::max == op) {
            return final_update(iop, nops, ReduceOpMax{}, dp, maxblocks, nblocks
#ifdef AMREX_USE_SYCL
                                , gh
#endif
                );
        } else if (ReduceOpType::min == op) {
            return final_update(iop, nops, ReduceOpMin{}, dp, maxblocks, nblocks
#ifdef AMREX_USE_SYCL
                                , gh
#endif
                );
        } else {
            return final_update(iop, nops, ReduceOpSum{}, dp, maxblocks, nblocks
#ifdef AMREX_USE_SYCL
                                , gh
#endif
                );
        }
    }

    template <typename I, typename T, typename F,
              std::enable_if_t<std::is_same_v<Reduce::detail::iterate_box,I>,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void mf_call_f (int iop, ReduceOpType op, T& dst, F const& f,
                    int ibox, int i, int j, int k, int)
    {
        auto pr = f(iop, ibox, i, j, k);
        local_update(op, dst, pr);
    }

    template <typename I, typename T, typename F,
              std::enable_if_t<std::is_same_v<Reduce::detail::iterate_box_comp,I>,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void mf_call_f (int iop, ReduceOpType op, T& dst, F const& f,
                    int ibox, int i, int j, int k, int ncomp)
    {
        for (int n = 0; n < ncomp; ++n) {
            auto pr = f(iop, ibox, i, j, k, n);
            local_update(op, dst, pr);
        }
    }

#else

    template <typename T>
    void reduceop_doit (ReduceOpType op, T& dst, T const& src)
    {
        if (ReduceOpType::max == op) {
            ReduceOpMax{}.local_update(dst, src);
        } else if (ReduceOpType::min == op) {
            ReduceOpMin{}.local_update(dst, src);
        } else {
            ReduceOpSum{}.local_update(dst, src);
        }
    }

#endif
}
/// \endcond

/// \ingroup amrex_collectives
/// Note that this cannot be reused.
template <typename T>
class Reducer
{
public:
    static_assert(std::is_trivially_copyable<T>(),
                  "Reducer's type must be trivially copyable");
    static_assert(std::is_trivially_destructible<T>(),
                  "Reducer's type must be trivially destructible");

    explicit Reducer (Vector<ReduceOpType> const& reduce_ops)
#ifdef AMREX_USE_GPU
        : m_max_blocks(Gpu::Device::maxBlocksPerLaunch()),
          m_host_data((T*)(The_Pinned_Arena()->alloc(sizeof(T)*reduce_ops.size()))),
          m_device_data((T*)(The_Arena()->alloc((AMREX_GPU_MAX_STREAMS)
                    * m_max_blocks * sizeof(T) * reduce_ops.size()))),
#else
        : m_data(OpenMP::in_parallel() ? 1 : OpenMP::get_max_threads()), // no support for nested OMP
#endif
          m_ops(reduce_ops.data(), reduce_ops.size())
    {
        // For now, let's not worry about land and lor, because they only
        // work on integral type.
        for (auto op : reduce_ops) {
            AMREX_ALWAYS_ASSERT(op != ReduceOpType::land &&
                                op != ReduceOpType::lor);
        }

#ifdef AMREX_USE_GPU
        for (int i = 0, N = int(m_ops.size()); i < N; ++i) {
            new (m_host_data+i) T();
        }
        m_nblocks.fill(0);
#else
        for (auto& v : m_data) {
            v.reserve(reduce_ops.size());
            for (auto op : reduce_ops) {
                v.push_back(reducer_detail::init_value<T>(op));
            }
        }
#endif
    }

    ~Reducer ()
    {
#ifdef AMREX_USE_GPU
        The_Pinned_Arena()->free(m_host_data);
        The_Arena()->free(m_device_data);
#endif
    }

    Reducer (Reducer const&) = delete;
    Reducer (Reducer &&) = delete;
    void operator= (Reducer const&) = delete;
    void operator= (Reducer &&) = delete;

    template <int dim, typename F>
    std::enable_if_t<IsCallable<F, int, int, int, int>::value ||
                     IsCallable<F, int, IntVectND<dim>>::value>
    eval (BoxND<dim> const& box, F const& f)
    {
        if (box.isEmpty()) { return; }

#ifdef AMREX_USE_GPU

        eval_box(Reduce::detail::iterate_box{}, box, 0, f);

#else

        auto* AMREX_RESTRICT dst = m_data[OpenMP::get_thread_num()].data();
        int nops = m_ops.size();
        auto const* AMREX_RESTRICT ops = m_ops.data();
        For(box, [&] (IntVectND<dim> const& iv)
        {
            for (int iop = 0; iop < nops; ++iop) {
                auto pr = reducer_detail::call_f_intvect
                    (std::make_index_sequence<dim>(), f, iop, iv);
                reducer_detail::reduceop_doit(ops[iop], dst[iop], pr);
            }
        });

#endif
    }

    template <int dim, typename F>
    std::enable_if_t<IsCallable<F, int, int, int, int, int>::value ||
                     IsCallable<F, int, IntVectND<dim>, int>::value>
    eval (BoxND<dim> const& box, int ncomp, F const& f)
    {
        if (box.isEmpty()) { return; }

#ifdef AMREX_USE_GPU

        eval_box(Reduce::detail::iterate_box_comp{}, box, ncomp, f);

#else

        auto* AMREX_RESTRICT dst = m_data[OpenMP::get_thread_num()].data();
        int nops = m_ops.size();
        auto const* AMREX_RESTRICT ops = m_ops.data();
        For(box, ncomp, [&] (IntVectND<dim> const& iv, int n)
        {
            for (int iop = 0; iop < nops; ++iop) {
                auto pr = reducer_detail::call_f_intvect_n
                    (std::make_index_sequence<dim>(), f, iop, iv, n);
                reducer_detail::reduceop_doit(ops[iop], dst[iop], pr);
            }
        });

#endif
    }

    template <typename MF, typename F>
    std::enable_if_t<IsFabArray<MF>::value &&
                     IsCallable<F, int, int, int, int, int>::value>
    eval (MF const& mf, IntVect const& nghost, F const& f)
    {
        const int nboxes = mf.local_size();
        if (nboxes == 0) {
            return;
        } else
#ifdef AMREX_USE_GPU
            if (!mf.isFusingCandidate())
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& b = mfi.growntilebox(nghost);
                int li = mfi.LocalIndex();
                this->eval(b, [=] AMREX_GPU_DEVICE (int op, int i, int j, int k)
                {
                    return f(op, li, i, j, k);
                });
            }
        }
#ifdef AMREX_USE_GPU
        else {
            eval_mf(Reduce::detail::iterate_box{}, mf, nghost, 0, f);
        }
#endif
    }

    template <typename MF, typename F>
    std::enable_if_t<IsFabArray<MF>::value &&
                     IsCallable<F, int, int, int, int, int, int>::value>
    eval (MF const& mf, IntVect const& nghost, int ncomp, F const& f)
    {
        const int nboxes = mf.local_size();
        if (nboxes == 0) {
            return;
        } else
#ifdef AMREX_USE_GPU
            if (!mf.isFusingCandidate())
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& b = mfi.growntilebox(nghost);
                int li = mfi.LocalIndex();
                this->eval(b, ncomp,
                    [=] AMREX_GPU_DEVICE (int op, int i, int j, int k, int n)
                {
                    return f(op, li, i, j, k, n);
                });
            }
        }
#ifdef AMREX_USE_GPU
        else {
            eval_mf(Reduce::detail::iterate_box_comp{}, mf, nghost, ncomp, f);
        }
#endif
    }

#ifdef AMREX_USE_GPU
/// \cond DOXYGEN_IGNORE

    template <typename I, int dim, typename F>
    void eval_box (I, BoxND<dim> const& box, int ncomp, F const& f)
    {
        if (box.isEmpty()) { return; }

        auto nops = int(m_ops.size());
        auto* ops = m_ops.data();
        const BoxIndexerND<dim> indexer(box);

        auto const& stream = Gpu::gpuStream();
        int sidx = Gpu::Device::streamIndex(stream);
        auto* dp = m_device_data + sidx*m_max_blocks*nops;
        int& nblocks = m_nblocks[sidx];
        auto maxblocks = m_max_blocks;
        constexpr int nitems_per_thread = 4;
        Long nblocks_ec = (box.numPts() + nitems_per_thread*AMREX_GPU_MAX_THREADS-1)
            / (nitems_per_thread*AMREX_GPU_MAX_THREADS);
        nblocks_ec = std::min<Long>(nblocks_ec, m_max_blocks);
        m_max_stream_index = std::max(m_max_stream_index, sidx);

#ifdef AMREX_USE_SYCL
        // device reduce needs local(i.e., shared) memory
        constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
        amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, shared_mem_bytes, stream,
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
        {
            Dim1 blockIdx {gh.blockIdx()};
            Dim1 threadIdx{gh.threadIdx()};
            Dim1 gridDim  {gh.gridDim()};
#else
        amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, 0, stream,
        [=] AMREX_GPU_DEVICE () noexcept
        {
#endif
            std::uint64_t icell = std::uint64_t(AMREX_GPU_MAX_THREADS)*blockIdx.x+threadIdx.x;
            std::uint64_t stride = std::uint64_t(AMREX_GPU_MAX_THREADS)*gridDim.x;

            for (int iop = 0; iop < nops; ++iop) {
                auto op = ops[iop];
                T* dst = dp + iop*maxblocks + blockIdx.x;
                if (threadIdx.x == 0 && static_cast<int>(blockIdx.x) >= nblocks) {
                    *dst = reducer_detail::init_value<T>(op);
                }
                reducer_detail::block_update<I>(iop, op, dst, icell, stride, indexer, ncomp, f
#ifdef AMREX_USE_SYCL
                                                , gh
#endif
                    );
            }
        });
        nblocks = amrex::max(nblocks, static_cast<int>(nblocks_ec));
    }

    template <typename I, typename MF, typename F>
    void eval_mf (I, MF const& mf, IntVect const& nghost, int ncomp, F const& f)
    {
        const int nboxes = mf.local_size();
        if (nboxes <= 0) { return; }

        auto nops = int(m_ops.size());
        auto* ops = m_ops.data();

        auto const& parforinfo = mf.getParForInfo(nghost);
        auto nblocks_per_box = parforinfo.getNBlocksPerBox(AMREX_GPU_MAX_THREADS);
        AMREX_ASSERT(Long(nblocks_per_box)*Long(nboxes) < Long(std::numeric_limits<int>::max()));
        const int nblocks = nblocks_per_box * nboxes;
        const BoxIndexer* dp_boxes = parforinfo.getBoxes();

        auto const& stream = Gpu::gpuStream();
        auto* dp = m_device_data;
        auto maxblocks = m_max_blocks;
        int nblocks_ec = std::min(nblocks, maxblocks);
        AMREX_ASSERT(Long(nblocks_ec)*2 <= Long(std::numeric_limits<int>::max()));
        m_nblocks[0] = nblocks_ec;
        m_max_stream_index = 0;

#ifdef AMREX_USE_SYCL
        // device reduce needs local(i.e., shared) memory
        constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
        amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, shared_mem_bytes, stream,
                      [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
        {
            Dim1 blockIdx {gh.blockIdx()};
            Dim1 threadIdx{gh.threadIdx()};
#else
        amrex::launch_global<AMREX_GPU_MAX_THREADS>
            <<<nblocks_ec, AMREX_GPU_MAX_THREADS, 0, stream>>>
            ([=] AMREX_GPU_DEVICE () noexcept
        {
#endif
            for (int iop = 0; iop < nops; ++iop) {
                auto op = ops[iop];
                T r = reducer_detail::init_value<T>(op);
                T* dst = dp + iop*maxblocks + blockIdx.x;
                if (threadIdx.x == 0) {
                    *dst = r;
                }
                for (int iblock = blockIdx.x; iblock < nblocks; iblock += nblocks_ec) {
                    int ibox = iblock / nblocks_per_box;
                    auto icell = std::uint64_t(iblock-ibox*nblocks_per_box)*AMREX_GPU_MAX_THREADS + threadIdx.x;

                    BoxIndexer const& indexer = dp_boxes[ibox];
                    if (icell < indexer.numPts()) {
                        auto [i, j, k] = indexer(icell);
                        reducer_detail::mf_call_f<I>(iop, op, r, f, ibox, i, j, k, ncomp);
                    }
                }

                reducer_detail::parallel_update(op, *dst, r
#ifdef AMREX_USE_SYCL
                                                , gh
#endif
                                                );
            }
        });
    }

/// \endcond
#endif

    template <typename N, typename F>
    std::enable_if_t<IsCallable<F, int, int>::value>
    eval (N n, F const& f)
    {
        if (n <= 0) { return; }

#ifdef AMREX_USE_GPU

        auto nops = int(m_ops.size());
        auto* ops = m_ops.data();

        auto const& stream = Gpu::gpuStream();
        int sidx = Gpu::Device::streamIndex(stream);
        auto* dp = m_device_data + sidx*m_max_blocks*nops;
        int& nblocks = m_nblocks[sidx];
        auto maxblocks = m_max_blocks;
        constexpr int nitems_per_thread = 4;
        int nblocks_ec = (n + nitems_per_thread*AMREX_GPU_MAX_THREADS-1)
            / (nitems_per_thread*AMREX_GPU_MAX_THREADS);
        nblocks_ec = std::min(nblocks_ec, m_max_blocks);
        m_max_stream_index = std::max(m_max_stream_index, sidx);

#ifdef AMREX_USE_SYCL
        // device reduce needs local(i.e., shared) memory
        constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
        amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, shared_mem_bytes, stream,
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
        {
            Dim1 blockIdx {gh.blockIdx()};
            Dim1 threadIdx{gh.threadIdx()};
            Dim1 gridDim  {gh.gridDim()};
#else
        amrex::launch<AMREX_GPU_MAX_THREADS>(nblocks_ec, 0, stream,
        [=] AMREX_GPU_DEVICE () noexcept
        {
#endif
            N i = N(AMREX_GPU_MAX_THREADS)*blockIdx.x+threadIdx.x;
            N stride = N(AMREX_GPU_MAX_THREADS)*gridDim.x;

            for (int iop = 0; iop < nops; ++iop) {
                auto op = ops[iop];
                T* dst = dp + iop*maxblocks + blockIdx.x;
                if (threadIdx.x == 0 && static_cast<int>(blockIdx.x) >= nblocks) {
                    *dst = reducer_detail::init_value<T>(op);
                }
#ifdef AMREX_USE_SYCL
                reducer_detail::block_update(iop, op, dst, i, n, stride, f, gh);
#else
                reducer_detail::block_update(iop, op, dst, i, n, stride, f);
#endif
            }
        });
        nblocks = amrex::max(nblocks, nblocks_ec);

#else

        auto* AMREX_RESTRICT dst = m_data[OpenMP::get_thread_num()].data();
        int nops = m_ops.size();
        auto const* AMREX_RESTRICT ops = m_ops.data();
        for (N i = 0; i < n; ++i) {
            for (int iop = 0; iop < nops; ++iop) {
                auto pr = f(iop, i);
                reducer_detail::reduceop_doit(ops[iop], dst[iop], pr);
            }
        }

#endif
    }

    Vector<T> getResults () const
    {
#ifdef AMREX_USE_GPU

        if (! m_results_ready)
        {
            auto nops = int(m_ops.size());
            auto* ops = m_ops.data();

            auto const& stream = Gpu::gpuStream();
            auto const* dp = m_device_data;
            auto* hp = m_host_data;
            auto maxblocks = m_max_blocks;
            auto nblocks = m_nblocks;

#ifdef AMREX_USE_SYCL
            // device reduce needs local(i.e., shared) memory
            constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
#ifndef AMREX_NO_SYCL_REDUCE_WORKAROUND
            // xxxxx SYCL todo: reduce bug workaround
            Gpu::DeviceVector<T> dtmp(m_ops.size());
            auto* presult = dtmp.data();
#else
            auto* presult = hp;
#endif
#else
            constexpr std::size_t shared_mem_bytes = 0;
#endif
            amrex::launch<AMREX_GPU_MAX_THREADS>(1, shared_mem_bytes, stream,
            [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_SYCL
                Gpu::Handler const& gh
#endif
                ) noexcept
            {
                for (int iop = 0; iop < nops; ++iop) {
                    auto r = reducer_detail::final_update<T>(iop, nops, ops[iop], dp, maxblocks, nblocks
#ifdef AMREX_USE_SYCL
                                                             , gh
#endif
                        );
#ifdef AMREX_USE_SYCL
                    if (gh.threadIdx() == 0) { presult[iop] = r; }
#else
                    if (threadIdx.x == 0) { hp[iop] = r; }
#endif
                }
            });
#if defined(AMREX_USE_SYCL) && !defined(AMREX_NO_SYCL_REDUCE_WORKAROUND)
            Gpu::dtoh_memcpy_async(hp, dtmp.data(), sizeof(T)*m_ops.size());
#endif
            Gpu::streamSynchronize();

            m_results_ready = true;
        }

        return Vector<T>(m_host_data, m_host_data+m_ops.size());

#else
        if (! m_results_ready) {
            for (int i = 1, N = int(m_data.size()); i < N; ++i) {
                for (int iop = 0, nops = int(m_ops.size()); iop < nops; ++iop) {
                    reducer_detail::reduceop_doit(m_ops[iop], m_data[0][iop], m_data[i][iop]);
                }
            }
            m_results_ready = true;
        }
        return m_data[0];
#endif
    }

private:

#ifdef AMREX_USE_GPU
    int m_max_blocks;
    int m_max_stream_index = 0;
    T* m_host_data = nullptr;
    T* m_device_data = nullptr;
    GpuArray<int,AMREX_GPU_MAX_STREAMS> m_nblocks;
#else
    mutable Vector<Vector<T>> m_data;
#endif

    Gpu::Buffer<ReduceOpType> m_ops;
    mutable bool m_results_ready = false;
};

}

#endif
