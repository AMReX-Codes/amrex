#ifndef AMREX_REDUCE_H_
#define AMREX_REDUCE_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_OpenMP.H>
#include <AMReX_MFIter.H>
#include <AMReX_ValLocPair.H>

#include <algorithm>
#include <functional>
#include <limits>

namespace amrex {

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
    void for_each_init (T& t)
    {
        P().init(amrex::get<I>(t));
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void for_each_init (T& t)
    {
        P().init(amrex::get<I>(t));
        for_each_init<I+1,T,P1,Ps...>(t);
    }
}

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

struct ReduceOpLogicalAnd
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral<T>::value>
    parallel_update (T& d, T s, Gpu::Handler const& h) const noexcept {
        T r = Gpu::blockReduceLogicalAnd(s,h);
        if (h.threadIdx() == 0) { d = d && r; }
    }
#else
    template <typename T, int MT=AMREX_GPU_MAX_THREADS>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral<T>::value>
    parallel_update (T& d, T s) const noexcept {
        T r = Gpu::blockReduceLogicalAnd<MT>(s);
        if (threadIdx.x == 0) { d = d && r; }
    }
#endif
#endif

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral<T>::value>
    local_update (T& d, T s) const noexcept { d = d && s; }

    template <typename T>
    constexpr std::enable_if_t<std::is_integral<T>::value>
    init (T& t) const noexcept { t = true; }
};

struct ReduceOpLogicalOr
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral<T>::value>
    parallel_update (T& d, T s, Gpu::Handler const& h) const noexcept {
        T r = Gpu::blockReduceLogicalOr(s,h);
        if (h.threadIdx() == 0) { d = d || r; }
    }
#else
    template <typename T, int MT=AMREX_GPU_MAX_THREADS>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral<T>::value>
    parallel_update (T& d, T s) const noexcept {
        T r = Gpu::blockReduceLogicalOr<MT>(s);
        if (threadIdx.x == 0) { d = d || r; }
    }
#endif
#endif

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::enable_if_t<std::is_integral<T>::value>
    local_update (T& d, T s) const noexcept { d = d || s; }

    template <typename T>
    constexpr std::enable_if_t<std::is_integral<T>::value>
    init (T& t) const noexcept { t = false; }
};

template <typename... Ps> class ReduceOps;

#ifdef AMREX_USE_GPU

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

namespace Reduce::detail {
    template <typename F>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f (F const& f, int i, int j, int k, IndexType)
        noexcept -> decltype(f(0,0,0))
    {
        return f(i,j,k);
    }

    template <typename F>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f (F const& f, int i, int j, int k, IndexType t)
        noexcept -> decltype(f(Box()))
    {
        AMREX_D_PICK(amrex::ignore_unused(j,k),amrex::ignore_unused(k),(void)0);
        return f(Box(IntVect(AMREX_D_DECL(i,j,k)),
                     IntVect(AMREX_D_DECL(i,j,k)),
                     t));
    }

    struct iterate_box {};
    struct iterate_box_comp {};

    template <typename I, typename F, typename T, typename... Ps,
              std::enable_if_t<std::is_same<iterate_box,I>::value,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void mf_call_f (F const& f, int ibox, int i, int j, int k, int, T& r) noexcept
    {
        auto const& pr = f(ibox,i,j,k);
        Reduce::detail::for_each_local<0, T, Ps...>(r, pr);
    }

    template <typename I, typename F, typename T, typename... Ps,
              std::enable_if_t<std::is_same<iterate_box_comp,I>::value,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void mf_call_f (F const& f, int ibox, int i, int j, int k, int ncomp, T& r) noexcept
    {
        for (int n = 0; n < ncomp; ++n) {
            auto const& pr = f(ibox,i,j,k,n);
            Reduce::detail::for_each_local<0, T, Ps...>(r, pr);
        }
    }
}

template <typename... Ps>
class ReduceOps
{
public:

    // This is public for CUDA
    template <typename I, typename MF, typename D, typename F>
    void eval_mf (I, MF const& mf, IntVect const& nghost, int ncomp, D& reduce_data, F&&f)
    {
        using ReduceTuple = typename D::Type;
        const int nboxes = mf.local_size();
        if (nboxes > 0) {
            auto const& parforinfo = mf.getParForInfo(nghost,AMREX_GPU_MAX_THREADS);
            auto par_for_blocks = parforinfo.getBlocks();
            const int nblocks = par_for_blocks.first[nboxes];
            const int block_0_size = par_for_blocks.first[1];
            const int* dp_nblocks = par_for_blocks.second;
            const Box* dp_boxes = parforinfo.getBoxes();

            auto const& stream = Gpu::gpuStream();
            auto pdst = reduce_data.devicePtr(stream);
            int nblocks_ec = std::min(nblocks, reduce_data.maxBlocks());
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
                    int ibox, icell;
                    if (dp_nblocks) {
                        ibox = amrex::bisect(dp_nblocks, 0, nboxes, iblock);
                        icell = (iblock-dp_nblocks[ibox])*AMREX_GPU_MAX_THREADS + threadIdx.x;
                    } else {
                        ibox = iblock / block_0_size;
                        icell = (iblock-ibox*block_0_size)*AMREX_GPU_MAX_THREADS + threadIdx.x;
                    }

                    Box const& b = dp_boxes[ibox];
                    int ncells = b.numPts();
                    if (icell < ncells) {
                        const auto len = amrex::length(b);
                        int k =  icell /   (len.x*len.y);
                        int j = (icell - k*(len.x*len.y)) /   len.x;
                        int i = (icell - k*(len.x*len.y)) - j*len.x;
                        AMREX_D_TERM(i += b.smallEnd(0);,
                                     j += b.smallEnd(1);,
                                     k += b.smallEnd(2););
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

    template <typename D, typename F>
    void eval (Box const& box, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        auto const& stream = Gpu::gpuStream();
        auto dp = reduce_data.devicePtr(stream);
        int& nblocks = reduce_data.nBlocks(stream);
        int ncells = box.numPts();
        const auto lo  = amrex::lbound(box);
        const auto len = amrex::length(box);
        const auto lenxy = len.x*len.y;
        const auto lenx = len.x;
        IndexType ixtype = box.ixType();
        constexpr int nitems_per_thread = 4;
        int nblocks_ec = (ncells + nitems_per_thread*AMREX_GPU_MAX_THREADS-1)
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
            Dim1 blockDim {gh.blockDim()};
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
            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < ncells; icell += stride) {
                int k =  icell /   lenxy;
                int j = (icell - k*lenxy) /   lenx;
                int i = (icell - k*lenxy) - j*lenx;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                auto pr = Reduce::detail::call_f(f,i,j,k,ixtype);
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
            }
#ifdef AMREX_USE_SYCL
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r, gh);
#else
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r);
#endif
        });
        nblocks = std::max(nblocks, nblocks_ec);
    }

    template <typename N, typename D, typename F,
              typename M=std::enable_if_t<std::is_integral<N>::value> >
    void eval (Box const& box, N ncomp, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        auto const& stream = Gpu::gpuStream();
        auto dp = reduce_data.devicePtr(stream);
        int& nblocks = reduce_data.nBlocks(stream);
        int ncells = box.numPts();
        const auto lo  = amrex::lbound(box);
        const auto len = amrex::length(box);
        const auto lenxy = len.x*len.y;
        const auto lenx = len.x;
        constexpr int nitems_per_thread = 4;
        int nblocks_ec = (ncells + nitems_per_thread*AMREX_GPU_MAX_THREADS-1)
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
            Dim1 blockDim {gh.blockDim()};
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
            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < ncells; icell += stride) {
                int k =  icell /   lenxy;
                int j = (icell - k*lenxy) /   lenx;
                int i = (icell - k*lenxy) - j*lenx;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                for (N n = 0; n < ncomp; ++n) {
                    auto pr = f(i,j,k,n);
                    Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
                }
            }
#ifdef AMREX_USE_SYCL
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r, gh);
#else
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r);
#endif
        });
        nblocks = std::max(nblocks, nblocks_ec);
    }

    template <typename N, typename D, typename F,
              typename M=std::enable_if_t<std::is_integral<N>::value> >
    void eval (N n, D & reduce_data, F&& f)
    {
        if (n <= 0) return;
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
            Dim1 blockDim {gh.blockDim()};
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
            for (N i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 i < n; i += stride) {
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
        using ReduceTuple = typename D::Type;
        auto const& stream = Gpu::gpuStream();
        auto hp = reduce_data.hostPtr();
        auto dp = reduce_data.devicePtr();
        auto const& nblocks = reduce_data.nBlocks();
#if defined(AMREX_USE_SYCL)
        if (reduce_data.maxStreamIndex() == 0 && nblocks[0] <= 4096) {
            const int N = nblocks[0];
            if (N == 0) {
                Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(*hp);
                return *hp;
            } else {
                Gpu::PinnedVector<ReduceTuple> tmp(N);
                Gpu::dtoh_memcpy_async(tmp.data(), dp, sizeof(ReduceTuple)*N);
                Gpu::streamSynchronize();
                for (int i = 1; i < N; ++i) {
                    Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(tmp[0], tmp[i]);
                }
                return tmp[0];
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
                    for (int i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                         i < nblocks[istream]; i += stride) {
                        Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, dp_stream[i]);
                    }
                }
                Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(dst, r);
                if (threadIdx.x == 0) { *hp = dst; }
            });
#endif
            Gpu::streamSynchronize();
            return *hp;
        }
    }
};

namespace Reduce {

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral<N>::value> >
T Sum (N n, T const* v, T init_val = 0)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {v[i]}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return amrex::get<0>(hv) + init_val;
}

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral<N>::value> >
T Sum (N n, F&& f, T init_val = 0)
{
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {f(i)}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return amrex::get<0>(hv) + init_val;
}

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral<N>::value> >
T Min (N n, T const* v, T init_val = std::numeric_limits<T>::max())
{
    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {v[i]}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return std::min(amrex::get<0>(hv),init_val);
}

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral<N>::value> >
T Min (N n, F&& f, T init_val = std::numeric_limits<T>::max())
{
    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {f(i)}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return std::min(amrex::get<0>(hv),init_val);
}

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral<N>::value> >
T Max (N n, T const* v, T init_val = std::numeric_limits<T>::lowest())
{
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {v[i]}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return std::max(amrex::get<0>(hv),init_val);
}

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral<N>::value> >
T Max (N n, F&& f, T init_val = std::numeric_limits<T>::lowest())
{
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<T> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(n, reduce_data, [=] AMREX_GPU_DEVICE (N i) -> ReduceTuple { return {f(i)}; });
    ReduceTuple hv = reduce_data.value(reduce_op);
    return std::max(amrex::get<0>(hv),init_val);
}

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral<N>::value> >
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

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral<N>::value> >
std::pair<T,T> MinMax (N n, F&& f)
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

template <typename T, typename N, typename P, typename M=std::enable_if_t<std::is_integral<N>::value> >
bool AnyOf (N n, T const* v, P&& pred)
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
            for (N i = gh.blockDim()*gh.blockIdx()+gh.threadIdx(), stride = gh.blockDim()*gh.gridDim();
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
        if (threadIdx.x == 0) has_any = *dp;
        __syncthreads();

        if (!has_any)
        {
            int r = false;
            for (N i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
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

template <typename P>
bool AnyOf (Box const& box, P&& pred)
{
    Gpu::LaunchSafeGuard lsg(true);
    Gpu::DeviceScalar<int> ds(0);
    int* dp = ds.dataPtr();
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto lenxy = len.x*len.y;
    const auto lenx = len.x;
    auto ec = Gpu::ExecutionConfig(ncells);
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
            for (int icell = gh.blockDim()*gh.blockIdx()+gh.threadIdx(), stride = gh.blockDim()*gh.gridDim();
                 icell < ncells && !r; icell += stride) {
                int k =  icell /   lenxy;
                int j = (icell - k*lenxy) /   lenx;
                int i = (icell - k*lenxy) - j*lenx;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                r = pred(i,j,k) ? 1 : 0;
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
        if (threadIdx.x == 0) has_any = *dp;
        __syncthreads();

        if (!has_any)
        {
            int r = false;
            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < ncells && !r; icell += stride) {
                int k =  icell /   lenxy;
                int j = (icell - k*lenxy) /   lenx;
                int i = (icell - k*lenxy) - j*lenx;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                r = pred(i,j,k) ? 1 : 0;
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

template <typename... Ps>
class ReduceOps
{
private:

    template <typename D, typename F>
    AMREX_FORCE_INLINE
    static auto call_f (Box const& box, typename D::Type & r, F const& f)
        noexcept -> std::enable_if_t<std::is_same<std::decay_t<decltype(f(0,0,0))>,
                                                  typename D::Type>::value>
    {
        using ReduceTuple = typename D::Type;
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, f(i,j,k));
        }}}
    }

    template <typename D, typename F>
    AMREX_FORCE_INLINE
    static auto call_f (Box const& box, typename D::Type & r, F const& f)
        noexcept -> std::enable_if_t<std::is_same<std::decay_t<decltype(f(Box()))>,
                                                  typename D::Type>::value>
    {
        using ReduceTuple = typename D::Type;
        Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, f(box));
    }

public:

    template <typename MF, typename D, typename F>
    std::enable_if_t<IsFabArray<MF>::value && IsCallable<F, int, int, int, int>::value>
    eval (MF const& mf, IntVect const& nghost, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
            Box const& b = mfi.growntilebox(nghost);
            const int li = mfi.LocalIndex();
            auto& rr = reduce_data.reference(OpenMP::get_thread_num());
            const auto lo = amrex::lbound(b);
            const auto hi = amrex::ubound(b);
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(rr, f(li,i,j,k));
            }}}
        }
    }

    template <typename MF, typename D, typename F>
    std::enable_if_t<IsFabArray<MF>::value && IsCallable<F, int, int, int, int, int>::value>
    eval (MF const& mf, IntVect const& nghost, int ncomp, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
            Box const& b = mfi.growntilebox(nghost);
            const int li = mfi.LocalIndex();
            auto& rr = reduce_data.reference(OpenMP::get_thread_num());
            const auto lo = amrex::lbound(b);
            const auto hi = amrex::ubound(b);
            for (int n = 0; n < ncomp; ++n) {
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(rr, f(li,i,j,k,n));
            }}}}
        }
    }

    template <typename D, typename F>
    void eval (Box const& box, D & reduce_data, F&& f)
    {
        auto& rr = reduce_data.reference(OpenMP::get_thread_num());
        call_f<D>(box, rr, f);
    }

    template <typename N, typename D, typename F,
              typename M=std::enable_if_t<std::is_integral<N>::value> >
    void eval (Box const& box, N ncomp, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        auto& rr = reduce_data.reference(OpenMP::get_thread_num());
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);
        for (N n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(rr, f(i,j,k,n));
        }}}}
    }

    template <typename N, typename D, typename F,
              typename M=std::enable_if_t<std::is_integral<N>::value> >
    void eval (N n, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        auto& rr = reduce_data.reference(OpenMP::get_thread_num());
        for (N i = 0; i < n; ++i) {
            Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(rr, f(i));
        }
    }

    template <typename D>
    typename D::Type value (D & reduce_data)
    {
        using ReduceTuple = typename D::Type;
        auto& rrv = reduce_data.reference();
        if (rrv.size() > 1) {
            for (int i = 1, N = rrv.size(); i < N; ++i) {
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(rrv[0], rrv[i]);
            }
        }
        return rrv[0];
    }
};

namespace Reduce {

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral<N>::value> >
T Sum (N n, F&& f, T init_val = 0)
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

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral<N>::value> >
T Sum (N n, T const* v, T init_val = 0)
{
    return Sum(n, [=] (N i) -> T { return v[i]; }, init_val);
}

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral<N>::value> >
T Min (N n, F&& f, T init_val = std::numeric_limits<T>::max())
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

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral<N>::value> >
T Min (N n, T const* v, T init_val = std::numeric_limits<T>::max())
{
    return Reduce::Min(n, [=] (N i) -> T { return v[i]; }, init_val);
}

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral<N>::value> >
T Max (N n, F&& f, T init_val = std::numeric_limits<T>::lowest())
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

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral<N>::value> >
T Max (N n, T const* v, T init_val = std::numeric_limits<T>::lowest())
{
    return Reduce::Max(n, [=] (N i) -> T { return v[i]; }, init_val);
}

template <typename T, typename N, typename F,
          typename M=std::enable_if_t<std::is_integral<N>::value> >
std::pair<T,T> Min (N n, F&& f)
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

template <typename T, typename N, typename M=std::enable_if_t<std::is_integral<N>::value> >
std::pair<T,T> MinMax (N n, T const* v)
{
    return Reduce::MinMax<T>(n, [=] (N i) -> T { return v[i]; });
}

template <typename T, typename N, typename P, typename M=std::enable_if_t<std::is_integral<N>::value> >
bool AnyOf (N n, T const* v, P&& pred)
{
    return std::any_of(v, v+n, pred);
}

template <typename P>
bool AnyOf (Box const& box, P&&pred)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        if (pred(i,j,k)) return true;
    }}}
    return false;
}

}

#endif

}

#endif
