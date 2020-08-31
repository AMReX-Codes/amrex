#ifndef AMREX_REDUCE_H_
#define AMREX_REDUCE_H_

#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>

#include <algorithm>

namespace amrex {

namespace Reduce { namespace detail {

#ifdef AMREX_USE_DPCPP
    template <std::size_t I, typename T, typename P>
    AMREX_GPU_DEVICE
    void for_each_parallel (T& d, T const& s, Gpu::Handler const& h)
    {
        P().parallel_update(amrex::get<I>(d), amrex::get<I>(s), h);
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_GPU_DEVICE
    void for_each_parallel (T& d, T const& s, Gpu::Handler const& h)
    {
        P().parallel_update(amrex::get<I>(d), amrex::get<I>(s), h);
        for_each_parallel<I+1,T,P1,Ps...>(d, s, h);
    }
#else
    template <std::size_t I, typename T, typename P>
    AMREX_GPU_DEVICE
    void for_each_parallel (T& d, T const& s)
    {
        P().parallel_update(amrex::get<I>(d), amrex::get<I>(s));
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_GPU_DEVICE
    void for_each_parallel (T& d, T const& s)
    {
        P().parallel_update(amrex::get<I>(d), amrex::get<I>(s));
        for_each_parallel<I+1,T,P1,Ps...>(d, s);
    }
#endif
    template <std::size_t I, typename T, typename P>
    AMREX_GPU_DEVICE
    void for_each_local (T& d, T const& s)
    {
        P().local_update(amrex::get<I>(d), amrex::get<I>(s));
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_GPU_DEVICE
    void for_each_local (T& d, T const& s)
    {
        P().local_update(amrex::get<I>(d), amrex::get<I>(s));
        for_each_local<I+1,T,P1,Ps...>(d, s);
    }

    template <std::size_t I, typename T, typename P>
    void for_each_init (T& t)
    {
        P().init(amrex::get<I>(t));
    }
    
    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    void for_each_init (T& t)
    {
        P().init(amrex::get<I>(t));
        for_each_init<I+1,T,P1,Ps...>(t);
    }
}}

struct ReduceOpSum
{
#ifdef AMREX_USE_DPCPP
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s, Gpu::Handler const& h) const noexcept {
        Gpu::deviceReduceSum(&d,s,h);
    }
#else
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s) const noexcept { Gpu::deviceReduceSum(&d,s); }
#endif

    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void local_update (T& d, T const& s) const noexcept { d += s; }

    template <typename T>
    void init (T& t) const noexcept { t = 0; }
};

struct ReduceOpMin
{
#ifdef AMREX_USE_DPCPP
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s, Gpu::Handler const& h) const noexcept {
        Gpu::deviceReduceMin(&d,s,h);
    }
#else
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s) const noexcept { Gpu::deviceReduceMin(&d,s); }
#endif

    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void local_update (T& d, T const& s) const noexcept { d = amrex::min(d,s); }

    template <typename T>
    void init (T& t) const noexcept { t = std::numeric_limits<T>::max(); }
};

struct ReduceOpMax
{
#ifdef AMREX_USE_DPCPP
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s, Gpu::Handler const& h) const noexcept {
        Gpu::deviceReduceMax(&d,s,h);
    }
#else
    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s) const noexcept { Gpu::deviceReduceMax(&d,s); }
#endif

    template <typename T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void local_update (T& d, T const& s) const noexcept { d = amrex::max(d,s); }

    template <typename T>
    void init (T& t) const noexcept { t = std::numeric_limits<T>::lowest(); }
};

struct ReduceOpLogicalAnd
{
#ifdef AMREX_USE_DPCPP
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (int& d, int s, Gpu::Handler const& h) const noexcept {
        Gpu::deviceReduceLogicalAnd(&d,s,h);
    }
#else
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (int& d, int s) const noexcept { Gpu::deviceReduceLogicalAnd(&d,s); }
#endif

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void local_update (int& d, int s) const noexcept { d = d && s; }

    void init (int& t) const noexcept { t = true; }
};

struct ReduceOpLogicalOr
{
#ifdef AMREX_USE_DPCPP
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (int& d, int s, Gpu::Handler const& h) const noexcept {
        Gpu::deviceReduceLogicalOr(&d,s,h);
    }
#else
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void parallel_update (int& d, int s) const noexcept { Gpu::deviceReduceLogicalOr(&d,s); }
#endif

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void local_update (int& d, int s) const noexcept { d = d || s; }

    void init (int& t) const noexcept { t = false; }
};

template <typename... Ps> class ReduceOps;

#ifdef AMREX_USE_GPU

namespace Reduce {
namespace detail {
#ifdef AMREX_USE_DPCPP
    // xxxxx DPCPP todo: the compiler seems to have a bug.  It cannot capture
    //                   GpuTuple<T> correctly even if it's standard layout.
    template <typename T>
    void
    init_tuple_on_device (T* p, T const& v)
    {
        amrex::single_task(Gpu::nullStream(), [=] AMREX_GPU_DEVICE () noexcept
        {
            new (p) T();
            new (p+1) T();
        });
        Array<T,2> tmp{v,v};
        Gpu::htod_memcpy(p, tmp.data(), sizeof(T)*2);
        Gpu::synchronize();
    }
#else
    template <typename T>
    void
    init_tuple_on_device (T* p, T const& v)
    {
        amrex::single_task(Gpu::nullStream(), [=] AMREX_GPU_DEVICE () noexcept
        {
            new (p) T(v);
            new (p+1) T(v);
        });
    }
#endif
}}

template <typename... Ts>
class ReduceData
{
public:
    using Type = GpuTuple<Ts...>;

    template <typename... Ps>
    explicit ReduceData (ReduceOps<Ps...> const&)
        : m_host_tuple(),
          m_device_tuple((Type*)(The_Device_Arena()->alloc(2*sizeof(m_host_tuple))))
    {
        static_assert(AMREX_IS_TRIVIALLY_COPYABLE(Type),
                      "ReduceData::Type must be trivially copyable");
        static_assert(std::is_trivially_destructible<Type>::value,
                      "ReduceData::Type must be trivially destructible");

        Reduce::detail::for_each_init<0, Type, Ps...>(m_host_tuple);
        Reduce::detail::init_tuple_on_device(m_device_tuple, m_host_tuple);
        Gpu::synchronize();
    }

    ~ReduceData () { The_Device_Arena()->free(m_device_tuple); }

    ReduceData (ReduceData<Ts...> const&) = delete;
    ReduceData (ReduceData<Ts...> &&) = delete;
    void operator= (ReduceData<Ts...> const&) = delete;
    void operator= (ReduceData<Ts...> &&) = delete;

    Type value ()
    {
        Gpu::dtoh_memcpy(&m_host_tuple, m_device_tuple, sizeof(m_host_tuple));
        return m_host_tuple;
    }

    Type* devicePtr () { return m_device_tuple; }

    Type& hostRef () { return m_host_tuple; }

private:
    Type m_host_tuple;
    Type* m_device_tuple;
};

namespace Reduce { namespace detail {
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
}}

template <typename... Ps>
class ReduceOps
{
public:

    template <typename D, typename F>
    void eval (Box const& box, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        auto dp = reduce_data.devicePtr();
        int ncells = box.numPts();
        const auto lo  = amrex::lbound(box);
        const auto len = amrex::length(box);
        IndexType ixtype = box.ixType();
        auto ec = Gpu::ExecutionConfig(ncells);
        ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
#ifdef AMREX_USE_DPCPP
        // device reduce needs local(i.e., shared) memory
        constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
        amrex::launch(ec.numBlocks.x, ec.numThreads.x, shared_mem_bytes, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
        {
            ReduceTuple r = *(dp+1);
            for (int icell = gh.item.get_global_id(0), stride = gh.item.get_global_range(0);
                 icell < ncells; icell += stride) {
                int k =  icell /   (len.x*len.y);
                int j = (icell - k*(len.x*len.y)) /   len.x;
                int i = (icell - k*(len.x*len.y)) - j*len.x;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                auto pr = Reduce::detail::call_f(f,i,j,k,ixtype);
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
            }
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(*dp,r, gh);
        });
#else
        amrex::launch(ec.numBlocks.x, ec.numThreads.x, 0, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept {
            ReduceTuple r = *(dp+1);
            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < ncells; icell += stride) {
                int k =  icell /   (len.x*len.y);
                int j = (icell - k*(len.x*len.y)) /   len.x;
                int i = (icell - k*(len.x*len.y)) - j*len.x;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                auto pr = Reduce::detail::call_f(f,i,j,k,ixtype);
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
            }
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(*dp,r);
        });
#endif
    }

    template <typename N, typename D, typename F,
              typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
    void eval (Box const& box, N ncomp, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        auto dp = reduce_data.devicePtr();
        int ncells = box.numPts();
        const auto lo  = amrex::lbound(box);
        const auto len = amrex::length(box);
        auto ec = Gpu::ExecutionConfig(ncells);
        ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
#ifdef AMREX_USE_DPCPP
        // device reduce needs local(i.e., shared) memory
        constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
        amrex::launch(ec.numBlocks.x, ec.numThreads.x, shared_mem_bytes, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
        {
            ReduceTuple r = *(dp+1);
            for (int icell = gh.item.get_global_id(0), stride = gh.item.get_global_range(0);
                 icell < ncells; icell += stride) {
                int k =  icell /   (len.x*len.y);
                int j = (icell - k*(len.x*len.y)) /   len.x;
                int i = (icell - k*(len.x*len.y)) - j*len.x;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                for (N n = 0; n < ncomp; ++n) {
                    auto pr = f(i,j,k,n);
                    Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
                }
            }
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(*dp,r, gh);
        });
#else
        amrex::launch(ec.numBlocks.x, ec.numThreads.x, 0, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept {
            ReduceTuple r = *(dp+1);
            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < ncells; icell += stride) {
                int k =  icell /   (len.x*len.y);
                int j = (icell - k*(len.x*len.y)) /   len.x;
                int i = (icell - k*(len.x*len.y)) - j*len.x;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                for (N n = 0; n < ncomp; ++n) {
                    auto pr = f(i,j,k,n);
                    Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
                }
            }
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(*dp,r);
        });
#endif
    }

    template <typename N, typename D, typename F,
              typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
    void eval (N n, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        auto dp = reduce_data.devicePtr();
        auto ec = Gpu::ExecutionConfig(n);
        ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
#ifdef AMREX_USE_DPCPP
        // device reduce needs local(i.e., shared) memory
        constexpr std::size_t shared_mem_bytes = sizeof(unsigned long long)*Gpu::Device::warp_size;
        amrex::launch(ec.numBlocks.x, ec.numThreads.x, shared_mem_bytes, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& gh) noexcept
        {
            ReduceTuple r = *(dp+1);
            for (N i = gh.item.get_global_id(0), stride = gh.item.get_global_range(0);
                 i < n; i += stride) {
                auto pr = f(i);
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r,pr);
            }
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(*dp,r, gh);
        });
#else
        amrex::launch(ec.numBlocks.x, ec.numThreads.x, 0, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept {
            auto r = *(dp+1);
            for (N i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 i < n; i += stride) {
                auto pr = f(i);
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r,pr);
            }
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(*dp,r);
        });
#endif
    }
};

namespace Reduce {

template <typename T, typename N, typename U, typename BOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Sum (N n, U const* v, T init_val, BOP bop)
{
    Gpu::LaunchSafeGuard lsg(true);
    Gpu::DeviceScalar<T> ds(init_val);
    T* dp = ds.dataPtr();
    amrex::VecReduce(n, init_val,
    [=] AMREX_GPU_DEVICE (N i, T* r) noexcept
    {
        *r = bop(*r, v[i]);
    },
#ifdef AMREX_USE_DPCPP
    [=] AMREX_GPU_DEVICE (T const& r, Gpu::Handler const& h) noexcept
    {
        Gpu::deviceReduceSum(dp, r, h);
    });
#else
    [=] AMREX_GPU_DEVICE (T const& r) noexcept
    {
        Gpu::deviceReduceSum(dp, r);
    });
#endif
    return ds.dataValue();
}

template <typename T, typename N, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Sum (N n, T const* v, T init_val = 0)
{
    return Reduce::Sum(n, v, init_val, amrex::Plus<T>());
}

template <typename T, typename N, typename U, typename BOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Min (N n, U const* v, T init_val, BOP bop)
{
    Gpu::LaunchSafeGuard lsg(true);
    Gpu::DeviceScalar<T> ds(init_val);
    T* dp = ds.dataPtr();
    amrex::VecReduce(n, init_val,
    [=] AMREX_GPU_DEVICE (N i, T* r) noexcept
    {
        *r = bop(*r, v[i]);
    },
#ifdef AMREX_USE_DPCPP
    [=] AMREX_GPU_DEVICE (T const& r, Gpu::Handler const& h) noexcept
    {
        Gpu::deviceReduceMin(dp, r, h);
    });
#else
    [=] AMREX_GPU_DEVICE (T const& r) noexcept
    {
        Gpu::deviceReduceMin(dp, r);
    });
#endif
    return ds.dataValue();
}

template <typename T, typename N, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Min (N n, T const* v, T init_val = std::numeric_limits<T>::max())
{
    return Reduce::Min(n, v, init_val, amrex::Less<T>());
}

template <typename T, typename N, typename U, typename BOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Max (N n, U const* v, T init_val, BOP bop)
{
    Gpu::LaunchSafeGuard lsg(true);
    Gpu::DeviceScalar<T> ds(init_val);
    T* dp = ds.dataPtr();
    amrex::VecReduce(n, init_val,
    [=] AMREX_GPU_DEVICE (N i, T* r) noexcept
    {
        *r = bop(*r, v[i]);
    },
#ifdef AMREX_USE_DPCPP
    [=] AMREX_GPU_DEVICE (T const& r, Gpu::Handler const& h) noexcept
    {
        Gpu::deviceReduceMax(dp, r, h);
    });
#else
    [=] AMREX_GPU_DEVICE (T const& r) noexcept
    {
        Gpu::deviceReduceMax(dp, r);
    });
#endif
    return ds.dataValue();
}

template <typename T, typename N, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Max (N n, T const* v, T init_val = std::numeric_limits<T>::lowest())
{
    return Reduce::Max(n, v, init_val, amrex::Greater<T>());
}

template <typename T, typename N, typename U, typename MINOP, typename MAXOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
std::pair<T,T> MinMax (N n, U const* v, MINOP minop, MAXOP maxop)
{
    Gpu::LaunchSafeGuard lsg(true);
    Array<T,2> hv{std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest()};
    T* dp = (T*)(The_Device_Arena()->alloc(2*sizeof(T)));
    Gpu::htod_memcpy(dp, hv.data(), 2*sizeof(T));
    typedef GpuArray<T,2> Real2;
    amrex::VecReduce(n, Real2{hv[0],hv[1]},
    [=] AMREX_GPU_DEVICE (N i, Real2* r) noexcept
    {
        (*r)[0] = minop((*r)[0], v[i]);
        (*r)[1] = maxop((*r)[1], v[i]);
    },
#ifdef AMREX_USE_DPCPP
    [=] AMREX_GPU_DEVICE (Real2 const& r, Gpu::Handler const& h) noexcept
    {
        Gpu::deviceReduceMin(dp  , r[0], h);
        Gpu::deviceReduceMax(dp+1, r[1], h);
    });
#else
    [=] AMREX_GPU_DEVICE (Real2 const& r) noexcept
    {
        Gpu::deviceReduceMin(dp  , r[0]);
        Gpu::deviceReduceMax(dp+1, r[1]);
    });
#endif
    Gpu::dtoh_memcpy(hv.data(), dp, 2*sizeof(T));
    The_Device_Arena()->free(dp);
    return std::make_pair(hv[0],hv[1]);
}

template <typename T, typename N, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
std::pair<T,T> MinMax (N n, T const* v)
{
    return Reduce::MinMax<T>(n, v, amrex::Less<T>(), amrex::Greater<T>());
}

template <typename T, typename N, typename P, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
bool AnyOf (N n, T const* v, P&& pred)
{
    Gpu::LaunchSafeGuard lsg(true);
    Gpu::DeviceScalar<int> ds(0);
    int* dp = ds.dataPtr();
#ifdef AMREX_USE_DPCPP
    // xxxxx DPCPP todo: Anyof N: better version
    amrex::ParallelFor(n, [=] (int i) noexcept
    {
        int r = pred(v[i]) ? 1 : 0;
        Gpu::Atomic::LogicalOr(dp, r);
    });
#else
    auto ec = Gpu::ExecutionConfig(n);
    ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
    amrex::launch(ec.numBlocks.x, ec.numThreads.x, 0, 0,
    [=] AMREX_GPU_DEVICE () noexcept {
        __shared__ int has_any;
        if (threadIdx.x == 0) has_any = *dp;
        __syncthreads();

        if (!has_any)
        {
            int r = false;
            for (N i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 i < n and !r; i += stride)
            {
                r = pred(v[i]) ? 1 : 0;
            }
            r = Gpu::blockReduce<Gpu::Device::warp_size>
                (r, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::Plus<int> >(), 0);
            if (threadIdx.x == 0 and r) *dp = 1;
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
#ifdef AMREX_USE_DPCPP
    // xxxxx DPCPP todo: Anyof Box: better version
    amrex::ParallelFor(box, [=] (int i, int j, int k) noexcept
    {
        int r = pred(i,j,k) ? 1 : 0;
        Gpu::Atomic::LogicalOr(dp, r);
    });
#else
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, 0,
    [=] AMREX_GPU_DEVICE () noexcept {
        __shared__ int has_any;
        if (threadIdx.x == 0) has_any = *dp;
        __syncthreads();

        if (!has_any)
        {
            int r = false;
            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < ncells and !r; icell += stride) {
                int k =  icell /   (len.x*len.y);
                int j = (icell - k*(len.x*len.y)) /   len.x;
                int i = (icell - k*(len.x*len.y)) - j*len.x;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                r = pred(i,j,k) ? 1 : 0;
            }
            r = Gpu::blockReduce<Gpu::Device::warp_size>
                (r, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::Plus<int> >(), 0);
            if (threadIdx.x == 0 and r) *dp = 1;
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
    explicit ReduceData (ReduceOps<Ps...> const&)
        : m_tuple()
    {
        Reduce::detail::for_each_init<0, Type, Ps...>(m_tuple);
    }

    ReduceData (ReduceData<Ts...> const&) = delete;
    ReduceData (ReduceData<Ts...> &&) = delete;
    void operator= (ReduceData<Ts...> const&) = delete;
    void operator= (ReduceData<Ts...> &&) = delete;

    Type value () const
    {
        return m_tuple;
    }

    Type& reference () { return m_tuple; }

private:
    Type m_tuple;
};

template <typename... Ps>
class ReduceOps
{
private:

    template <typename D, typename F>
    AMREX_FORCE_INLINE
    static auto call_f (Box const& box, D & /*reduce_data*/, F const& f)
        noexcept -> decltype(f(0,0,0))
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple r;
        Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(r);
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            auto pr = f(i,j,k);
            Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
        }}}
        return r;
    }

    template <typename D, typename F>
    AMREX_FORCE_INLINE
    static auto call_f (Box const& box, D&, F const& f)
        noexcept -> decltype(f(Box()))
    {
        return f(box);
    }

public:

    template <typename D, typename F>
    void eval (Box const& box, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple& rr = reduce_data.reference();
        auto r = call_f(box, reduce_data, f);
        Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(rr,r);
    }

    template <typename N, typename D, typename F,
              typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
    void eval (Box const& box, N ncomp, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple r;
        Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(r);
        ReduceTuple& rr = reduce_data.reference();
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);
        for (N n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            auto pr = f(i,j,k,n);
            Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
        }}}}
        Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(rr,r);
    }

    template <typename N, typename D, typename F,
              typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
    void eval (N n, D & reduce_data, F&& f)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple r;
        Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(r);
        ReduceTuple& rr = reduce_data.reference();
        for (N i = 0; i < n; ++i) {
            auto pr = f(i);
            Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
        }
        Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(rr,r);
    }
};

namespace Reduce {

template <typename N, typename T, typename U, typename BOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Sum (N n, U const* v, T init_val, BOP bop)
{
    T sum = init_val;
    T* dp = &sum;
    amrex::VecReduce(n, init_val,
    [=] (N i, T* r) noexcept
    {
        *r = bop(*r, v[i]);
    },
    [=] (T r) noexcept
    {
        Gpu::deviceReduceSum(dp, r);
    });
    return sum;
}

template <typename N, typename T, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Sum (N n, T const* v, T init_val = 0)
{
    return Reduce::Sum(n, v, init_val, amrex::Plus<T>());
}

template <typename N, typename T, typename U, typename BOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Min (N n, U const* v, T init_val, BOP bop)
{
    T mn = init_val;
    T* dp = &mn;
    amrex::VecReduce(n, init_val,
    [=] (N i, T* r) noexcept
    {
        *r = bop(*r, v[i]);
    },
    [=] (T r) noexcept
    {
        Gpu::deviceReduceMin(dp, r);
    });
    return mn;
}

template <typename N, typename T, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Min (N n, T const* v, T init_val = std::numeric_limits<T>::max())
{
    return Reduce::Min(n, v, init_val, amrex::Less<T>());
}

template <typename N, typename T, typename U, typename BOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Max (N n, U const* v, T init_val, BOP bop)
{
    T mx = init_val;
    T* dp = &mx;
    amrex::VecReduce(n, init_val,
    [=] (N i, T* r) noexcept
    {
        *r = bop(*r, v[i]);
    },
    [=] (T r) noexcept
    {
        Gpu::deviceReduceMax(dp, r);
    });
    return mx;
}

template <typename N, typename T, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Max (N n, T const* v, T init_val = std::numeric_limits<T>::lowest())
{
    return Reduce::Max(n, v, init_val, amrex::Greater<T>());
}

template <typename T, typename N, typename U, typename MINOP, typename MAXOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
std::pair<T,T> MinMax (N n, U const* v, MINOP minop, MAXOP maxop)
{
    Array<T,2> hv{std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest()};
    T* dp = hv.data();
    typedef GpuArray<T,2> Real2;
    amrex::VecReduce(n, Real2{hv[0],hv[1]},
    [=] (N i, Real2* r) noexcept
    {
        (*r)[0] = minop((*r)[0], v[i]);
        (*r)[1] = maxop((*r)[1], v[i]);
    },
    [=] (Real2 const& r) noexcept
    {
        Gpu::deviceReduceMin(dp  , r[0]);
        Gpu::deviceReduceMax(dp+1, r[1]);
    });
    return std::make_pair(hv[0],hv[1]);
}

template <typename T, typename N, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
std::pair<T,T> MinMax (N n, T const* v)
{
    return Reduce::MinMax<T>(n, v, amrex::Less<T>(), amrex::Greater<T>());
}

template <typename T, typename N, typename P, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
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
