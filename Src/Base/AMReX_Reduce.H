#ifndef AMREX_REDUCE_H_
#define AMREX_REDUCE_H_

#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>

namespace amrex {

namespace Reduce { namespace detail {
    template <std::size_t I, typename T, typename P>
    AMREX_GPU_HOST_DEVICE
    void for_each_parallel (T& d, T const& s)
    {
        P().parallel_update(amrex::get<I>(d), amrex::get<I>(s));
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_GPU_HOST_DEVICE
    void for_each_parallel (T& d, T const& s)
    {
        P().parallel_update(amrex::get<I>(d), amrex::get<I>(s));
        for_each_parallel<I+1,T,P1,Ps...>(d, s);
    }

    template <std::size_t I, typename T, typename P>
    AMREX_GPU_HOST_DEVICE
    void for_each_local (T& d, T const& s)
    {
        P().local_update(amrex::get<I>(d), amrex::get<I>(s));
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_GPU_HOST_DEVICE
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
    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s) const noexcept { Gpu::deviceReduceSum(&d,s); }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void local_update (T& d, T const& s) const noexcept { d += s; }

    template <typename T>
    void init (T& t) const noexcept { t = 0; }
};

struct ReduceOpMin
{
    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s) const noexcept { Gpu::deviceReduceMin(&d,s); }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void local_update (T& d, T const& s) const noexcept { d = amrex::min(d,s); }

    template <typename T>
    void init (T& t) const noexcept { t = std::numeric_limits<T>::max(); }
};

struct ReduceOpMax
{
    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void parallel_update (T& d, T const& s) const noexcept { Gpu::deviceReduceMax(&d,s); }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void local_update (T& d, T const& s) const noexcept { d = amrex::max(d,s); }

    template <typename T>
    void init (T& t) const noexcept { t = std::numeric_limits<T>::lowest(); }
};

template <typename... Ps> class ReduceOps;

#ifdef AMREX_USE_GPU

namespace Reduce {
namespace detail {
    template <typename T>
    void init_tuple_on_device (T* p, T const& v)
    {
        amrex::launch_global<<<1,1,0,Gpu::gpuStream()>>>(
        [=] AMREX_GPU_DEVICE () noexcept
        {
            new (p) T(v);
        });
    }
}}

template <typename... Ts>
class ReduceData
{
public:
    using Type = GpuTuple<Ts...>;

    template <typename... Ps>
    explicit ReduceData (ReduceOps<Ps...> const& ops)
        : m_host_tuple(),
          m_device_tuple((Type*)(The_Device_Arena()->alloc(sizeof(m_host_tuple))))
    {
        static_assert(AMREX_IS_TRIVIALLY_COPYABLE(Type),
                      "ReduceData::Type must be trivially copyable");
        static_assert(std::is_trivially_destructible<Type>::value,
                      "ReduceData::Type must be trivially destructible");

        Reduce::detail::for_each_init<0, Type, Ps...>(m_host_tuple);

        Reduce::detail::init_tuple_on_device(m_device_tuple, m_host_tuple);
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

    Type initialValue () const { return m_host_tuple; }

private:
    Type m_host_tuple;
    Type* m_device_tuple;
};

template <typename... Ps>
class ReduceOps
{
public:

    template <typename D, typename F>
    void eval (Box const& box, D & reduce_data, F&& f, std::size_t shared_mem_bytes = 0)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple init_val = reduce_data.initialValue();
        auto dp = reduce_data.devicePtr();
        int ncells = box.numPts();
        const auto lo  = amrex::lbound(box);
        const auto len = amrex::length(box);
        auto ec = Gpu::ExecutionConfig(ncells);
        ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
        std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
        amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, amrex::Gpu::gpuStream()>>>(
        [=] AMREX_GPU_DEVICE () noexcept {
            ReduceTuple r = init_val;
            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < ncells; icell += stride) {
                int k =  icell /   (len.x*len.y);
                int j = (icell - k*(len.x*len.y)) /   len.x;
                int i = (icell - k*(len.x*len.y)) - j*len.x;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                auto pr = f(i,j,k);
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
            }
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(*dp,r);
        });
        AMREX_GPU_ERROR_CHECK();
    }

    template <typename N, typename D, typename F,
              typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
    void eval (Box const& box, N ncomp, D & reduce_data, F&& f, std::size_t shared_mem_bytes = 0)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple init_val = reduce_data.initialValue();
        auto dp = reduce_data.devicePtr();
        int ncells = box.numPts();
        const auto lo  = amrex::lbound(box);
        const auto len = amrex::length(box);
        auto ec = Gpu::ExecutionConfig(ncells);
        ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
        std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
        amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, amrex::Gpu::gpuStream()>>>(
        [=] AMREX_GPU_DEVICE () noexcept {
            ReduceTuple r = init_val;
            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < ncells; icell += stride) {
                int k =  icell /   (len.x*len.y);
                int j = (icell - k*(len.x*len.y)) /   len.x;
                int i = (icell - k*(len.x*len.y)) - j*len.x;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                for (N n = 0; n < ncomp; ++n) {
                    auto pr = f(i,j,k);
                    Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
                }
            }
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(*dp,r);
        });
        AMREX_GPU_ERROR_CHECK();
    }

    template <typename N, typename D, typename F,
              typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
    void eval (N n, D & reduce_data, F&& f, std::size_t shared_mem_bytes = 0)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple init_val = reduce_data.initialValue();
        auto dp = reduce_data.devicePtr();
        auto ec = Gpu::ExecutionConfig(n);
        ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
        std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
        amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, Gpu::gpuStream()>>>(
        [=] AMREX_GPU_DEVICE () noexcept {
            auto r = init_val;
            for (N i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 i < n; i += stride) {
                auto pr = f(i);
                Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r,pr);
            }
            Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(*dp,r);
        });
        AMREX_GPU_ERROR_CHECK();
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
    [=] AMREX_GPU_DEVICE (T const& r) noexcept
    {
        Gpu::deviceReduceSum(dp, r);
    });
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
    [=] AMREX_GPU_DEVICE (T const& r) noexcept
    {
        Gpu::deviceReduceMin(dp, r);
    });
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
    [=] AMREX_GPU_DEVICE (T const& r) noexcept
    {
        Gpu::deviceReduceMax(dp, r);
    });
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
    [=] AMREX_GPU_DEVICE (Real2 const& r) noexcept
    {
        Gpu::deviceReduceMin(dp  , r[0]);
        Gpu::deviceReduceMax(dp+1, r[1]);
    });
    Gpu::dtoh_memcpy(hv.data(), dp, 2*sizeof(T));
    The_Device_Arena()->free(dp);
    return std::make_pair(hv[0],hv[1]);
}

template <typename T, typename N, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
std::pair<T,T> MinMax (N n, T const* v)
{
    return Reduce::MinMax<T>(n, v, amrex::Less<T>(), amrex::Greater<T>());
}

}

#else

template <typename... Ts>
class ReduceData
{
public:
    using Type = GpuTuple<Ts...>;

    template <typename... Ps>
    explicit ReduceData (ReduceOps<Ps...> const& ops)
        : m_init_val(), m_tuple()
    {
        Reduce::detail::for_each_init<0, Type, Ps...>(m_init_val);
        m_tuple = m_init_val;
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

    Type initialValue () const { return m_init_val; }

private:
    Type m_init_val;
    Type m_tuple;
};

template <typename... Ps>
class ReduceOps
{
public:

    template <typename D, typename F>
    void eval (Box const& box, D & reduce_data, F&& f, std::size_t shared_mem_bytes = 0)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple r = reduce_data.initialValue();
        ReduceTuple& rr = reduce_data.reference();
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            auto pr = f(i,j,k);
            Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
        }}}
        Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(rr,r);
    }

    template <typename N, typename D, typename F,
              typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
    void eval (Box const& box, N ncomp, D & reduce_data, F&& f, std::size_t shared_mem_bytes = 0)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple r = reduce_data.initialValue();
        ReduceTuple& rr = reduce_data.reference();
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);
        for (N n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            auto pr = f(i,j,k);
            Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(r, pr);
        }}}}
        Reduce::detail::for_each_parallel<0, ReduceTuple, Ps...>(rr,r);
    }

    template <typename N, typename D, typename F,
              typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
    void eval (N n, D & reduce_data, F&& f, std::size_t shared_mem_bytes = 0)
    {
        using ReduceTuple = typename D::Type;
        ReduceTuple r = reduce_data.initialValue();
        ReduceTuple& rr = reduce_data.reference();
        AMREX_PRAGMA_SIMD
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
    [=] AMREX_GPU_DEVICE (N i, T* r) noexcept
    {
        *r = bop(*r, v[i]);
    },
    [=] AMREX_GPU_DEVICE (T r) noexcept
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
    [=] AMREX_GPU_DEVICE (N i, T* r) noexcept
    {
        *r = bop(*r, v[i]);
    },
    [=] AMREX_GPU_DEVICE (T r) noexcept
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
    [=] AMREX_GPU_DEVICE (N i, Real2* r) noexcept
    {
        (*r)[0] = minop((*r)[0], v[i]);
        (*r)[1] = maxop((*r)[1], v[i]);
    },
    [=] AMREX_GPU_DEVICE (Real2 const& r) noexcept
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

}

#endif

}

#endif
