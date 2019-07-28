#ifndef AMREX_GPU_LAUNCH_FUNCTS_G_H_
#define AMREX_GPU_LAUNCH_FUNCTS_G_H_

namespace amrex {

template<typename T, typename L>
void launch (T const& n, L f, std::size_t shared_mem_bytes=0) noexcept
{
    const auto ec = Gpu::ExecutionConfig(n);
    std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (auto const i : Gpu::Range(n)) {
            f(i);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void For (T n, L f, std::size_t shared_mem_bytes=0) noexcept
{
    const auto ec = Gpu::ExecutionConfig(n);
    std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (T i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             i < n; i += stride) {
            f(i);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void ParallelFor (T n, L f, std::size_t shared_mem_bytes=0) noexcept
{
    For(n,f,shared_mem_bytes);
}

template <typename L>
void For (Box const& box, L f, std::size_t shared_mem_bytes=0) noexcept
{
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Gpu::ExecutionConfig(ncells);
    std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            f(i,j,k);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename L>
void ParallelFor (Box const& box, L f, std::size_t shared_mem_bytes=0) noexcept
{
    For(box,f,shared_mem_bytes);
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void For (Box const& box, T ncomp, L f, std::size_t shared_mem_bytes=0) noexcept
{
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Gpu::ExecutionConfig(ncells);
    std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            for (T n = 0; n < ncomp; ++n) {
                f(i,j,k,n);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void ParallelFor (Box const& box, T ncomp, L f, std::size_t shared_mem_bytes=0) noexcept
{
    For(box,ncomp,f,shared_mem_bytes);
}

template <typename T, typename L1, typename L2>
void FabReduce (Box const& box, T const& init_val,
                L1 f1, L2 f2, std::size_t shared_mem_bytes=0) noexcept
{
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
    std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        auto r = init_val;
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            f1(i,j,k,&r);
        }
        f2(r);
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void FabReduce (Box const& box, N ncomp, T const& init_val,
                L1 f1, L2 f2, std::size_t shared_mem_bytes=0) noexcept
{
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
    std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        auto r = init_val;
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            for (N n = 0; n < ncomp; ++n) {
                f1(i,j,k,n,&r);
            }
        }
        f2(r);
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void VecReduce (N n, T const& init_val,
                L1&& f1, L2&& f2, std::size_t shared_mem_bytes=0) noexcept
{
    auto ec = Gpu::ExecutionConfig(n);
    ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
    std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        auto r = init_val;
        for (N i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             i < n; i += stride) {
            f1(i,&r);
        }
        f2(r);
    });
    AMREX_GPU_ERROR_CHECK();

}

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
        Gpu::ReduceSum(dp, r);
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
        Gpu::ReduceMin(dp, r);
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
        Gpu::ReduceMax(dp, r);
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
        Gpu::ReduceMin(dp  , r[0]);
        Gpu::ReduceMax(dp+1, r[1]);
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

}

#endif
