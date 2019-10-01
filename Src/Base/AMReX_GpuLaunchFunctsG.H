#ifndef AMREX_GPU_LAUNCH_FUNCTS_G_H_
#define AMREX_GPU_LAUNCH_FUNCTS_G_H_

namespace amrex {

template<typename T, typename L>
void launch (T const& n, L&& f, std::size_t shared_mem_bytes=0) noexcept
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
void For (T n, L&& f, std::size_t shared_mem_bytes=0) noexcept
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
void ParallelFor (T n, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    For(n,std::forward<L>(f),shared_mem_bytes);
}

template <typename L>
void For (Box const& box, L&& f, std::size_t shared_mem_bytes=0) noexcept
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
void ParallelFor (Box const& box, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    For(box,std::forward<L>(f),shared_mem_bytes);
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void For (Box const& box, T ncomp, L&& f, std::size_t shared_mem_bytes=0) noexcept
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
void ParallelFor (Box const& box, T ncomp, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    For(box,ncomp,std::forward<L>(f),shared_mem_bytes);
}

template <typename L1, typename L2>
void For (Box const& box1, Box const& box2, L1&& f1, L2&& f2,
          std::size_t shared_mem_bytes=0) noexcept
{
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells = amrex::max(ncells1, ncells2);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto ec = Gpu::ExecutionConfig(ncells);
    std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < ncells1) {
                int k =  icell /   (len1.x*len1.y);
                int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                i += lo1.x;
                j += lo1.y;
                k += lo1.z;
                f1(i,j,k);
            }
            if (icell < ncells2) {
                int k =  icell /   (len2.x*len2.y);
                int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                i += lo2.x;
                j += lo2.y;
                k += lo2.z;
                f2(i,j,k);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename L1, typename L2, typename L3>
void For (Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3,
          std::size_t shared_mem_bytes=0) noexcept
{
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells3 = box3.numPts();
    int ncells = amrex::max(ncells1, ncells2, ncells3);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto lo3  = amrex::lbound(box3);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto len3 = amrex::length(box3);
    const auto ec = Gpu::ExecutionConfig(ncells);
    std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < ncells1) {
                int k =  icell /   (len1.x*len1.y);
                int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                i += lo1.x;
                j += lo1.y;
                k += lo1.z;
                f1(i,j,k);
            }
            if (icell < ncells2) {
                int k =  icell /   (len2.x*len2.y);
                int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                i += lo2.x;
                j += lo2.y;
                k += lo2.z;
                f2(i,j,k);
            }
            if (icell < ncells3) {
                int k =  icell /   (len3.x*len3.y);
                int j = (icell - k*(len3.x*len3.y)) /   len3.x;
                int i = (icell - k*(len3.x*len3.y)) - j*len3.x;
                i += lo3.x;
                j += lo3.y;
                k += lo3.z;
                f3(i,j,k);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2,
          std::size_t shared_mem_bytes=0) noexcept
{
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells = amrex::max(ncells1, ncells2);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto ec = Gpu::ExecutionConfig(ncells);
    std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < ncells1) {
                int k =  icell /   (len1.x*len1.y);
                int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                i += lo1.x;
                j += lo1.y;
                k += lo1.z;
                for (T1 n = 0; n < ncomp1; ++n) {
                    f1(i,j,k,n);
                }
            }
            if (icell < ncells2) {
                int k =  icell /   (len2.x*len2.y);
                int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                i += lo2.x;
                j += lo2.y;
                k += lo2.z;
                for (T2 n = 0; n < ncomp2; ++n) {
                    f2(i,j,k,n);
                }
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value>,
          typename M3=amrex::EnableIf_t<std::is_integral<T3>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2,
          Box const& box3, T3 ncomp3, L3&& f3,
          std::size_t shared_mem_bytes=0) noexcept
{
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells3 = box3.numPts();
    int ncells = amrex::max(ncells1, ncells2, ncells3);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto lo3  = amrex::lbound(box3);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto len3 = amrex::length(box3);
    const auto ec = Gpu::ExecutionConfig(ncells);
    std::size_t sm = std::max(ec.sharedMem, shared_mem_bytes);
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, sm, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < ncells1) {
                int k =  icell /   (len1.x*len1.y);
                int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                i += lo1.x;
                j += lo1.y;
                k += lo1.z;
                for (T1 n = 0; n < ncomp1; ++n) {
                    f1(i,j,k,n);
                }
            }
            if (icell < ncells2) {
                int k =  icell /   (len2.x*len2.y);
                int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                i += lo2.x;
                j += lo2.y;
                k += lo2.z;
                for (T2 n = 0; n < ncomp2; ++n) {
                    f2(i,j,k,n);
                }
            }
            if (icell < ncells3) {
                int k =  icell /   (len3.x*len3.y);
                int j = (icell - k*(len3.x*len3.y)) /   len3.x;
                int i = (icell - k*(len3.x*len3.y)) - j*len3.x;
                i += lo3.x;
                j += lo3.y;
                k += lo3.z;
                for (T3 n = 0; n < ncomp3; ++n) {
                    f3(i,j,k,n);
                }
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename L1, typename L2>
void ParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2,
                  std::size_t shared_mem_bytes=0) noexcept
{
    For(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2),shared_mem_bytes);
}

template <typename L1, typename L2, typename L3>
void ParallelFor (Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3,
                  std::size_t shared_mem_bytes=0) noexcept
{
    For(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3),shared_mem_bytes);
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2,
                  std::size_t shared_mem_bytes=0) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2),shared_mem_bytes);
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value>,
          typename M3=amrex::EnableIf_t<std::is_integral<T3>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2,
                  Box const& box3, T3 ncomp3, L3&& f3,
                  std::size_t shared_mem_bytes=0) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L2>(f3),shared_mem_bytes);
}

template <typename T, typename L1, typename L2>
void FabReduce (Box const& box, T const& init_val,
                L1&& f1, L2&& f2, std::size_t shared_mem_bytes=0) noexcept
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
                L1&& f1, L2&& f2, std::size_t shared_mem_bytes=0) noexcept
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

}

#endif
