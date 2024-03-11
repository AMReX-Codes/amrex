#ifndef AMREX_TAG_PARALLELFOR_H_
#define AMREX_TAG_PARALLELFOR_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Vector.H>
#include <limits>
#include <utility>

namespace amrex {

template <class T>
struct Array4PairTag {
    Array4<T      > dfab;
    Array4<T const> sfab;
    Box dbox;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return dbox; }
};

template <class T0, class T1=T0>
struct Array4CopyTag {
    Array4<T0      > dfab;
    Array4<T1 const> sfab;
    Box dbox;
    Dim3 offset; // sbox.smallEnd() - dbox.smallEnd()

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return dbox; }
};

template <class T0, class T1=T0>
struct Array4MaskCopyTag {
    Array4<T0      > dfab;
    Array4<T1 const> sfab;
    Array4<int     > mask;
    Box dbox;
    Dim3 offset; // sbox.smallEnd() - dbox.smallEnd()

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return dbox; }
};

template <class T>
struct Array4Tag {
    Array4<T> dfab;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box box () const noexcept { return Box(dfab); }
};

template <class T>
struct Array4BoxTag {
    Array4<T> dfab;
    Box       dbox;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return dbox; }
};

template <class T>
struct Array4BoxValTag {
    Array4<T> dfab;
    Box       dbox;
    T          val;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return dbox; }
};

template <class T>
struct VectorTag {
    T* p;
    Long m_size;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Long size () const noexcept { return m_size; }
};

#ifdef AMREX_USE_GPU

namespace detail {

template <typename T>
std::enable_if_t<std::is_same<std::decay_t<decltype(std::declval<T>().box())>, Box>::value,
                 Long>
get_tag_size (T const& tag) noexcept
{
    AMREX_ASSERT(tag.box().numPts() < Long(std::numeric_limits<int>::max()));
    return static_cast<int>(tag.box().numPts());
}

template <typename T>
std::enable_if_t<std::is_integral<std::decay_t<decltype(std::declval<T>().size())> >::value,
                 Long>
get_tag_size (T const& tag) noexcept
{
    AMREX_ASSERT(tag.size() < Long(std::numeric_limits<int>::max()));
    return tag.size();
}

template <typename T, typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::enable_if_t<std::is_same<std::decay_t<decltype(std::declval<T>().box())>, Box>::value>
tagparfor_call_f (
#ifdef AMREX_USE_SYCL
    sycl::nd_item<1> const& item,
#endif
    int icell, T const& tag, F&& f) noexcept
{
    int ncells = tag.box().numPts();
    const auto len = amrex::length(tag.box());
    const auto lo  = amrex::lbound(tag.box());
    int k =  icell /   (len.x*len.y);
    int j = (icell - k*(len.x*len.y)) /   len.x;
    int i = (icell - k*(len.x*len.y)) - j*len.x;
    i += lo.x;
    j += lo.y;
    k += lo.z;
#ifdef AMREX_USE_SYCL
    f(item, icell, ncells, i, j, k, tag);
#else
    f(      icell, ncells, i, j, k, tag);
#endif
}

template <typename T, typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::enable_if_t<std::is_integral<std::decay_t<decltype(std::declval<T>().size())> >::value>
tagparfor_call_f (
#ifdef AMREX_USE_SYCL
    sycl::nd_item<1> const& item,
#endif
    int i, T const& tag, F&& f) noexcept
{
    int N = tag.size();
#ifdef AMREX_USE_SYCL
    f(item, i, N, tag);
#else
    f(      i, N, tag);
#endif
}

template <class TagType, class F>
void
ParallelFor_doit (Vector<TagType> const& tags, F && f)
{
    const int ntags = tags.size();
    if (ntags == 0) { return; }

    Long l_ntotwarps = 0;
    int ntotwarps = 0;
    Vector<int> nwarps;
    nwarps.reserve(ntags+1);
    for (int i = 0; i < ntags; ++i)
    {
        auto& tag = tags[i];
        nwarps.push_back(ntotwarps);
        auto nw = (get_tag_size(tag) + Gpu::Device::warp_size-1) / Gpu::Device::warp_size;
        l_ntotwarps += nw;
        ntotwarps += static_cast<int>(nw);
    }
    nwarps.push_back(ntotwarps);

    std::size_t sizeof_tags = ntags*sizeof(TagType);
    std::size_t offset_nwarps = Arena::align(sizeof_tags);
    std::size_t sizeof_nwarps = (ntags+1)*sizeof(int);
    std::size_t total_buf_size = offset_nwarps + sizeof_nwarps;

    char* h_buffer = (char*)The_Pinned_Arena()->alloc(total_buf_size);
    char* d_buffer = (char*)The_Arena()->alloc(total_buf_size);

    std::memcpy(h_buffer, tags.data(), sizeof_tags);
    std::memcpy(h_buffer+offset_nwarps, nwarps.data(), sizeof_nwarps);
    Gpu::htod_memcpy_async(d_buffer, h_buffer, total_buf_size);

    auto d_tags = reinterpret_cast<TagType*>(d_buffer);
    auto d_nwarps = reinterpret_cast<int*>(d_buffer+offset_nwarps);

    constexpr int nthreads = 256;
    constexpr int nwarps_per_block = nthreads/Gpu::Device::warp_size;
    int nblocks = (ntotwarps + nwarps_per_block-1) / nwarps_per_block;

    amrex::ignore_unused(l_ntotwarps);
    AMREX_ASSERT(l_ntotwarps+nwarps_per_block-1 < Long(std::numeric_limits<int>::max()));

    amrex::launch(nblocks, nthreads, Gpu::gpuStream(),
#ifdef AMREX_USE_SYCL
    [=] AMREX_GPU_DEVICE (sycl::nd_item<1> const& item) noexcept
    [[sycl::reqd_work_group_size(1,1,nthreads)]]
    [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
#else
    [=] AMREX_GPU_DEVICE () noexcept
#endif
    {
#ifdef AMREX_USE_SYCL
        std::size_t g_tid = item.get_global_id(0);
#else
        auto g_tid = std::size_t(blockDim.x)*blockIdx.x + threadIdx.x;
#endif
        auto g_wid = int(g_tid / Gpu::Device::warp_size);
        if (g_wid >= ntotwarps) { return; }

        int tag_id = amrex::bisect(d_nwarps, 0, ntags, g_wid);

        int b_wid = g_wid - d_nwarps[tag_id]; // b_wid'th warp on this box
#ifdef AMREX_USE_SYCL
        int lane = item.get_local_id(0) % Gpu::Device::warp_size;
#else
        int lane = threadIdx.x % Gpu::Device::warp_size;
#endif
        int icell = b_wid*Gpu::Device::warp_size + lane;

#ifdef AMREX_USE_SYCL
        tagparfor_call_f(item, icell, d_tags[tag_id], f);
#else
        tagparfor_call_f(      icell, d_tags[tag_id], f);
#endif
    });

    Gpu::streamSynchronize();
    The_Pinned_Arena()->free(h_buffer);
    The_Arena()->free(d_buffer);
}

}

template <class TagType, class F>
std::enable_if_t<std::is_same<std::decay_t<decltype(std::declval<TagType>().box())>,
                              Box>::value>
ParallelFor (Vector<TagType> const& tags, int ncomp, F && f)
{
    detail::ParallelFor_doit(tags,
        [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_SYCL
            sycl::nd_item<1> const& /*item*/,
#endif
            int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
        {
            if (icell < ncells) {
                for (int n = 0; n < ncomp; ++n) {
                    f(i,j,k,n,tag);
                }
            }
        });
}

template <class TagType, class F>
std::enable_if_t<std::is_same<std::decay_t<decltype(std::declval<TagType>().box())>, Box>::value>
ParallelFor (Vector<TagType> const& tags, F && f)
{
    detail::ParallelFor_doit(tags,
        [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_SYCL
            sycl::nd_item<1> const& /*item*/,
#endif
            int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
        {
            if (icell < ncells) {
                f(i,j,k,tag);
            }
        });
}

template <class TagType, class F>
std::enable_if_t<std::is_integral<std::decay_t<decltype(std::declval<TagType>().size())> >::value>
ParallelFor (Vector<TagType> const& tags, F && f)
{
    detail::ParallelFor_doit(tags,
        [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_SYCL
            sycl::nd_item<1> const& /*item*/,
#endif
            int icell, int ncells, TagType const& tag) noexcept
        {
            if (icell < ncells) {
                f(icell,tag);
            }
        });
}

#endif

}

#endif
