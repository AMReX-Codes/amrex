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
    int              dindex;
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
struct Array4BoxOrientationTag {
    Array4<T> fab;
    Box bx;
    Orientation face;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box() const noexcept { return bx; }
};

template <class T>
struct Array4BoxOffsetTag {
    Array4<T> fab;
    Box bx;
    Dim3 offset;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box() const noexcept { return bx; }
};

template <class T>
struct VectorTag {
    T* p;
    Long m_size;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Long size () const noexcept { return m_size; }
};

template <class T>
struct CommRecvBufTag { // for unpacking recv buffer
    Array4<T> dfab;
    std::ptrdiff_t poff;
    Box bx;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return bx; }
};

template <class T>
struct CommSendBufTag { // for packing send buffer
    Array4<T const> sfab;
    std::ptrdiff_t poff;
    Box bx;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return bx; }
};

namespace detail {

    template <typename T>
    std::enable_if_t<std::is_same_v<std::decay_t<decltype(std::declval<T>().box())>, Box>, Long>
    get_tag_size (T const& tag) noexcept
    {
        AMREX_ASSERT(tag.box().numPts() < Long(std::numeric_limits<int>::max()));
        return static_cast<int>(tag.box().numPts());
    }

    template <typename T>
    std::enable_if_t<std::is_integral_v<std::decay_t<decltype(std::declval<T>().size())> >, Long>
    get_tag_size (T const& tag) noexcept
    {
        AMREX_ASSERT(tag.size() < Long(std::numeric_limits<int>::max()));
        return tag.size();
    }

    //! \cond

    template <typename T>
    constexpr
    std::enable_if_t<std::is_same_v<std::decay_t<decltype(std::declval<T>().box())>, Box>, bool>
    is_box_tag (T const&) { return true; }

    template <typename T>
    constexpr
    std::enable_if_t<std::is_integral_v<std::decay_t<decltype(std::declval<T>().size())> >, bool>
    is_box_tag (T const&) { return false; }

    //! \endcond

}

template <class TagType>
struct TagVector {

    char* h_buffer = nullptr;
    char* d_buffer = nullptr;
    TagType* d_tags =  nullptr;
    int* d_nwarps = nullptr;
    int ntags = 0;
    int ntotwarps = 0;
    int nblocks = 0;
    bool defined = false;
    static constexpr int nthreads = 256;

    TagVector () = default;

    TagVector (Vector<TagType> const& tags) {
        define(tags);
    }

    ~TagVector () {
        if (defined) {
            undefine();
        }
    }

    TagVector (const TagVector& other) = delete;
    TagVector& operator= (const TagVector& other) = delete;
    TagVector (TagVector&& other) noexcept
        : h_buffer{other.h_buffer},
          d_buffer{other.d_buffer},
          d_tags{other.d_tags},
          d_nwarps{other.d_nwarps},
          ntags{other.ntags},
          ntotwarps{other.ntotwarps},
          nblocks{other.nblocks},
          defined{other.defined}
    {
        other.h_buffer = nullptr;
        other.d_buffer = nullptr;
        other.d_tags = nullptr;
        other.d_nwarps = nullptr;
        other.ntags = 0;
        other.ntotwarps = 0;
        other.nblocks = 0;
        other.defined = false;
    }
    TagVector& operator= (TagVector&& other) noexcept {
        if (this == &other) {
            return *this;
        }
        undefine();
        h_buffer = other.h_buffer;
        other.h_buffer = nullptr;
        d_buffer = other.d_buffer;
        other.d_buffer = nullptr;
        d_tags = other.d_tags;
        other.d_tags = nullptr;
        d_nwarps = other.d_nwarps;
        other.d_nwarps = nullptr;
        ntags = other.ntags;
        other.ntags = 0;
        ntotwarps = other.ntotwarps;
        other.ntotwarps = 0;
        nblocks = other.nblocks;
        other.nblocks = 0;
        defined = other.defined;
        other.defined = false;
        return *this;
    }

    [[nodiscard]] bool is_defined () const { return defined; }

    void define (Vector<TagType> const& tags) {
        if (defined) {
            undefine();
        }

        ntags = tags.size();
        if (ntags == 0) {
            defined = true;
            return;
        }

#ifdef AMREX_USE_GPU
        Long l_ntotwarps = 0;
        ntotwarps = 0;
        Vector<int> nwarps;
        nwarps.reserve(ntags+1);
        for (int i = 0; i < ntags; ++i)
        {
            auto& tag = tags[i];
            nwarps.push_back(ntotwarps);
            auto nw = (detail::get_tag_size(tag) + Gpu::Device::warp_size-1) /
                Gpu::Device::warp_size;
            l_ntotwarps += nw;
            ntotwarps += static_cast<int>(nw);
        }
        nwarps.push_back(ntotwarps);

        std::size_t sizeof_tags = ntags*sizeof(TagType);
        std::size_t offset_nwarps = Arena::align(sizeof_tags);
        std::size_t sizeof_nwarps = (ntags+1)*sizeof(int);
        std::size_t total_buf_size = offset_nwarps + sizeof_nwarps;

        h_buffer = (char*)The_Pinned_Arena()->alloc(total_buf_size);
        d_buffer = (char*)The_Arena()->alloc(total_buf_size);

        std::memcpy(h_buffer, tags.data(), sizeof_tags);
        std::memcpy(h_buffer+offset_nwarps, nwarps.data(), sizeof_nwarps);
        Gpu::htod_memcpy_async(d_buffer, h_buffer, total_buf_size);

        d_tags = reinterpret_cast<TagType*>(d_buffer);
        d_nwarps = reinterpret_cast<int*>(d_buffer+offset_nwarps);

        constexpr int nwarps_per_block = nthreads/Gpu::Device::warp_size;
        nblocks = (ntotwarps + nwarps_per_block-1) / nwarps_per_block;

        defined = true;

        amrex::ignore_unused(l_ntotwarps);
        AMREX_ALWAYS_ASSERT(l_ntotwarps+nwarps_per_block-1 < Long(std::numeric_limits<int>::max()));
#else
        std::size_t sizeof_tags = ntags*sizeof(TagType);
        h_buffer = (char*)The_Pinned_Arena()->alloc(sizeof_tags);

        std::memcpy(h_buffer, tags.data(), sizeof_tags);

        d_tags = reinterpret_cast<TagType*>(h_buffer);

        defined = true;
#endif
    }

    void undefine () {
        if (defined) {
            Gpu::streamSynchronize();
            The_Pinned_Arena()->free(h_buffer);
            The_Arena()->free(d_buffer);
            h_buffer = nullptr;
            d_buffer = nullptr;
            d_tags = nullptr;
            d_nwarps = nullptr;
            ntags = 0;
            ntotwarps = 0;
            nblocks = 0;
            defined = false;
        }
    }
};

namespace detail {

#ifdef AMREX_USE_GPU

//! \cond

template <typename T, typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::enable_if_t<std::is_same_v<std::decay_t<decltype(std::declval<T>().box())>, Box>, void>
tagparfor_call_f (int icell, T const& tag, F&& f) noexcept
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
    f(icell, ncells, i, j, k, tag);
}

template <typename T, typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::enable_if_t<std::is_integral_v<std::decay_t<decltype(std::declval<T>().size())> >, void>
tagparfor_call_f (int i, T const& tag, F&& f) noexcept
{
    int N = tag.size();
    f(i, N, tag);
}

//! \endcond

template <class TagType, class F>
void
ParallelFor_doit (TagVector<TagType> const& tv, F const& f)
{
    AMREX_ALWAYS_ASSERT(tv.is_defined());

    if (tv.ntags == 0) { return; }

    const auto d_tags = tv.d_tags;
    const auto d_nwarps = tv.d_nwarps;
    const auto ntags = tv.ntags;
    const auto ntotwarps = tv.ntotwarps;
    constexpr auto nthreads = TagVector<TagType>::nthreads;

    amrex::launch<nthreads>(tv.nblocks, Gpu::gpuStream(),
#ifdef AMREX_USE_SYCL
    [=] AMREX_GPU_DEVICE (sycl::nd_item<1> const& item) noexcept
    [[sycl::reqd_work_group_size(nthreads)]]
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

        tagparfor_call_f(icell, d_tags[tag_id], f);
    });
}

#else // ifdef AMREX_USE_GPU

template <class TagType, class F>
void
ParallelFor_doit (TagVector<TagType> const& tv, F const& f)
{
    // Note: this CPU version may not have optimal performance:
    // The loop over ncomp is the innermost instead of the outermost
    // There is no load-balancing or splitting of tags
    AMREX_ALWAYS_ASSERT(tv.is_defined());

    constexpr bool tag_type = is_box_tag(TagType{});

    if (tv.ntags == 0) { return; }

    const auto d_tags = tv.d_tags;
    const auto ntags = tv.ntags;

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int itag = 0; itag < ntags; ++itag) {

        const auto& t = d_tags[itag];

        if constexpr (tag_type) {
            const auto lo = amrex::lbound(t.box());
            const auto hi = amrex::ubound(t.box());

            for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                        f(0, 1, i, j, k, t);
                    }
                }
            }
        } else {
            const auto size = t.size();

            AMREX_PRAGMA_SIMD
            for (int i = 0; i < size; ++i) {
                f(i, size, t);
            }
        }
    }
}

#endif

}

template <class TagType, class F>
std::enable_if_t<std::is_same_v<std::decay_t<decltype(std::declval<TagType>().box())>, Box>>
ParallelFor (TagVector<TagType> const& tv, int ncomp, F const& f)
{
    detail::ParallelFor_doit(tv,
        [=] AMREX_GPU_DEVICE (
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
std::enable_if_t<std::is_same_v<std::decay_t<decltype(std::declval<TagType>().box())>, Box>, void>
ParallelFor (TagVector<TagType> const& tv, F const& f)
{
    detail::ParallelFor_doit(tv,
        [=] AMREX_GPU_DEVICE (
            int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
        {
            if (icell < ncells) {
                f(i,j,k,tag);
            }
        });
}

template <class TagType, class F>
std::enable_if_t<std::is_integral_v<std::decay_t<decltype(std::declval<TagType>().size())> >, void>
ParallelFor (TagVector<TagType> const& tv, F const& f)
{
    detail::ParallelFor_doit(tv,
        [=] AMREX_GPU_DEVICE (
            int icell, int ncells, TagType const& tag) noexcept
        {
            if (icell < ncells) {
                f(icell,tag);
            }
        });
}

template <class TagType, class F>
std::enable_if_t<std::is_same_v<std::decay_t<decltype(std::declval<TagType>().box())>, Box>>
ParallelFor (Vector<TagType> const& tags, int ncomp, F && f)
{
    TagVector<TagType> tv{tags};
    ParallelFor(tv, ncomp, std::forward<F>(f));
}

template <class TagType, class F>
std::enable_if_t<std::is_same_v<std::decay_t<decltype(std::declval<TagType>().box())>, Box>, void>
ParallelFor (Vector<TagType> const& tags, F && f)
{
    TagVector<TagType> tv{tags};
    ParallelFor(tv, std::forward<F>(f));
}

template <class TagType, class F>
std::enable_if_t<std::is_integral_v<std::decay_t<decltype(std::declval<TagType>().size())> >, void>
ParallelFor (Vector<TagType> const& tags, F && f)
{
    TagVector<TagType> tv{tags};
    ParallelFor(tv, std::forward<F>(f));
}

}

#endif
