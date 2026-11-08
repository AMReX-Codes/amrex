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


#ifdef AMREX_USE_LOOM
/// \cond DOXYGEN_IGNORE
// LOOM device code cannot load pointers from memory. Tags are stored on the
// device with every pointer replaced by its byte offset from one base pointer
// passed as a kernel argument, and rebuilt by value inside the kernel.
namespace loomtag {

struct A4 {
    Long off;
    detail::Stride<3> stride;
    IntVectND<4> begin;
    IntVectND<4> end;
};

template <class T>
A4 store (Array4<T> const& a, char const* base) noexcept {
    return {a.p ? Long((char const*)a.p - base) : 0, a.stride, a.begin, a.end};
}

template <class T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Array4<T> load (A4 const& s, char* base) noexcept {
    return Array4<T>(reinterpret_cast<T*>(base + s.off), s.stride, s.begin, s.end);
}

template <class T>
void note (Array4<T> const& a, char const*& base) noexcept {
    auto* q = (char const*)a.p;
    if (q && (!base || q < base)) { base = q; }
}

// Default: the original layout (tags without a mapping keep their pointers;
// kernels that load them are dropped by the device pass).
template <class Tag> struct Traits {
    using Stored = Tag;
    static void bases (Tag const&, char const*&) {}
    static Stored store (Tag const& t, char const*) { return t; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag const& load (Stored const& s, char*) { return s; }
};
// Box (or size) of a stored tag, read without rebuilding its pointers.
template <class S>
requires requires (S const& s) { s.dbox; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Box stored_box (S const& s) { return s.dbox; }
template <class S>
requires (requires (S const& s) { s.bx; } && !requires (S const& s) { s.dbox; })
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Box stored_box (S const& s) { return s.bx; }
template <class S>
requires (requires (S const& s) { s.box(); } && !requires (S const& s) { s.bx; } &&
          !requires (S const& s) { s.dbox; })
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Box stored_box (S const& s) { return s.box(); }

template <class T0, class T1>
struct Traits<Array4CopyTag<T0,T1>> {
    using Tag = Array4CopyTag<T0,T1>;
    struct Stored { A4 dfab; int dindex; A4 sfab; Box dbox; Dim3 offset; };
    static void bases (Tag const& t, char const*& b) { note(t.dfab,b); note(t.sfab,b); }
    static Stored store (Tag const& t, char const* b) {
        return {loomtag::store(t.dfab,b), t.dindex, loomtag::store(t.sfab,b), t.dbox, t.offset};
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag load (Stored const& s, char* b) {
        return Tag{loomtag::load<T0>(s.dfab,b), s.dindex, loomtag::load<T1 const>(s.sfab,b),
                   s.dbox, s.offset};
    }
};

template <class T0, class T1>
struct Traits<Array4MaskCopyTag<T0,T1>> {
    using Tag = Array4MaskCopyTag<T0,T1>;
    struct Stored { A4 dfab; A4 sfab; A4 mask; Box dbox; Dim3 offset; };
    static void bases (Tag const& t, char const*& b) { note(t.dfab,b); note(t.sfab,b); note(t.mask,b); }
    static Stored store (Tag const& t, char const* b) {
        return {loomtag::store(t.dfab,b), loomtag::store(t.sfab,b), loomtag::store(t.mask,b),
                t.dbox, t.offset};
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag load (Stored const& s, char* b) {
        return Tag{loomtag::load<T0>(s.dfab,b), loomtag::load<T1 const>(s.sfab,b),
                   loomtag::load<int>(s.mask,b), s.dbox, s.offset};
    }
};

template <class T>
struct Traits<Array4PairTag<T>> {
    using Tag = Array4PairTag<T>;
    struct Stored { A4 dfab; A4 sfab; Box dbox; };
    static void bases (Tag const& t, char const*& b) { note(t.dfab,b); note(t.sfab,b); }
    static Stored store (Tag const& t, char const* b) {
        return {loomtag::store(t.dfab,b), loomtag::store(t.sfab,b), t.dbox};
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag load (Stored const& s, char* b) {
        return Tag{loomtag::load<T>(s.dfab,b), loomtag::load<T const>(s.sfab,b), s.dbox};
    }
};

template <class T>
struct Traits<Array4Tag<T>> {
    using Tag = Array4Tag<T>;
    struct Stored {
        A4 dfab;
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Box box () const noexcept {
            return Box(IntVect(dfab.begin[0], dfab.begin[1], dfab.begin[2]),
                       IntVect(dfab.end[0]-1, dfab.end[1]-1, dfab.end[2]-1));
        }
    };
    static void bases (Tag const& t, char const*& b) { note(t.dfab,b); }
    static Stored store (Tag const& t, char const* b) { return {loomtag::store(t.dfab,b)}; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag load (Stored const& s, char* b) { return Tag{loomtag::load<T>(s.dfab,b)}; }
};

template <class T>
struct Traits<Array4BoxTag<T>> {
    using Tag = Array4BoxTag<T>;
    struct Stored { A4 dfab; Box dbox; };
    static void bases (Tag const& t, char const*& b) { note(t.dfab,b); }
    static Stored store (Tag const& t, char const* b) { return {loomtag::store(t.dfab,b), t.dbox}; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag load (Stored const& s, char* b) { return Tag{loomtag::load<T>(s.dfab,b), s.dbox}; }
};

template <class T>
struct Traits<Array4BoxValTag<T>> {
    using Tag = Array4BoxValTag<T>;
    struct Stored { A4 dfab; Box dbox; T val; };
    static void bases (Tag const& t, char const*& b) { note(t.dfab,b); }
    static Stored store (Tag const& t, char const* b) { return {loomtag::store(t.dfab,b), t.dbox, t.val}; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag load (Stored const& s, char* b) { return Tag{loomtag::load<T>(s.dfab,b), s.dbox, s.val}; }
};

template <class T>
struct Traits<Array4BoxOrientationTag<T>> {
    using Tag = Array4BoxOrientationTag<T>;
    struct Stored { A4 fab; Box bx; Orientation face; };
    static void bases (Tag const& t, char const*& b) { note(t.fab,b); }
    static Stored store (Tag const& t, char const* b) { return {loomtag::store(t.fab,b), t.bx, t.face}; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag load (Stored const& s, char* b) { return Tag{loomtag::load<T>(s.fab,b), s.bx, s.face}; }
};

template <class T>
struct Traits<Array4BoxOffsetTag<T>> {
    using Tag = Array4BoxOffsetTag<T>;
    struct Stored { A4 fab; Box bx; Dim3 offset; };
    static void bases (Tag const& t, char const*& b) { note(t.fab,b); }
    static Stored store (Tag const& t, char const* b) { return {loomtag::store(t.fab,b), t.bx, t.offset}; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag load (Stored const& s, char* b) { return Tag{loomtag::load<T>(s.fab,b), s.bx, s.offset}; }
};

template <class T>
struct Traits<CommRecvBufTag<T>> {
    using Tag = CommRecvBufTag<T>;
    struct Stored { A4 dfab; std::ptrdiff_t poff; Box bx; };
    static void bases (Tag const& t, char const*& b) { note(t.dfab,b); }
    static Stored store (Tag const& t, char const* b) { return {loomtag::store(t.dfab,b), t.poff, t.bx}; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag load (Stored const& s, char* b) { return Tag{loomtag::load<T>(s.dfab,b), s.poff, s.bx}; }
};

template <class T>
struct Traits<CommSendBufTag<T>> {
    using Tag = CommSendBufTag<T>;
    struct Stored { A4 sfab; std::ptrdiff_t poff; Box bx; };
    static void bases (Tag const& t, char const*& b) { note(t.sfab,b); }
    static Stored store (Tag const& t, char const* b) { return {loomtag::store(t.sfab,b), t.poff, t.bx}; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag load (Stored const& s, char* b) { return Tag{loomtag::load<T const>(s.sfab,b), s.poff, s.bx}; }
};

template <class T>
struct Traits<VectorTag<T>> {
    using Tag = VectorTag<T>;
    struct Stored { Long off; Long m_size; };
    static void bases (Tag const& t, char const*& b) {
        auto* q = (char const*)t.p;
        if (q && (!b || q < b)) { b = q; }
    }
    static Stored store (Tag const& t, char const* b) {
        return {t.p ? Long((char const*)t.p - b) : 0, t.m_size};
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Tag load (Stored const& s, char* b) { return Tag{reinterpret_cast<T*>(b + s.off), s.m_size}; }
};

}
/// \endcond
#endif

/// \cond DOXYGEN_IGNORE
namespace detail {

    template <typename T>
    requires (std::same_as<std::decay_t<decltype(std::declval<T>().box())>, Box>)
    Long get_tag_size (T const& tag) noexcept
    {
        AMREX_ASSERT(tag.box().numPts() < Long(std::numeric_limits<int>::max()));
        return static_cast<int>(tag.box().numPts());
    }

    template <typename T>
    requires (std::is_integral_v<std::decay_t<decltype(std::declval<T>().size())> >)
    Long get_tag_size (T const& tag) noexcept
    {
        AMREX_ASSERT(tag.size() < Long(std::numeric_limits<int>::max()));
        return tag.size();
    }

    template <typename T>
    requires (std::same_as<std::decay_t<decltype(std::declval<T>().box())>, Box>)
    constexpr bool is_box_tag (T const&) { return true; }

    template <typename T>
    requires (std::is_integral_v<std::decay_t<decltype(std::declval<T>().size())> >)
    constexpr bool is_box_tag (T const&) { return false; }

}
/// \endcond

template <class TagType>
struct TagVector {

    char* h_buffer = nullptr;
    char* d_buffer = nullptr;
    TagType* d_tags =  nullptr;
    int* d_nwarps = nullptr;
#ifdef AMREX_USE_LOOM
    char* loom_base = nullptr; // tags hold pointer offsets from this base
#endif
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
#ifdef AMREX_USE_LOOM
          loom_base{other.loom_base},
#endif
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
#ifdef AMREX_USE_LOOM
        loom_base = other.loom_base;
#endif
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

#ifdef AMREX_USE_LOOM
        using LoomTraits = loomtag::Traits<TagType>;
        using LoomStored = typename LoomTraits::Stored;
        char const* lbase = nullptr;
        for (auto const& tag : tags) { LoomTraits::bases(tag, lbase); }
        loom_base = const_cast<char*>(lbase);
        Vector<LoomStored> stored;
        stored.reserve(ntags);
        for (auto const& tag : tags) { stored.push_back(LoomTraits::store(tag, lbase)); }
        std::size_t sizeof_tags = ntags*sizeof(LoomStored);
#else
        std::size_t sizeof_tags = ntags*sizeof(TagType);
#endif
        std::size_t offset_nwarps = Arena::align(sizeof_tags);
        std::size_t sizeof_nwarps = (ntags+1)*sizeof(int);
        std::size_t total_buf_size = offset_nwarps + sizeof_nwarps;

        h_buffer = (char*)The_Pinned_Arena()->alloc(total_buf_size);
        d_buffer = (char*)The_Arena()->alloc(total_buf_size);

#ifdef AMREX_USE_LOOM
        std::memcpy(h_buffer, stored.data(), sizeof_tags);
#else
        std::memcpy(h_buffer, tags.data(), sizeof_tags);
#endif
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

/// \cond DOXYGEN_IGNORE
namespace detail {

#ifdef AMREX_USE_GPU

template <typename T, typename F>
requires (std::same_as<std::decay_t<decltype(std::declval<T>().box())>, Box>)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void tagparfor_call_f (int icell, T const& tag, F&& f) noexcept
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
requires (std::is_integral_v<std::decay_t<decltype(std::declval<T>().size())> >)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void tagparfor_call_f (int i, T const& tag, F&& f) noexcept
{
    int N = tag.size();
    f(i, N, tag);
}

#ifdef AMREX_USE_LOOM
template <class Traits, class S, class F>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void loom_tagparfor_call_f (int icell, S const& s, char* base, F const& f) noexcept
{
    if constexpr (requires { loomtag::stored_box(s); }) {
        const Box bx = loomtag::stored_box(s);
        int ncells = bx.numPts();
        const auto len = amrex::length(bx);
        const auto lo  = amrex::lbound(bx);
        int k =  icell /   (len.x*len.y);
        int j = (icell - k*(len.x*len.y)) /   len.x;
        int i = (icell - k*(len.x*len.y)) - j*len.x;
        i += lo.x;
        j += lo.y;
        k += lo.z;
        f(icell, ncells, i, j, k, Traits::load(s, base));
    } else {
        int N = static_cast<int>(s.m_size);
        f(icell, N, Traits::load(s, base));
    }
}
#endif

template <class TagType, class F>
void
ParallelFor_doit (TagVector<TagType> const& tv, F const& f)
{
    AMREX_ALWAYS_ASSERT(tv.is_defined());

    if (tv.ntags == 0 || tv.ntotwarps == 0) { return; }

#ifdef AMREX_USE_LOOM
    using LoomTraits = loomtag::Traits<TagType>;
    const auto* d_tags = reinterpret_cast<typename LoomTraits::Stored const*>(tv.d_tags);
    char* loom_base = tv.loom_base;
#else
    const auto d_tags = tv.d_tags;
#endif
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

#ifdef AMREX_USE_LOOM
        loom_tagparfor_call_f<LoomTraits>(icell, d_tags[tag_id], loom_base, f);
#else
        tagparfor_call_f(icell, d_tags[tag_id], f);
#endif
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
/// \endcond

template <class TagType, class F>
requires (std::same_as<std::decay_t<decltype(std::declval<TagType>().box())>, Box>)
void ParallelFor (TagVector<TagType> const& tv, int ncomp, F const& f)
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
requires (std::same_as<std::decay_t<decltype(std::declval<TagType>().box())>, Box>)
void ParallelFor (TagVector<TagType> const& tv, F const& f)
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
requires (std::is_integral_v<std::decay_t<decltype(std::declval<TagType>().size())> >)
void ParallelFor (TagVector<TagType> const& tv, F const& f)
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
requires (std::same_as<std::decay_t<decltype(std::declval<TagType>().box())>, Box>)
void ParallelFor (Vector<TagType> const& tags, int ncomp, F && f)
{
    TagVector<TagType> tv{tags};
    ParallelFor(tv, ncomp, std::forward<F>(f));
}

template <class TagType, class F>
requires (std::same_as<std::decay_t<decltype(std::declval<TagType>().box())>, Box>)
void ParallelFor (Vector<TagType> const& tags, F && f)
{
    TagVector<TagType> tv{tags};
    ParallelFor(tv, std::forward<F>(f));
}

template <class TagType, class F>
requires (std::is_integral_v<std::decay_t<decltype(std::declval<TagType>().size())> >)
void ParallelFor (Vector<TagType> const& tags, F && f)
{
    TagVector<TagType> tv{tags};
    ParallelFor(tv, std::forward<F>(f));
}

}

#endif
