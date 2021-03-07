#ifndef AMREX_NONLOCAL_BC_IMPL_H_
#define AMREX_NONLOCAL_BC_IMPL_H_
#include <AMReX_Config.H>
#include <AMReX_TypeTraits.H>

#ifndef AMREX_STATIC_ASSERT_NO_MESSAGE
#  define AMREX_STATIC_ASSERT_NO_MESSAGE(expr) static_assert((expr), #expr)
#  define AMREX_HAS_STATIC_ASSERT_NO_MESSAGE 0
#else
#  define AMREX_HAS_STATIC_ASSERT_NO_MESSAGE 1
#endif

namespace amrex { namespace NonLocalBC {

/// Tests if a given type `IndexMap` is usable as an index mapping between
/// two index based coordinate systems.
template <class IndexMap>
struct IsIndexMapping : Conjunction<IsCallableR<IntVect, IndexMap&, const IntVect&>,
                                    IsCallableR<Dim3, IndexMap&, const Dim3&>,
                                    IsCallableR<Box, IndexMap&, const Box&>> {};

struct Rotate90ClockWise {
    AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(iv[1], -1-iv[0], iv[2])};
    }

    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{a.y, -1-a.x, a.z};
    }

    Box operator() (Box const& box) const noexcept {
        return Box(operator()(IntVect{AMREX_D_DECL(box.bigEnd  (0),
                                                   box.smallEnd(1),
                                                   box.smallEnd(2))}),
                   operator()(IntVect{AMREX_D_DECL(box.smallEnd(0),
                                                   box.bigEnd  (1),
                                                   box.bigEnd  (2))}));
    }
};
static_assert(IsIndexMapping<Rotate90ClockWise>(),
              "Rotate90ClockWise is expected to satisfy IndexMapping");

struct Rotate90CounterClockWise {
    AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(-1-iv[1], iv[0], iv[2])};
    }

    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{-1-a.y, a.x, a.z};
    }

    Box operator() (Box const& box) const noexcept {
        return Box(operator()(IntVect{AMREX_D_DECL(box.smallEnd(0),
                                                   box.bigEnd  (1),
                                                   box.smallEnd(2))}),
                   operator()(IntVect{AMREX_D_DECL(box.bigEnd  (0),
                                                   box.smallEnd(1),
                                                   box.bigEnd  (2))}));
    }
};
static_assert(IsIndexMapping<Rotate90CounterClockWise>(),
              "Rotate90CounterClockWise is expected to satisfy IndexMapping");

struct Rotate90DstToSrc
{
    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        if (a.x < 0) {
            return Rotate90ClockWise()(a);
        } else {
            return Rotate90CounterClockWise()(a);
        }
    }
};

template <typename Base, bool> struct InverseMappingFacade {
    AMREX_NO_UNIQUE_ADDRESS Base base;
};
template<typename Base> struct InverseMappingFacade<Base, true> {
    constexpr InverseMappingFacade() = default;
    constexpr InverseMappingFacade(const Base& b) : base(b) {}
    constexpr InverseMappingFacade(Base&& b) : base(std::move(b)) {}

    AMREX_NO_UNIQUE_ADDRESS Base base;

    static constexpr bool is_noexcept_invertible = noexcept(std::declval<Base const&>().Inverse(Dim3{}));

    AMREX_GPU_HOST_DEVICE
    IntVect Inverse(IntVect const& iv) const noexcept(is_noexcept_invertible) {
        Dim3 i{};
        AMREX_D_DECL(i.x = iv[0], i.y = iv[1], i.z = iv[2]);
        i = base.Inverse(i);
        return {AMREX_D_DECL(i.x, i.y, i.z)};
    }
    
    Box InverseImage(const Box& box) const noexcept(is_noexcept_invertible) {
        IntVect mapped_smallEnd = Inverse(box.smallEnd());
        IntVect mapped_bigEnd = Inverse(box.bigEnd());
        IntVect smallEnd;
        IntVect bigEnd;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            smallEnd[d] = std::min(mapped_smallEnd[d], mapped_bigEnd[d]);
            bigEnd[d] = std::max(mapped_smallEnd[d], mapped_bigEnd[d]);
        }
        return Box{smallEnd, bigEnd};
    }
};

template <typename I, typename... Args>
using Inverse_t = decltype(std::declval<I>().Inverse(std::declval<Args>()...));

/// This facade lifts a given base class from being Dim3 -> Dim3 callable to an  index map
/// that satisfies our concept.
template <typename Base>
struct IndexMappingFacade : InverseMappingFacade<Base, IsDetectedExact<Dim3, Inverse_t, Base, Dim3>::value> {
    using BaseFacade = InverseMappingFacade<Base, IsDetectedExact<Dim3, Inverse_t, Base, Dim3>::value>;
    
    using InverseMappingFacade<Base, IsDetectedExact<Dim3, Inverse_t, Base, Dim3>::value>::InverseMappingFacade;
    
    static_assert(IsCallableR<Dim3, Base, Dim3>(),
                  "Base class needs to be a mapping: Dim3 -> Dim3");

    static constexpr bool is_noexcept_callable = noexcept(std::declval<Base const&>()(Dim3{}));

    AMREX_GPU_HOST_DEVICE
    IntVect operator()(IntVect const& iv) const noexcept(is_noexcept_callable) {
        Dim3 i{};
        AMREX_D_DECL(i.x = iv[0], i.y = iv[1], i.z = iv[2]);
        i = BaseFacade::base(i);
        return {AMREX_D_DECL(i.x, i.y, i.z)};
    }

    AMREX_GPU_HOST_DEVICE
    constexpr Dim3 operator()(Dim3 i) const noexcept(is_noexcept_callable) { return BaseFacade::base(i); }

    //! Performs the index mapping on each member of the box and returns a valid
    //! box.
    //!
    //! \return Returns the smallest box that contains each mapped index.
    Box operator()(const Box& box) const noexcept(is_noexcept_callable) {
        IntVect mapped_smallEnd = this->operator()(box.smallEnd());
        IntVect mapped_bigEnd = this->operator()(box.bigEnd());
        IntVect smallEnd;
        IntVect bigEnd;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            smallEnd[d] = std::min(mapped_smallEnd[d], mapped_bigEnd[d]);
            bigEnd[d] = std::max(mapped_smallEnd[d], mapped_bigEnd[d]);
        }
        return Box{smallEnd, bigEnd};
    }
};
using Rotate90IndexMapping = IndexMappingFacade<Rotate90DstToSrc>;
static_assert(IsIndexMapping<Rotate90IndexMapping>(),
              "Rotate90IndexMapping is expected to satisfy IndexMapping");

struct Rotate180Fn {
    int Ly;

    AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(-1-iv[0], Ly-1-iv[1], iv[2])};
    }

    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{-1-a.x, Ly-1-a.y, a.z};
    }

    Box operator() (Box const& box) const noexcept {
        return Box(operator()(IntVect{AMREX_D_DECL(box.bigEnd  (0),
                                                   box.bigEnd  (1),
                                                   box.smallEnd(2))}),
                   operator()(IntVect{AMREX_D_DECL(box.smallEnd(0),
                                                   box.smallEnd(1),
                                                   box.bigEnd  (2))}));
    }
};
static_assert(IsIndexMapping<Rotate180Fn>(),
              "Rotate180Fn is expected to satisfy IndexMapping");

struct PolarFn {
    int Lx, Ly;

    AMREX_GPU_HOST_DEVICE
    int i_index (int i) const noexcept {
        return (i < Lx/2) ? -1-i : 2*Lx-1-i;
    }

    AMREX_GPU_HOST_DEVICE
    int j_index (int j) const noexcept {
        return (j < Ly/2) ? j+Ly/2 : j-Ly/2;
    }

    AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(i_index(iv[0]), j_index(iv[1]), iv[2])};
    }

    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{i_index(a.x), j_index(a.y), a.z};
    }

    Box operator() (Box const& box) const noexcept {
        return Box(operator()(IntVect{AMREX_D_DECL(box.bigEnd  (0),
                                                   box.smallEnd(1),
                                                   box.smallEnd(2))}),
                   operator()(IntVect{AMREX_D_DECL(box.smallEnd(0),
                                                   box.bigEnd  (1),
                                                   box.bigEnd  (2))}));
    }
};
static_assert(IsIndexMapping<PolarFn>(),
              "PolarFn is expected to satisfy IndexMapping");

struct PolarFn2 { // for the x-y corners
    int Lx, Ly;

    AMREX_GPU_HOST_DEVICE
    int i_index (int i) const noexcept {
        return (i < Lx/2) ? -1-i : 2*Lx-1-i;
    }

    AMREX_GPU_HOST_DEVICE
    int j_index (int j) const noexcept {
        if (j < 0) {
            return j+Ly/2;
        } else if (j >= Ly) {
            return j-Ly/2;
        } else if (j < Ly/2) {
            return j-Ly/2;
        } else {
            return j+Ly/2;
        }
    }

    AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(i_index(iv[0]), j_index(iv[1]), iv[2])};
    }

    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{i_index(a.x), j_index(a.y), a.z};
    }

    Box operator() (Box const& box) const noexcept {
        return Box(operator()(IntVect{AMREX_D_DECL(box.bigEnd  (0),
                                                   box.smallEnd(1),
                                                   box.smallEnd(2))}),
                   operator()(IntVect{AMREX_D_DECL(box.smallEnd(0),
                                                   box.bigEnd  (1),
                                                   box.bigEnd  (2))}));
    }
};
static_assert(IsIndexMapping<PolarFn2>(),
              "PolarFn2 is expected to satisfy IndexMapping");


//! This struct describes a bijective index transformation for two coordinate
//! systems.
struct MultiBlockDestToSrc {
    IntVect permutation{AMREX_D_DECL(0, 1, 2)};
    IntVect offset{AMREX_D_DECL(0, 0, 0)};
    IntVect sign{AMREX_D_DECL(1, 1, 1)};

    //! Applies the index mapping on (i, j, k) and returns an index in the source
    //! space.
    AMREX_GPU_HOST_DEVICE Dim3 operator()(Dim3 i) const noexcept {
        int iv[3]{i.x, i.y, i.z};
        int iv_new[3]{};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            iv_new[d] = sign[d] * (iv[permutation[d]] - offset[d]);
        }
        return {iv_new[0], iv_new[1], iv_new[2]};
    }

    Dim3 Inverse(Dim3 i) const noexcept {
        int iv_new[3]{i.x, i.y, i.z};
        int iv[3]{};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            iv[permutation[d]] = iv_new[d] * sign[d] + offset[d];
        }
        return {iv[0], iv[1], iv[2]};
    }
};
using MultiBlockIndexMapping = IndexMappingFacade<MultiBlockDestToSrc>;
static_assert(IsIndexMapping<MultiBlockIndexMapping>(),
              "MultiBlockIndexMapping is expected to satisfy IndexMapping");

template <typename P, typename FAB>
struct IsFabProjection
    : IsCallableR<typename FAB::value_type, P, Array4<const typename FAB::value_type>, Dim3, int> 
{};

struct IdentityDestToSrc {
    constexpr Dim3 operator()(Dim3 i) const noexcept { return i; }
    constexpr Dim3 Inverse(Dim3 i) const noexcept { return i; }
};

/// This type acts as a default no-op operator 
struct Identity : IndexMappingFacade<IdentityDestToSrc> {
    using IndexMappingFacade<IdentityDestToSrc>::operator();

    template <typename T>
    constexpr T operator()(Array4<const T> array, Dim3 i, int comp = 0) const
        noexcept(noexcept(array(i.x, i.y, i.z, comp))) {
        return array(i.x, i.y, i.z, comp);
    }

    constexpr int operator()(int i) const noexcept { return i; }
};
static constexpr Identity identity{};

AMREX_STATIC_ASSERT_NO_MESSAGE((sizeof(Identity) == 1));
AMREX_STATIC_ASSERT_NO_MESSAGE((std::is_trivially_default_constructible<Identity>::value));
AMREX_STATIC_ASSERT_NO_MESSAGE((std::is_trivially_copy_assignable<Identity>::value));
AMREX_STATIC_ASSERT_NO_MESSAGE((std::is_trivially_copy_constructible<Identity>::value));
AMREX_STATIC_ASSERT_NO_MESSAGE((IsIndexMapping<Identity>()));
AMREX_STATIC_ASSERT_NO_MESSAGE((IsFabProjection<Identity, FArrayBox>()));

template <typename Base, typename Map = Identity> struct MapComponents {
    static_assert(IsCallable<Base, Array4<const Real>, Dim3, int>::value,
                  "Base needs to be a callable function: (Array4<const T>, Dim3, i) -> auto.");

    static_assert(IsCallableR<int, Map, int>::value,
                  "Map needs to be a callable function: int -> int.");

    AMREX_NO_UNIQUE_ADDRESS Base base;
    AMREX_NO_UNIQUE_ADDRESS Map map;

    template <typename T,
              typename = EnableIf_t<IsCallable<Base, Array4<const T>, Dim3, int>::value>,
              typename = EnableIf_t<IsCallableR<int, Map, int>::value>>
    constexpr decltype(auto) operator()(Array4<const T> array, Dim3 i, int comp) const
        noexcept(noexcept(base(array, i, map(comp)))) {
        return base(array, i, map(comp));
    }
};
#ifdef AMREX_HAS_NO_UNIQUE_ADDRESS
AMREX_STATIC_ASSERT_NO_MESSAGE((sizeof(MapComponents<Identity, Identity>) == 1));
#endif
AMREX_STATIC_ASSERT_NO_MESSAGE((std::is_trivially_copy_assignable<MapComponents<Identity>>()));
AMREX_STATIC_ASSERT_NO_MESSAGE((std::is_trivially_copy_constructible<MapComponents<Identity>>()));

#ifdef AMREX_DEFINE_SWAP_INDICES_OPERATOR
#error "The AMREX_DEFINE_SWAP_INDICES_OPERATOR macro is already defined"
#endif
#define AMREX_DEFINE_SWAP_INDICES_OPERATOR                                                         \
    constexpr int operator()(int i) const noexcept {                                               \
        const int map[2] = {I, J};                                                                 \
        return i == I || i == J ? map[std::size_t(i == I)] : i;                                    \
    }

template <int I, int J> struct SwapIndices { AMREX_DEFINE_SWAP_INDICES_OPERATOR };

template <int I> struct SwapIndices<I, -1> {
    int J;
    AMREX_DEFINE_SWAP_INDICES_OPERATOR
};

template <int J> struct SwapIndices<-1, J> {
    int I;
    AMREX_DEFINE_SWAP_INDICES_OPERATOR
};

template <> struct SwapIndices<-1, -1> {
    int I;
    int J;
    AMREX_DEFINE_SWAP_INDICES_OPERATOR
};
#undef AMREX_DEFINE_SWAP_INDICES_OPERATOR

template <int I, int J> static constexpr SwapIndices<I, J> swap_indices{};
using SwapComponents = SwapIndices<-1, -1>;

AMREX_STATIC_ASSERT_NO_MESSAGE(sizeof(SwapIndices<0, 1>) == 1);
AMREX_STATIC_ASSERT_NO_MESSAGE(sizeof(SwapComponents) == 2 * sizeof(int));
AMREX_STATIC_ASSERT_NO_MESSAGE(sizeof(SwapIndices<0, -1>) == sizeof(int));
AMREX_STATIC_ASSERT_NO_MESSAGE(sizeof(SwapIndices<-1, 1>) == sizeof(int));
#ifdef AMREX_HAS_NO_UNIQUE_ADDRESS
AMREX_STATIC_ASSERT_NO_MESSAGE(sizeof(MapComponents<Identity, SwapIndices<0, 1>>) == sizeof(SwapIndices<0, 1>));
#endif
AMREX_STATIC_ASSERT_NO_MESSAGE((std::is_trivially_default_constructible<MapComponents<Identity, SwapIndices<0, 1>>>()));
AMREX_STATIC_ASSERT_NO_MESSAGE((std::is_trivially_copy_assignable<MapComponents<Identity, SwapIndices<0, 1>>>()));
AMREX_STATIC_ASSERT_NO_MESSAGE((std::is_trivially_copy_constructible<MapComponents<Identity, SwapIndices<0, 1>>>()));
AMREX_STATIC_ASSERT_NO_MESSAGE((IsFabProjection<MapComponents<Identity, SwapIndices<0, 1>>, FArrayBox>()));

AMREX_STATIC_ASSERT_NO_MESSAGE((swap_indices<0, 1>(0) == 1));
AMREX_STATIC_ASSERT_NO_MESSAGE((swap_indices<0, 1>(1) == 0));
AMREX_STATIC_ASSERT_NO_MESSAGE((swap_indices<0, 1>(2) == 2));
AMREX_STATIC_ASSERT_NO_MESSAGE((SwapComponents{0, 1}(0) == 1));
AMREX_STATIC_ASSERT_NO_MESSAGE((SwapComponents{0, 1}(1) == 0));
AMREX_STATIC_ASSERT_NO_MESSAGE((SwapComponents{0, 1}(2) == 2));

template <class FAB, class DTOS = Identity, class Proj = Identity>
EnableIf_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
local_copy_cpu (FabArray<FAB>& dest, const FabArray<FAB>& src, int dcomp, int scomp, int ncomp,
                FabArrayBase::CopyComTagsContainer const& local_tags, DTOS dtos = DTOS{},
                Proj proj = Proj{}) noexcept {
    const int N_locs = local_tags.size();
    if (N_locs == 0) return;
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int itag = 0; itag < N_locs; ++itag) {
        const auto& tag = LocTags[itag];
        auto const& sfab = mf.const_array(tag.srcIndex);
        auto const& dfab = mf.array      (tag.dstIndex);
        amrex::LoopConcurrentOnCpu(tag.dbox, ncomp, [=] (int i, int j, int k, int n) noexcept
        {
            auto const si = dtos(Dim3{i,j,k});
            dfab(i,j,k,dcomp+n) = proj(sfab,si,scomp+n);
        });
    }
}

template <class FAB, class DTOS = Identity, class Proj = Identity>
EnableIf_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
unpack_recv_buffer_cpu (FabArray<FAB>& mf, int dcomp, int ncomp, Vector<char*> const& recv_data,
                        Vector<std::size_t> const& recv_size,
                        Vector<FabArrayBase::CopyComTagsContainer const*> const& recv_cctc,
                        DTOS dtos = DTOS{}, Proj proj = Proj{}) noexcept {
    amrex::ignore_unused(recv_size);

    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) return;

    using T = typename FAB::value_type;
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int ircv = 0; ircv < N_rcvs; ++ircv) {
        const char* dptr = recv_data[ircv];
        auto const& cctc = *recv_cctc[ircv];
        for (auto const& tag : cctc) {
            auto const& dfab = mf.array(tag.dstIndex);
            auto const& sfab = amrex::makeArray4((T const*)(dptr), tag.sbox, ncomp);
            amrex::LoopConcurrentOnCpu(tag.dbox, ncomp, [=](int i, int j, int k, int n) noexcept {
                auto const si = dtos(Dim3{i, j, k});
                dfab(i, j, k, dcomp + n) = proj(sfab, si, n);
            });
            dptr += tag.sbox.numPts() * ncomp * sizeof(T);
            AMREX_ASSERT(dptr <= recv_data[ircv] + recv_size[ircv]);
        }
    }
}

#ifdef AMREX_USE_GPU
template <class T>
struct Array4Array4Box {
    Array4<T      > dfab;
    Array4<T const> sfab;
    Box dbox;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return dbox; }
};

template <class FAB, class DTOS = Identity, class Proj = Identity>
EnableIf_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
local_copy_gpu (FabArray<FAB>& dest, const FabArray<FAB>& src, int dcomp, int scomp, int ncomp,
                FabArrayBase::CopyComTagsContainer const& local_tags, DTOS dtos = DTOS{},
                Proj proj = Proj{}) noexcept {
    int N_locs = local_tags.size();
    if (N_locs == 0) return;

    using T = typename FAB::value_type;
    Vector<Array4Array4Box<T> > loc_copy_tags;
    loc_copy_tags.reserve(N_locs);
    for (auto const& tag : local_tags) {
        loc_copy_tags.push_back({dest.array(tag.dstIndex), src.const_array(tag.srcIndex), tag.dbox});
    }

    ParallelFor(loc_copy_tags, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                                                            Array4Array4Box<T> const& tag) noexcept
    {
        auto const si = dtos(Dim3{i,j,k});
        tag.dfab(i,j,k,dcomp+n) = proj(tag.sfab, si, scomp+n);
    });
}

template <class FAB, class DTOS = Identity, class Proj = Identity>
EnableIf_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
unpack_recv_buffer_gpu (FabArray<FAB>& mf, int scomp, int ncomp,
                        Vector<char*> const& recv_data,
                        Vector<std::size_t> const& recv_size,
                        Vector<FabArrayBase::CopyComTagsContainer const*> const& recv_cctc,
                        DTOS dtos = DTOS{}, Proj proj = Proj{})
{
    amrex::ignore_unused(recv_size);

    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) return;

    char* pbuffer = recv_data[0];
#if 0
    std::size_t szbuffer = 0;
    // For linear solver test on summit, this is slower than writing to
    // pinned memory directly on device.
    if (not ParallelDescriptor::UseGpuAwareMpi()) {
        // Memory in recv_data is pinned.
        szbuffer = (recv_data[N_rcvs-1]-recv_data[0]) + recv_size[N_rcvs-1];
        pbuffer = (char*)The_Arena()->alloc(szbuffer);
        Gpu::copyAsync(Gpu::hostToDevice,recv_data[0],recv_data[0]+szbuffer,pbuffer);
        Gpu::synchronize();
    }
#endif

    using T = typename FAB::value_type;
    using TagType = Array4Array4Box<T>;
    Vector<TagType> tags;
    tags.reserve(N_rcvs);

    for (int k = 0; k < N_rcvs; ++k)
    {
        std::size_t offset = recv_data[k]-recv_data[0];
        const char* dptr = pbuffer + offset;
        auto const& cctc = *recv_cctc[k];
        for (auto const& tag : cctc)
        {
            tags.emplace_back(TagType{mf.array(tag.dstIndex),
                                      amrex::makeArray4((T const*)dptr, tag.sbox, ncomp),
                                      tag.dbox});
            dptr += tag.dbox.numPts() * ncomp * sizeof(T);
            BL_ASSERT(dptr <= pbuffer + offset + recv_size[k]);
        }
    }

    ParallelFor(tags, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                                                   Array4Array4Box<T> const& tag) noexcept
    {
        auto const si = dtos(Dim3{i,j,k});
        tag.dfab(i,j,k,scomp+n) = proj(tag.sfab, si ,n);
    });

    // There is Gpu::synchronize in ParalleFor above

    if (pbuffer != recv_data[0]) {
        The_Arena()->free(pbuffer);
    }
}
#endif

struct CommData {
#ifdef AMREX_USE_MPI
    TheFaArenaPointer the_data = nullptr;
    Vector<int> rank{};
    Vector<char*> data{};
    Vector<std::size_t> offset{};
    Vector<std::size_t> size{};
    Vector<MPI_Request> request{};
    Vector<MPI_Status> stats{};
    Vector<const FabArrayBase::CopyComTagsContainer*> cctc{};
#endif
};

struct CommHandler {
#ifdef AMREX_USE_MPI
    int mpi_tag{};
    CommData recv{};
    CommData send{};
#endif
};

#ifdef AMREX_USE_MPI
void PostRecvs(CommData& comm, int mpi_tag);
void PostSends(CommData& comm, int mpi_tag);
#endif

////////////////////////////////////////////////////////////////////////////////////
//                                                             [concept.DataPacking]
//
template <typename... Args>
using PrepareSendBuffers_t = decltype(PrepareSendBuffers(std::declval<Args>()...));

template <typename... Args>
using PrepareRecvBuffers_t = decltype(PrepareRecvBuffers(std::declval<Args>()...));

template <typename... Args>
using PackSendBuffers_t = decltype(PackSendBuffers(std::declval<Args>()...));

template <typename... Args>
using UnpackRecvBuffers_t = decltype(UnpackRecvBuffers(std::declval<Args>()...));

template <typename... Args>
using LocalCopy_t = decltype(LocalCopy(std::declval<Args>()...));

template <typename DP, typename FAB>
struct IsDataPacking : 
    Conjunction<
        IsDetected<LocalCopy_t, DP&, FabArray<FAB>&, const FabArray<FAB>&, const FabArrayBase::CopyComTagsContainer&>
#ifdef AMREX_USE_MPI
       ,IsDetected<PrepareSendBuffers_t, DP&, FabArray<FAB>&, const FabArray<FAB>&, CommData&, const FabArrayBase::MapOfCopyComTagContainers&>,
        IsDetected<PrepareRecvBuffers_t, DP&, FabArray<FAB>&, const FabArray<FAB>&, CommData&, const FabArrayBase::MapOfCopyComTagContainers&>,
        IsDetected<PackSendBuffers_t, DP&, const FabArray<FAB>&, CommData&>,
        IsDetected<UnpackRecvBuffers_t, DP&, FabArray<FAB>&, CommData&>
#endif
    > {};

////////////////////////////////////////////////////////////////////////////////////
//                                                      [DataPacking.PackComponents]
//
struct PackComponents {
    int dest_component{0};
    int src_component{0};
    int n_components{0};
};

template <typename FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
LocalCopy (const PackComponents& components, FabArray<FAB>& dest, const FabArray<FAB>& src,
           const FabArrayBase::CopyComTagsContainer& local_tags) {
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        local_copy_gpu(dest, src, components.dest_component, components.src_component,
                       components.n_components, local_tags);
    } else
#endif
    {
        local_copy_cpu(dest, src, components.dest_component, components.src_component,
                       components.n_components, local_tags);
    }
}

#ifdef AMREX_USE_MPI
void PrepareCommBuffers(CommData& comm, const PackComponents& components, 
                        const FabArrayBase::MapOfCopyComTagContainers& cctc,
                        std::size_t object_size, std::size_t align);

template <typename FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
PrepareSendBuffers (const PackComponents& components, FabArray<FAB>&, const FabArray<FAB>& src, CommData& comm,
                    const FabArrayBase::MapOfCopyComTagContainers& cctc) {
    using T = typename FAB::value_type;
    PrepareCommBuffers(comm, components, cctc, sizeof(T), alignof(T));
}

template <typename FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
PrepareRecvBuffers (const PackComponents& components, FabArray<FAB>&, const FabArray<FAB>&, CommData& comm,
                    const FabArrayBase::MapOfCopyComTagContainers& cctc) {
    using T = typename FAB::value_type;
    PrepareCommBuffers(comm, components, cctc, sizeof(T), alignof(T));
}

template <typename FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value> 
PackSendBuffers (const PackComponents& components, const FabArray<FAB>& src, CommData& send) {
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        FabArray<FAB>::pack_send_buffer_gpu(src, components.src_component, components.n_components,
                                            send.data, send.size, send.cctc);
    } else
#endif
    {
        FabArray<FAB>::pack_send_buffer_cpu(src, components.src_component, components.n_components,
                                            send.data, send.size, send.cctc);
    }
}

template <typename FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
UnpackRecvBuffers (const PackComponents& components, FabArray<FAB>& dest, const CommData& recv) {
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        unpack_recv_buffer_gpu(dest, components.dest_component, components.n_components, recv.data,
                               recv.size, recv.cctc);
    } else
#endif
    {
        unpack_recv_buffer_cpu(dest, components.dest_component, components.n_components, recv.data,
                               recv.size, recv.cctc);
    }
}
#endif

static_assert(IsDataPacking<PackComponents, FArrayBox>(),
              "PackComponents is expected to satisfy the concept DataPacking.");

////////////////////////////////////////////////////////////////////////////////////
//                                    [DataPacking.ApplyDtosAndProjectionOnReciever]
//
template <typename DTOS = Identity, typename FabProj = Identity>
struct ApplyDtosAndProjectionOnReciever : PackComponents {
    constexpr ApplyDtosAndProjectionOnReciever() = default;
    constexpr ApplyDtosAndProjectionOnReciever(const PackComponents& components, DTOS dtos_ = DTOS{}, FabProj proj_ = FabProj{})
        : PackComponents(components), dtos(std::move(dtos_)), proj(std::move(proj_)) {}

    AMREX_NO_UNIQUE_ADDRESS DTOS dtos;
    AMREX_NO_UNIQUE_ADDRESS FabProj proj;

    static_assert(IsCallableR<Dim3, DTOS, Dim3>(), "DTOS needs to be a callable: Dim3 -> Dim3");
    static_assert(IsFabProjection<FabProj, FArrayBox>(), "FabProj needs to be at least a projection on FArrayBox.");
};

template <typename FAB, typename DTOS, typename FabProj>
amrex::EnableIf_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<FabProj, FAB>()>
LocalCopy (const ApplyDtosAndProjectionOnReciever<DTOS, FabProj>& packing, FabArray<FAB>& dest,
           const FabArray<FAB>& src, const FabArrayBase::CopyComTagsContainer& local_tags) {
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        local_copy_gpu(dest, src, packing.dest_component, packing.src_component,
                       packing.n_components, local_tags, packing.dtos, packing.proj);
    } else
#endif
    {
        local_copy_cpu(dest, src, packing.dest_component, packing.src_component,
                       packing.n_components, local_tags, packing.dtos, packing.proj);
    }
}

#ifdef AMREX_USE_MPI
template <typename FAB, typename DTOS, typename FabProj>
amrex::EnableIf_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<FabProj, FAB>()>
UnpackRecvBuffers (const ApplyDtosAndProjectionOnReciever<DTOS, FabProj>& packing,
                   FabArray<FAB>& dest, const CommData& recv) {
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        unpack_recv_buffer_gpu(dest, packing.dest_component, packing.n_components, recv.data,
                               recv.size, recv.cctc, packing.dtos, packing.proj);
    } else
#endif
    {
        unpack_recv_buffer_cpu(dest, packing.dest_component, packing.n_components, recv.data,
                               recv.size, recv.cctc, packing.dtos, packing.proj);
    }
}
#endif

static_assert(IsDataPacking<ApplyDtosAndProjectionOnReciever<>, FArrayBox>(),
              "ApplyDtosAndProjectionOnReciever<> is expected to satisfy the DataPacking concept.");

/// Initiate recv and send calls for MPI and immediately return without doing any work.
///
/// DataPacking is a customization point object to control the behaviour of packing and unpacking
/// send or recv data buffers. It is used to perform interpolation or data transformations on either
/// sender or reciever side.
///
/// This function performs a data packing on sender side and we expect a call to Parallel_finish
/// that performs data unpacking on the reciever side.
///
/// \return Returns a CommHandler object that owns context and memory buffers for the whole life
/// time of the MPI transaction.
template <typename FAB, typename DataPacking, 
          typename = EnableIf_t<IsBaseFab<FAB>::value>,
          typename = EnableIf_t<IsDataPacking<DataPacking, FAB>::value>>
#ifdef AMREX_USE_MPI
AMREX_NODISCARD CommHandler 
ParallelCopy_nowait (FabArray<FAB>& dest, const FabArray<FAB>& src,
                     const FabArrayBase::CommMetaData& cmd, const DataPacking& data_packing) {
    CommHandler handler{};
    if (ParallelContext::NProcsSub() == 1) {
        return handler;
    }
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    handler.mpi_tag = ParallelDescriptor::SeqNum();

    const int N_rcvs = cmd.m_RcvTags ? cmd.m_RcvTags->size() : 0;
    if (N_rcvs > 0) {
        PrepareRecvBuffers(data_packing, dest, src, handler.recv, *cmd.m_RcvTags);
        PostRecvs(handler.recv, handler.mpi_tag);
    }

    const int N_snds = cmd.m_SndTags ? cmd.m_SndTags->size() : 0;
    if (N_snds > 0) {
        PrepareSendBuffers(data_packing, dest, src, handler.send, *cmd.m_SndTags);
        PackSendBuffers(data_packing, src, handler.send);
        PostSends(handler.send, handler.mpi_tag);
    }
    return handler;
}
#else
CommHandler ParallelCopy_nowait (FabArray<FAB>&, const FabArray<FAB>&,
                                 const FabArrayBase::CommMetaData&, const DataPacking&) {
    return CommHandler{};
}
#endif

static constexpr struct NoLocalCopy {} no_local_copy{};

template <typename FAB, typename DataPacking>
EnableIf_t<IsBaseFab<FAB>() && IsDataPacking<DataPacking, FAB>()>
#ifdef AMREX_USE_MPI
ParallelCopy_finish (NoLocalCopy, FabArray<FAB>& dest, const FabArray<FAB>& src, CommHandler handler,
                     const FabArrayBase::CommMetaData& cmd, const DataPacking& data_packing) {
    // If any FabArray is empty we have nothing to do.
    if (dest.size() == 0) {
        return;
    }
    // Return if nothing do
    if (ParallelContext::NProcsSub() == 1) {
        return;
    }
    // Unpack recieves
    const int N_rcvs = cmd.m_RcvTags ? cmd.m_RcvTags->size() : 0;
    if (N_rcvs > 0) {
        ParallelDescriptor::Waitall(handler.recv.request, handler.recv.stats);
#ifdef AMREX_D_DEBUG
        if (!CheckRcvStats(stats, handler.recv_size, handler.mpi_tag)) {
            amrex::Abort("NonLocalPC::ParallelCopy_finish failed with wrong message size");
        }
#endif
        UnpackRecvBuffers(data_packing, dest, handler.recv);
    }

    // Wait for all sends to be done
    const int N_snds = cmd.m_SndTags ? cmd.m_SndTags->size() : 0;
    if (N_snds > 0) {
        ParallelDescriptor::Waitall(handler.send.request, handler.send.stats);
    }
}
#else
ParallelCopy_finish (NoLocalCopy, FabArray<FAB>&, const FabArray<FAB>&, CommHandler,
                     const FabArrayBase::CommMetaData&, const DataPacking&) {}
#endif

template <typename FAB, typename DataPacking>
EnableIf_t<IsBaseFab<FAB>() && IsDataPacking<DataPacking, FAB>()>
ParallelCopy_finish (FabArray<FAB>& dest, const FabArray<FAB>& src, CommHandler handler,
                     const FabArrayBase::CommMetaData& cmd, const DataPacking& data_packing) {
    amrex::ignore_unused(handler);
    // If any FabArray is empty we have nothing to do.
    if (dest.size() == 0) {
        return;
    }
    // Eagerly do the local work and hope for some overlap with communication
    if (cmd.m_LocTags && cmd.m_LocTags->size() > 0) {
        LocalCopy(data_packing, dest, src, *cmd.m_LocTags);
    }
    ParallelCopy_finish(no_local_copy, dest, src, std::move(handler), cmd, data_packing);
}

template <typename T, typename... Args>
using InverseImage_t = decltype(std::declval<T>().InverseImage(std::declval<Args>()...));

template <typename T>
struct HasInverseImage : IsDetectedExact<Box, InverseImage_t, T&, Box> {};

//! This class stores data dependencies for an inter-block communication.
//!
//! In communication between two blocks one might need to do an index
//! transformation from one block to antoher.
struct MultiBlockCommMetaData : FabArrayBase::CommMetaData {
    enum { not_first_only = 0, first_only = 1 };

  //! \name Constructors

  //! @{
  //! Build global meta data that is being used to identify send and recv
  //! dependencies in communication routines.
  //!
  //! This call is quadratic in the number of boxes, i.e. it is in
  //! O(dstba.size() * srcba.size()). Therefore, it might be wise to cache the
  //! construction of this object to minimize its computation.
    template <typename DTOS, 
              typename = EnableIf_t<IsIndexMapping<DTOS>::value>,
              typename = EnableIf_t<HasInverseImage<DTOS>::value>>
    MultiBlockCommMetaData(const FabArrayBase& dst, const Box& dstbox, const FabArrayBase& src,
                           const IntVect& ngrow, DTOS dtos)
        : MultiBlockCommMetaData(dst.boxArray(), dst.DistributionMap(), dstbox, src.boxArray(),
                                 src.DistributionMap(), ngrow, dtos) {}

    template <typename DTOS, 
              typename = EnableIf_t<IsIndexMapping<DTOS>::value>,
              typename = EnableIf_t<HasInverseImage<DTOS>::value>>
    MultiBlockCommMetaData(const BoxArray& dstba, const DistributionMapping& dstdm,
                           const Box& dstbox, const BoxArray& srcba,
                           const DistributionMapping& srcdm, const IntVect& ngrow, DTOS dtos) {
        m_LocTags = std::make_unique<FabArrayBase::CopyComTagsContainer>();
        m_SndTags = std::make_unique<FabArrayBase::MapOfCopyComTagContainers>();
        m_RcvTags = std::make_unique<FabArrayBase::MapOfCopyComTagContainers>();
        const int myproc = ParallelDescriptor::MyProc();
        for (int i = 0, N = dstba.size(); i < N; ++i) {
            const int dest_owner = dstdm[i];
            const Box partial_dstbox = grow(dstba[i], ngrow) & dstbox;
            if (partial_dstbox.isEmpty()) {
                continue;
            }
            const Box partial_dstbox_mapped_in_src = dtos(partial_dstbox);
            std::vector<std::pair<int, Box>> boxes_from_src =
                srcba.intersections(partial_dstbox_mapped_in_src, not_first_only, ngrow);
            for (std::pair<int, Box> counted_box : boxes_from_src) {
                const int k = counted_box.first;
                const Box src_box = counted_box.second;
                AMREX_ASSERT(k < srcba.size());
                const int src_owner = srcdm[k];
                if (dest_owner == myproc || src_owner == myproc) {
                    if (src_owner == dest_owner) {
                        const BoxList tilelist(src_box, FabArrayBase::comm_tile_size);
                        for (const Box& tilebox : tilelist) {
                            const Box inverse_image = dtos.InverseImage(tilebox);
                            if ((inverse_image & partial_dstbox).ok()) {
                                m_LocTags->emplace_back(inverse_image, tilebox, i, k);
                            }
                        }
                    } else {
                        const Box inverse_image = dtos.InverseImage(src_box);
                        if ((inverse_image & partial_dstbox).ok()) {
                            FabArrayBase::CopyComTagsContainer& copy_tags =
                                (src_owner == myproc) ? (*m_SndTags)[dest_owner]
                                                    : (*m_RcvTags)[src_owner];
                            copy_tags.emplace_back(inverse_image, src_box, i, k);
                        }
                    }
                }
            }
        }
  }
  //! @}
};

template <typename FAB, typename DTOS = Identity, typename Proj = Identity>
EnableIf_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
ParallelCopy (FabArray<FAB>& dest, const FabArray<FAB>& src, const FabArrayBase::CommMetaData& cmd,
              int destcomp, int srccomp, int numcomp, DTOS dtos = DTOS{}, Proj proj = Proj{}) {
    BL_PROFILE("NonLocalBC::ParallelCopy_with_MultiBlockCommMetaData");
    ApplyDtosAndProjectionOnReciever<DTOS, Proj> packing{PackComponents{destcomp, srccomp, numcomp},
                                                         std::move(dtos), std::move(proj)};
    CommHandler handler = ParallelCopy_nowait(dest, src, cmd, packing);
    ParallelCopy_finish(dest, src, std::move(handler), cmd, packing);
}

template <typename FAB, typename DTOS = Identity, typename Proj = Identity>
EnableIf_t<IsBaseFab<FAB>() && IsIndexMapping<DTOS>() && HasInverseImage<DTOS>() && IsFabProjection<Proj, FAB>()>
ParallelCopy (FabArray<FAB>& dest, const Box& destbox, const FabArray<FAB>& src, int destcomp,
              int srccomp, int numcomp, const IntVect& ngrow, DTOS dtos = DTOS{}, Proj proj = Proj{}) {
    BL_PROFILE("NonLocalBC::ParallelCopy");
    MultiBlockCommMetaData cmd(dest, destbox, src, ngrow, dtos);
    ParallelCopy(dest, src, cmd, destcomp, srccomp, numcomp, std::move(dtos), std::move(proj));
}

template <class FAB, class DTOS>
#ifdef AMREX_USE_MPI
AMREX_NODISCARD
#endif
CommHandler
Comm_nowait (FabArray<FAB>& mf, int scomp, int ncomp, FabArrayBase::CommMetaData const& cmd,
             DTOS dtos)
{
#ifdef AMREX_USE_MPI
    if (ParallelContext::NProcsSub() == 1)
#endif
    {
        int N_locs = (*cmd.m_LocTags).size();
        if (N_locs == 0) return CommHandler{};
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            local_copy_gpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos);
        } else
#endif
        {
            local_copy_cpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos);
        }
        return CommHandler{};
    }

#ifdef AMREX_USE_MPI
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum = ParallelDescriptor::SeqNum();

    const int N_locs = cmd.m_LocTags->size();
    const int N_rcvs = cmd.m_RcvTags->size();
    const int N_snds = cmd.m_SndTags->size();

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0) {
        // No work to do.
        return CommHandler{};
    }

    CommHandler handler{};
    handler.mpi_tag = SeqNum;

    if (N_rcvs > 0) {
        handler.recv.the_data = mf.PostRcvs(*cmd.m_RcvTags, handler.recv.data, handler.recv.size,
                                            handler.recv.rank, handler.recv.request, ncomp, SeqNum);
    }

    if (N_snds > 0) {
        handler.send.the_data =
            mf.PrepareSendBuffers(*cmd.m_SndTags, handler.send.data, handler.send.size,
                                  handler.send.rank, handler.send.request, handler.send.cctc, ncomp);

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            FabArray<FAB>::pack_send_buffer_gpu(mf, scomp, ncomp, handler.send.data,
                                                handler.send.size, handler.send.cctc);
        } else
#endif
        {
            FabArray<FAB>::pack_send_buffer_cpu(mf, scomp, ncomp, handler.send.data,
                                                handler.send.size, handler.send.cctc);
        }

        FabArray<FAB>::PostSnds(handler.send.data, handler.send.size, handler.send.rank, handler.send.request, SeqNum);
    }

    if (N_locs > 0)
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            local_copy_gpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos);
        } else
#endif
        {
            local_copy_cpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos);
        }
    }

    return handler;
#endif
}

#ifdef AMREX_USE_MPI
template <class FAB, class DTOS>
void
Comm_finish (FabArray<FAB>& mf, int scomp, int ncomp, FabArrayBase::CommMetaData const& cmd,
             CommHandler handler, DTOS dtos)
{
    if (ParallelContext::NProcsSub() == 1) return;

    const int N_rcvs = cmd.m_RcvTags->size();
    if (N_rcvs > 0)
    {
        handler.recv.cctc.resize(N_rcvs, nullptr);
        for (int k = 0; k < N_rcvs; ++k) {
            auto const& cctc = cmd.m_RcvTags->at(handler.recv.rank[k]);
            handler.recv.cctc[k] = &cctc;
        }
        handler.recv.stats.resize(handler.recv.request.size());
        ParallelDescriptor::Waitall(handler.recv.request, handler.recv.stats);
#ifdef AMREX_D_DEBUG
        if (!CheckRcvStats(stats, handler.recv.size, handler.mpi_tag)) {
            amrex::Abort("NonLocalBC::Comm_finish failed with wrong message size");
        }
#endif

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            unpack_recv_buffer_gpu(mf, scomp, ncomp, handler.recv.data,
                                   handler.recv.size, handler.recv.cctc, dtos);
        } else
#endif
        {
            unpack_recv_buffer_cpu(mf, scomp, ncomp, handler.recv.data,
                                   handler.recv.size, handler.recv.cctc, dtos);
        }
    }

    const int N_snds = cmd.m_SndTags->size();
    if (N_snds > 0) {
        handler.send.stats.resize(handler.send.request.size());
        ParallelDescriptor::Waitall(handler.send.request, handler.send.stats);
    }
}
#endif

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate90 (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain)
{
    BL_PROFILE("Rotate90");

    AMREX_ASSERT(domain.cellCentered());
    AMREX_ASSERT(domain.smallEnd() == 0);
    AMREX_ASSERT(domain.length(0) == domain.length(1));
    AMREX_ASSERT(mf.is_cell_centered());
    AMREX_ASSERT(scomp < mf.nComp() && scomp+ncomp <= mf.nComp());
    AMREX_ASSERT(nghost.allLE(mf.nGrowVect()) && nghost[0] == nghost[1]);

    if (nghost[0] <= 0) return;

    const FabArrayBase::RB90& TheRB90 = mf.getRB90(nghost, domain);

    auto handler = Comm_nowait(mf, scomp, ncomp, TheRB90,Rotate90DstToSrc{});

    Box corner(-nghost, IntVect{AMREX_D_DECL(-1,-1,domain.bigEnd(2)+nghost[2])});
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        Box const& bx = corner & mfi.fabbox();
        if (bx.ok()) {
            auto const& fab = mf.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx,ncomp,i,j,k,n,
            {
                fab(i,j,k,n) = fab(-i-1,-j-1,k,n);
            });
        }
    }

#ifdef AMREX_USE_MPI
    Comm_finish(mf, scomp, ncomp, TheRB90, std::move(handler), Rotate90DstToSrc{});
#else
    amrex::ignore_unused(handler);
#endif
}

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate90 (FabArray<FAB>& mf, Box const& domain)
{
    Rotate90(mf, 0, mf.nComp(), mf.nGrowVect(), domain);
}

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate180 (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain)
{
    BL_PROFILE("Rotate180");

    AMREX_ASSERT(domain.cellCentered());
    AMREX_ASSERT(domain.smallEnd() == 0);
    AMREX_ASSERT(domain.length(1) % 2 == 0);
    AMREX_ASSERT(mf.is_cell_centered());
    AMREX_ASSERT(scomp < mf.nComp() && scomp+ncomp <= mf.nComp());
    AMREX_ASSERT(nghost.allLE(mf.nGrowVect()));

    if (nghost[0] <= 0) return;

    const FabArrayBase::RB180& TheRB180 = mf.getRB180(nghost, domain);

    auto handler = Comm_nowait(mf, scomp, ncomp, TheRB180, Rotate180Fn{domain.length(1)});

#ifdef AMREX_USE_MPI
    Comm_finish(mf, scomp, ncomp, TheRB180, std::move(handler), Rotate180Fn{domain.length(1)});
#else
    amrex::ignore_unused(handler);
#endif
}

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate180 (FabArray<FAB>& mf, Box const& domain)
{
    Rotate180(mf, 0, mf.nComp(), mf.nGrowVect(), domain);
}

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
FillPolar (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain)
{
    BL_PROFILE("FillPolar");

    AMREX_ASSERT(domain.cellCentered());
    AMREX_ASSERT(domain.smallEnd() == 0);
    AMREX_ASSERT(domain.length(1) % 2 == 0);
    AMREX_ASSERT(mf.is_cell_centered());
    AMREX_ASSERT(scomp < mf.nComp() && scomp+ncomp <= mf.nComp());
    AMREX_ASSERT(nghost.allLE(mf.nGrowVect()));

    if (nghost[0] <= 0) return;

    const FabArrayBase::PolarB& ThePolarB = mf.getPolarB(nghost, domain);

    auto handler = Comm_nowait(mf, scomp, ncomp, ThePolarB,
                               PolarFn{domain.length(0), domain.length(1)});

#ifdef AMREX_USE_MPI
    Comm_finish(mf, scomp, ncomp, ThePolarB, std::move(handler),
                PolarFn{domain.length(0), domain.length(1)});
#else
    amrex::ignore_unused(handler);
#endif
}

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
FillPolar (FabArray<FAB>& mf, Box const& domain)
{
    FillPolar(mf, 0, mf.nComp(), mf.nGrowVect(), domain);
}

extern template
void Rotate90(FabArray<FArrayBox>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain);

extern template
void Rotate90(FabArray<FArrayBox>& mf, Box const& domain);

extern template
void Rotate180(FabArray<FArrayBox>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain);

extern template 
void Rotate180(FabArray<FArrayBox>& mf, Box const& domain);

extern template
void FillPolar(FabArray<FArrayBox>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain);

extern template
void FillPolar (FabArray<FArrayBox>& mf, Box const& domain);

extern template
void ParallelCopy(FabArray<FArrayBox>& dest, const Box& destbox, const FabArray<FArrayBox>& src, int destcomp,
             int srccomp, int numcomp, const IntVect& ngrow, Identity, Identity);

extern template
void ParallelCopy(FabArray<FArrayBox>& dest, const Box& destbox, const FabArray<FArrayBox>& src, int destcomp,
             int srccomp, int numcomp, const IntVect& ngrow, MultiBlockIndexMapping, Identity);

}}


#if AMREX_HAS_STATIC_ASSERT_NO_MESSAGE == 0
#  undef AMREX_STATIC_ASSERT_NO_MESSAGE
#endif

#endif
