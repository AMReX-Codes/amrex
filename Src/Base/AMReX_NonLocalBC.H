#ifndef AMREX_NONLOCAL_BC_H_
#define AMREX_NONLOCAL_BC_H_
#include <AMReX_Config.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_FabArray.H>
#include <AMReX_FArrayBox.H>

namespace amrex::NonLocalBC {

////////////////////////////////////////////////////////////////////////////////////
//                                                            [concept.IndexMapping]
//

//! \brief Return type of an InverseImage class member function.
template <typename T, typename... Args>
using Inverse_t = decltype(std::declval<T>().Inverse(std::declval<Args>()...));

//! \brief Type trait that tests if T has an InverseImage class member function.
template <typename T>
struct HasInverseMemFn : IsDetectedExact<Dim3, Inverse_t, T, Dim3> {};

//! \brief Tests if a given type IndexMap is usable as an index mapping between two index based
//! coordinate systems.
template <class IndexMap>
struct IsIndexMapping
    : Conjunction<IsCallableR<Dim3, IndexMap&, const Dim3&>,
                  HasInverseMemFn<const IndexMap&>> {};

////////////////////////////////////////////////////////////////////////////////////
//                                             [IndexMapping.MultiBlockIndexMapping]
//

//! \brief This struct describes an affine index transformation for two coordinate systems.
//!
//! This DTOS is used to map indices from one block to another. It respects a permutation of
//! (x,y,z) coordinates, an offset and a change in orientation.
struct MultiBlockIndexMapping {
    //! \brief This vector needs to be a valid permutation
    IntVect permutation{AMREX_D_DECL(0, 1, 2)};
    //! \brief The offset in the source index space.
    IntVect offset{AMREX_D_DECL(0, 0, 0)};
    //! \brief A vector of 1 and -1 describing the orientation in each component.
    IntVect sign{AMREX_D_DECL(1, 1, 1)};

    //! \brief Applies this mapping on the index from destination space and returns an index in the
    //!        source space.
    //!
    //! \param[in] i  The index that lives in the destination space.
    //!
    //! \return Returns an index in the source space that is given by this affine transformation.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE Dim3 operator()(Dim3 i) const noexcept {
        int iv[3]{i.x, i.y, i.z};
        int iv_new[3]{};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            iv_new[d] = sign[d] * (iv[permutation[d]] - offset[d]);
        }
        return {iv_new[0], iv_new[1], iv_new[2]};
    }

    //! \brief The inverse function is given by rearringing all above terms.
    //!
    //! \param[in] i  The index that lives in the source space.
    //!
    //! \return Returns an index in the destination space that is given by this affine
    //! transformation.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE Dim3 Inverse(Dim3 i) const noexcept {
        int iv_new[3]{i.x, i.y, i.z};
        int iv[3]{};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            AMREX_ASSERT(sign[d] == 1 || sign[d] == -1);
            iv[permutation[d]] = iv_new[d] * sign[d] + offset[d];
        }
        return {iv[0], iv[1], iv[2]};
    }

    [[nodiscard]] IndexType operator()(IndexType it) const noexcept {
        return IndexType{IntVect{AMREX_D_DECL(it[permutation[0]], it[permutation[1]], it[permutation[2]])}};
    }

    [[nodiscard]] IndexType Inverse(IndexType it) const noexcept {
        IntVect inverse_permutation;
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            inverse_permutation[permutation[i]] = i;
        }
        return IndexType{IntVect{AMREX_D_DECL(it[inverse_permutation[0]], it[inverse_permutation[1]], it[inverse_permutation[2]])}};
    }
};

//! \brief Applies the Dim3 to Dim3 mapping onto IntVects.
//!
//! This is used to map indices from the dest index space into the source index space.
//!
//! IntVect is being embedded in Dim3 by trailing zeros.
//! Dim3 is being projected to IntVect by dopping z for AMREX_SPACEDIM = 2 or z and y components
//! for AMREX_SPACEDIM = 1.
//!
//! \param[in] iv The IntVect that lives in the destination index space.
//!
//! \return Returns IntVect{dtos(Dim3{iv})}
template <typename DTOS>
std::enable_if_t<IsCallableR<Dim3, DTOS, Dim3>::value, IntVect>
Apply(DTOS dtos, const IntVect& iv)
{
    Dim3 i = dtos(iv.dim3());
    return IntVect{AMREX_D_DECL(i.x, i.y, i.z)};
}

//! \brief Applies the Dim3 to Dim3 mapping onto Boxes but does not change the index type.
//!
//! This function assumes monotonicity of dtos in each component.
//!
//! \param[in] box The box that lives the destination index space
//!
//! \return Returns the smallest Box in the source index space that contains images of
//! Apply(dtos, box.smallEnd()) and Apply(dtos, box.bigEnd()).
template <typename DTOS>
std::enable_if_t<IsCallableR<Dim3, DTOS, Dim3>::value && !IsCallableR<IndexType, DTOS, IndexType>::value, Box>
Image(DTOS dtos, const Box& box)
{
    IntVect mapped_smallEnd = Apply(dtos, box.smallEnd());
    IntVect mapped_bigEnd = Apply(dtos, box.bigEnd());
    IntVect smallEnd;
    IntVect bigEnd;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        smallEnd[d] = std::min(mapped_smallEnd[d], mapped_bigEnd[d]);
        bigEnd[d] = std::max(mapped_smallEnd[d], mapped_bigEnd[d]);
    }
    return Box{smallEnd, bigEnd, box.ixType()};
}

//! \brief Applies the Dim3 to Dim3 mapping onto Boxes and maps the index type.
//!
//! This function assumes monotonicity of dtos in each component.
//!
//! \param[in] box The box that lives the destination index space
//!
//! \return Returns the smallest Box in the source index space that contains images of
//! Apply(dtos, box.smallEnd()) and Apply(dtos, box.bigEnd()).
template <typename DTOS>
std::enable_if_t<IsCallableR<Dim3, DTOS, Dim3>::value && IsCallableR<IndexType, DTOS, IndexType>::value, Box>
Image(DTOS dtos, const Box& box)
{
    // "Forget" the index type mapping and invoke Image without changing the index type.
    Box srcbox = Image([&dtos](Dim3 d) { return dtos(d); }, box);
    // Fix the index type of the resulting box
    srcbox.setType(dtos(box.ixType()));
    return srcbox;
}

//! \brief Applies the Dim3 to Dim3 invserse mapping onto IntVects.
//!
//! This is used to map indices from the src index space into the dest index space.
//!
//! IntVect is being embedded in Dim3 by trailing zeros.
//! Dim3 is being projected to IntVect by dopping z for AMREX_SPACEDIM = 2 or z and y components
//! for AMREX_SPACEDIM = 1.
//!
//! \param[in] iv The IntVect that lives in the src index space.
//!
//! \return Returns IntVect{dtos.Inverse(Dim3{iv})}
template <typename DTOS>
std::enable_if_t<HasInverseMemFn<DTOS>::value, IntVect>
ApplyInverse(DTOS dtos, const IntVect& iv)
{
    return Apply([&dtos](Dim3 i) { return dtos.Inverse(i); }, iv);
}

//! \brief Applies the inverse Dim3 to Dim3 mapping onto Boxes without changing the index type.
//!
//! This function assumes monotonicity in each component.
//!
//! \param[in] box The box that lives the source index space
//!
//! \return Returns the smallest Box in the destination index space that contains images of
//! ApplyInverse(box.smallEnd()) and ApplyInverse(box.bigEnd()).
template <typename DTOS>
std::enable_if_t<HasInverseMemFn<DTOS>::value && !IsCallableR<IndexType, DTOS, IndexType>::value, Box>
InverseImage(DTOS dtos, const Box& box)
{
    return Image([&dtos](Dim3 i) { return dtos.Inverse(i); }, box);
}

//! \brief Applies the inverse Dim3 to Dim3 mapping onto Boxes
//!
//! This function assumes monotonicity in each component.
//!
//! \param[in] box The box that lives the source index space
//!
//! \return Returns the smallest Box in the destination index space that contains images of
//! ApplyInverse(box.smallEnd()) and ApplyInverse(box.bigEnd()).
template <typename DTOS>
std::enable_if_t<HasInverseMemFn<DTOS>::value && IsCallableR<IndexType, DTOS, IndexType>::value, Box>
InverseImage(DTOS dtos, const Box& box)
{
    return Image([&dtos](auto&& i) { return dtos.Inverse(i); }, box);
}

//! \brief This is the index mapping based on the DTOS MultiBlockDestToSrc.
static_assert(IsIndexMapping<MultiBlockIndexMapping>(), // NOLINT(bugprone-throw-keyword-missing)
              "MultiBlockIndexMapping is expected to satisfy IndexMapping");

////////////////////////////////////////////////////////////////////////////////////
//                                                    [class.MultiBlockCommMetaData]
//

//! \brief This class stores data dependencies for an inter-block communication.
//!
//! In communication between two blocks one might need to do an index
//! transformation from one block to antoher.
struct MultiBlockCommMetaData : FabArrayBase::CommMetaData {
    //! \name Constructors

    //! \brief Build global meta data by calling the define() member function.
    //!
    //! \see MultiBlockCommMetaData::define
    template <typename DTOS,
              typename = std::enable_if_t<IsIndexMapping<DTOS>::value>>
    MultiBlockCommMetaData(const FabArrayBase& dst, const Box& dstbox, const FabArrayBase& src,
                           const IntVect& ngrow, DTOS dtos);

    //! \brief Build global meta data by calling the define() member function.
    //!
    //! \see MultiBlockCommMetaData::define
    template <typename DTOS,
              typename = std::enable_if_t<IsIndexMapping<DTOS>::value>>
    MultiBlockCommMetaData(const BoxArray& dstba, const DistributionMapping& dstdm,
                           const Box& dstbox, const BoxArray& srcba,
                           const DistributionMapping& srcdm, const IntVect& ngrow, DTOS dtos);

    //! \name Manipulators

    //! \brief Build global meta data that is being used to identify send and recv
    //! dependencies in communication routines.
    //!
    //! This call is quadratic in the number of boxes, i.e. it is in
    //! O(dstba.size() * srcba.size()). Therefore, it might be wise to cache the
    //! construction of this object to minimize its computation.
    //!
    //! \tparam DTOS This parameter needs to satisfy IsIndexMapping and HasInverseImage.
    //!
    //! \param[in] dstba  The destination box array.
    //!
    //! \param[in] dstdm  The destination distribution mapping.
    //!
    //! \param[in] dstbox The box that will be filled and that lives in the destination space.
    //!
    //! \param[in] srcba  The source box array.
    //!
    //! \param[in] srcdm  The source distribution mapping.
    //!
    //! \param[in] ngrow  The number of ghost cells that shall be considered when filling data.
    //!
    //! \param[in] dtos   The dest to source index mapping that has an inverse image.
    //!
    //! \return Nothing.
    template <typename DTOS>
    std::enable_if_t<IsIndexMapping<DTOS>::value>
    define(const BoxArray& dstba, const DistributionMapping& dstdm, const Box& dstbox,
           const BoxArray& srcba, const DistributionMapping& srcdm, const IntVect& ngrow,
           DTOS dtos);
};

////////////////////////////////////////////////////////////////////////////////////
//                                                           [concept.FabProjection]
//

//! \brief This type trait tests if a type P is a projetion for FAB.
template <typename P, typename FAB>
struct IsFabProjection
    : IsCallableR<typename FAB::value_type, P, Array4<const typename FAB::value_type>, Dim3, int>
{};

////////////////////////////////////////////////////////////////////////////////////
//                                                          [FabProjection.Identity]
//                                                           [IndexMapping.Identity]

//! \brief This class acts as a default no-op operator.
//!
//! This class satisfies IndexMapping and FabProjection.
struct Identity  {
    //! The identity function for Dim3
    constexpr Dim3 operator()(Dim3 i) const noexcept { return i; }
    //! The identity function for Dim3
    [[nodiscard]] static constexpr Dim3 Inverse(Dim3 i) noexcept { return i; }

    //! \return Returns array(i.x,i.y,i.z,comp)
    template <typename T>
    constexpr T operator()(Array4<const T> array, Dim3 i, int comp = 0) const
        noexcept(noexcept(array(i.x, i.y, i.z, comp))) {
        return array(i.x, i.y, i.z, comp);
    }

    //! The identity for int
    constexpr int operator()(int i) const noexcept { return i; }
};
static constexpr Identity identity{};

static_assert(sizeof(Identity) == 1 );
static_assert(std::is_trivially_default_constructible<Identity>::value );
static_assert(std::is_trivially_copy_assignable<Identity>::value );
static_assert(std::is_trivially_copy_constructible<Identity>::value );
static_assert(IsIndexMapping<Identity>() ); // NOLINT(bugprone-throw-keyword-missing)
static_assert(IsFabProjection<Identity, FArrayBox>() ); // NOLINT(bugprone-throw-keyword-missing)

////////////////////////////////////////////////////////////////////////////////////
//                                                     [FabProjection.MapComponents]

//! \brief This class takes a projection and a component map and combines them to form a new
//! projection.
//!
//! We use this to apply a permutations on components FABs.
//!
//! \see SwapComponents
template <typename Base, typename Map = Identity> struct MapComponents {
    static_assert(IsCallable<Base, Array4<const Real>, Dim3, int>::value,
                  "Base needs to be a callable function: (Array4<const T>, Dim3, i) -> auto.");

    static_assert(IsCallableR<int, Map, int>::value,
                  "Map needs to be a callable function: int -> int.");

    Base base;
    Map map;

    template <typename T,
              typename = std::enable_if_t<IsCallable<Base, Array4<const T>, Dim3, int>::value>,
              typename = std::enable_if_t<IsCallableR<int, Map, int>::value>>
    constexpr decltype(auto) operator()(Array4<const T> array, Dim3 i, int comp) const
        noexcept(noexcept(base(array, i, map(comp)))) {
        return base(array, i, map(comp));
    }
};

static_assert(std::is_trivially_copy_assignable<MapComponents<Identity>>() ); // NOLINT(bugprone-throw-keyword-missing)
static_assert(std::is_trivially_copy_constructible<MapComponents<Identity>>() ); // NOLINT(bugprone-throw-keyword-missing)
static_assert(IsFabProjection<MapComponents<Identity>, FArrayBox>() ); // NOLINT(bugprone-throw-keyword-missing)

////////////////////////////////////////////////////////////////////////////////////
//                                      [FabProjection.MapComponents.SwapComponents]

//! \brief This is a permutation where only two components are swapped.
//!
//! This class is used with MapComponents to define fab projections that swap components on FABs,
//! e.g. for velocity components. The template arguments I, J should be greater or equal to -1. The
//! value -1 indicates run-time values and increase the object size of this function object.
template <int I, int J> struct SwapComponents {
    static_assert(I >= 0 && J >= 0, "I >= 0 && J >= 0");

    //! \brief Swaps indices I and J.
    //!
    //! \return If i == I then return J. If i == J then return I. Otherwise returns i.
    constexpr int operator()(int i) const noexcept {
        const int map[2] = {I, J};
        return i == I || i == J ? map[std::size_t(i == I)] : i;
    }
};

template <int I> struct SwapComponents<I, -1> {
    static_assert(I >= 0, "I >= 0");

    int J;
    constexpr int operator()(int i) const noexcept {
        const int map[2] = {I, J};
        return i == I || i == J ? map[std::size_t(i == I)] : i;
    }
};

template <int J> struct SwapComponents<-1, J> {
    static_assert(J >= 0, "J >= 0");

    int I;
    constexpr int operator()(int i) const noexcept {
        const int map[2] = {I, J};
        return i == I || i == J ? map[std::size_t(i == I)] : i;
    }
};

template <> struct SwapComponents<-1, -1> {
    int I;
    int J;
    constexpr int operator()(int i) const noexcept {
        const int map[2] = {I, J};
        return i == I || i == J ? map[std::size_t(i == I)] : i;
    }
};

using DynamicSwapComponents = SwapComponents<-1, -1>;

template <int I, int J> static constexpr SwapComponents<I, J> swap_indices{};

static_assert(sizeof(SwapComponents<0, 1>) == 1 );
static_assert(sizeof(DynamicSwapComponents) == 2 * sizeof(int) );
static_assert(sizeof(SwapComponents<0, -1>) == sizeof(int) );
static_assert(sizeof(SwapComponents<-1, 1>) == sizeof(int) );
static_assert(std::is_trivially_default_constructible<MapComponents<Identity, SwapComponents<0, 1>>>() ); // NOLINT(bugprone-throw-keyword-missing)
static_assert(std::is_trivially_copy_assignable<MapComponents<Identity, SwapComponents<0, 1>>>() ); // NOLINT(bugprone-throw-keyword-missing)
static_assert(std::is_trivially_copy_constructible<MapComponents<Identity, SwapComponents<0, 1>>>() ); // NOLINT(bugprone-throw-keyword-missing)
static_assert(IsFabProjection<MapComponents<Identity, SwapComponents<0, 1>>, FArrayBox>() ); // NOLINT(bugprone-throw-keyword-missing)

static_assert(swap_indices<0, 1>(0) == 1 );
static_assert(swap_indices<0, 1>(1) == 0 );
static_assert(swap_indices<0, 1>(2) == 2 );
static_assert(DynamicSwapComponents{0, 1}(0) == 1 );
static_assert(DynamicSwapComponents{0, 1}(1) == 0 );
static_assert(DynamicSwapComponents{0, 1}(2) == 2 );

////////////////////////////////////////////////////////////////////////////////////
//                                                                  [class.CommData]
//                                                               [class.CommHandler]

//! \brief This class holds data buffers for either immediate MPI send or recv calls.
//!
//! The data buffers need to be kept alive until all MPI_Request have been completed.
//! All member variables are resized a corresponding call to PrepareCommBuffers and
//! the data transactions are initiated either by calls to PostRecvs or PostSends.
struct CommData {
#ifdef AMREX_USE_MPI
    //! \brief Holds 'em all in one data pointer.
    TheFaArenaPointer the_data = nullptr;
    //! \brief Stores MPI ranks. For recvs it is the 'from' rank and for sends it is the 'to' rank.
    Vector<int> rank{};
    //! \brief Pointers to the_data that can be used for each single data transaction.
    Vector<char*> data{};
    //! \brief All offsets of data in the_data.
    Vector<std::size_t> offset{};
    //! \brief The size in bytes for each data transaction.
    Vector<std::size_t> size{};
    //! \brief The associated MPI_Request for each data transaction.
    Vector<MPI_Request> request{};
    //! \brief For each request the corresponding MPI_status, used for debugging.
    Vector<MPI_Status> stats{};
    //! \brief For each request the copy comm tags for the corresponding data FABs.
    Vector<const FabArrayBase::CopyComTagsContainer*> cctc{};
#endif
};

#ifdef AMREX_USE_MPI
//! \brief Fill all class member variables of comm but the request and the stats vector.
void PrepareCommBuffers(CommData& comm, const FabArrayBase::MapOfCopyComTagContainers& cctc,
                        int n_components, std::size_t object_size, std::size_t align);

//! \brief Initiate all recvieves with MPI_Irecv calls associated with tag mpi_tag.
void PostRecvs(CommData& recv, int mpi_tag);

//! \brief Initiate all sends with MPI_Isend calls associated with tag mpi_tag.
void PostSends(CommData& send, int mpi_tag);
#endif


//! \brief This class stores both recv and send buffers with an associated MPI tag.
//!
//! This handler stores data for a whole ParallelCopy MPI transaction.
struct CommHandler {
#ifdef AMREX_USE_MPI
    int mpi_tag{};
    CommData recv{};
    CommData send{};
#endif
};

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

//! \brief This type trait tests if a given type DP satisfies the DataPacking concept for type FAB.
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

template <class FAB, class DTOS = Identity, class Proj = Identity>
std::enable_if_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
local_copy_cpu (FabArray<FAB>& dest, const FabArray<FAB>& src, int dcomp, int scomp, int ncomp,
                FabArrayBase::CopyComTagsContainer const& local_tags, DTOS dtos = DTOS{},
                Proj proj = Proj{}) noexcept;

template <class FAB, class DTOS = Identity, class Proj = Identity>
std::enable_if_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
unpack_recv_buffer_cpu (FabArray<FAB>& mf, int dcomp, int ncomp, Vector<char*> const& recv_data,
                        Vector<std::size_t> const& recv_size,
                        Vector<FabArrayBase::CopyComTagsContainer const*> const& recv_cctc,
                        DTOS dtos = DTOS{}, Proj proj = Proj{}) noexcept;

#ifdef AMREX_USE_GPU
template <class FAB, class DTOS = Identity, class Proj = Identity>
std::enable_if_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
local_copy_gpu (FabArray<FAB>& dest, const FabArray<FAB>& src, int dcomp, int scomp, int ncomp,
                FabArrayBase::CopyComTagsContainer const& local_tags, DTOS dtos = DTOS{},
                Proj proj = Proj{}) noexcept;

template <class FAB, class DTOS = Identity, class Proj = Identity>
std::enable_if_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
unpack_recv_buffer_gpu (FabArray<FAB>& mf, int scomp, int ncomp,
                        Vector<char*> const& recv_data,
                        Vector<std::size_t> const& recv_size,
                        Vector<FabArrayBase::CopyComTagsContainer const*> const& recv_cctc,
                        DTOS dtos = DTOS{}, Proj proj = Proj{});
#endif

////////////////////////////////////////////////////////////////////////////////////
//                                                      [DataPacking.PackComponents]
//
// PackComponents is the simplest data packing policy.
// This provides us with sane default behaviour that we can use when defining new
// data packing policies.

//! \brief Contains information about which components take part of the data transaction.
struct PackComponents {
    int dest_component{0};
    int src_component{0};
    int n_components{0};
};

//! \brief Dispatch local copies to the default behaviour that knows no DTOS nor projection.
template <typename FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
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
//! \brief Calls PrepareComBuffers.
template <typename FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
PrepareSendBuffers (const PackComponents& components, FabArray<FAB>& dest, const FabArray<FAB>& src,
                    CommData& comm, const FabArrayBase::MapOfCopyComTagContainers& cctc) {
    using T = typename FAB::value_type;
    ignore_unused(dest, src);
    PrepareCommBuffers(comm, cctc, components.n_components, sizeof(T), alignof(T));
}

//! \brief Calls PrepareComBuffers.
template <typename FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
PrepareRecvBuffers (const PackComponents& components, FabArray<FAB>& dest, const FabArray<FAB>& src,
                    CommData& comm, const FabArrayBase::MapOfCopyComTagContainers& cctc) {
    using T = typename FAB::value_type;
    ignore_unused(dest, src);
    PrepareCommBuffers(comm, cctc, components.n_components, sizeof(T), alignof(T));
}

//! \brief Serializes FAB data without any knowledge of a DTOS nor a projection.
template <typename FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
PackSendBuffers (const PackComponents& components, const FabArray<FAB>& src, CommData& send) {
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        FabArray<FAB>::pack_send_buffer_gpu(src, components.src_component, components.n_components,
                                            send.data, send.size, send.cctc);
    } else
#endif // AMREX_USE_GPU
    {
        FabArray<FAB>::pack_send_buffer_cpu(src, components.src_component, components.n_components,
                                            send.data, send.size, send.cctc);
    }
}

//! \brief De-serializes FAB data without any knowledge of a DTOS nor a projection.
template <typename FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
UnpackRecvBuffers (const PackComponents& components, FabArray<FAB>& dest, const CommData& recv) {
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        unpack_recv_buffer_gpu(dest, components.dest_component, components.n_components, recv.data,
                               recv.size, recv.cctc);
    } else
#endif // AMREX_USE_GPU
    {
        unpack_recv_buffer_cpu(dest, components.dest_component, components.n_components, recv.data,
                               recv.size, recv.cctc);
    }
}
#endif // AMREX_USE_MPI

static_assert(IsDataPacking<PackComponents, FArrayBox>(), // NOLINT(bugprone-throw-keyword-missing)
              "PackComponents is expected to satisfy the concept DataPacking.");

////////////////////////////////////////////////////////////////////////////////////
//                                    [DataPacking.ApplyDtosAndProjectionOnReciever]
//
//! \brief This class specializes behaviour on local copies and unpacking receive buffers.
//!
//! It takes a DTOS and FabProjection to apply them on the receiver side.
template <typename DTOS = Identity, typename FabProj = Identity>
struct ApplyDtosAndProjectionOnReciever : PackComponents {
    constexpr ApplyDtosAndProjectionOnReciever() = default;
    constexpr ApplyDtosAndProjectionOnReciever(const PackComponents& components, DTOS dtos_ = DTOS{}, FabProj proj_ = FabProj{})
        : PackComponents(components), dtos(std::move(dtos_)), proj(std::move(proj_)) {}

    DTOS dtos;
    FabProj proj;

    static_assert(IsCallableR<Dim3, DTOS, Dim3>(), "DTOS needs to be a callable: Dim3 -> Dim3");
    static_assert(IsFabProjection<FabProj, FArrayBox>(), "FabProj needs to be at least a projection on FArrayBox.");
};

//! \brief Do local copies of FABs using DTOS and projection.
template <typename FAB, typename DTOS, typename FabProj>
std::enable_if_t<IsBaseFab<FAB>::value>
LocalCopy (const ApplyDtosAndProjectionOnReciever<DTOS, FabProj>& packing, FabArray<FAB>& dest,
           const FabArray<FAB>& src, const FabArrayBase::CopyComTagsContainer& local_tags) {
    static_assert(IsFabProjection<FabProj, FAB>(), "FabProj needs to be a projection for given FAB type.");
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
//! \brief Copy from received data in the buffer to destination FABs using DTOS and projection.
template <typename FAB, typename DTOS, typename FabProj>
std::enable_if_t<IsBaseFab<FAB>::value>
UnpackRecvBuffers (const ApplyDtosAndProjectionOnReciever<DTOS, FabProj>& packing,
                   FabArray<FAB>& dest, const CommData& recv) {
    // If FAB is not FArrayBox we have not checked for the correct types yet.
    static_assert(IsFabProjection<FabProj, FAB>(), "FabProj needs to be a projection for given FAB type.");
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        unpack_recv_buffer_gpu(dest, packing.dest_component, packing.n_components, recv.data,
                               recv.size, recv.cctc, packing.dtos, packing.proj);
    } else
#endif // AMREX_USE_GPU
    {
        unpack_recv_buffer_cpu(dest, packing.dest_component, packing.n_components, recv.data,
                               recv.size, recv.cctc, packing.dtos, packing.proj);
    }
}
#endif // AMREX_USE_MPI

static_assert(IsDataPacking<ApplyDtosAndProjectionOnReciever<>, FArrayBox>(), // NOLINT(bugprone-throw-keyword-missing)
              "ApplyDtosAndProjectionOnReciever<> is expected to satisfy the DataPacking concept.");

////////////////////////////////////////////////////////////////////////////////////
//                                                             [ParallelCopy_nowait]

static constexpr struct NoLocalCopy {} no_local_copy{};
static constexpr struct DoLocalCopy {} do_local_copy{};

//! Initiate recv and send calls for MPI and immediately return without doing any work.
//!
//! DataPacking is a customization point object to control the behaviour of packing and unpacking
//! send or recv data buffers. It is used to perform interpolation or data transformations on either
//! sender or receiver side.
//!
//! This function performs a data packing on sender side and we expect a call to Parallel_finish
//! that performs data unpacking on the receiver side.
//!
//! \param[out] dest     The Multifab that is going to be filled with received data.
//!
//! \param[in]  src      The Multifab that is used to fill the send buffers.
//!
//! \param[in]  cmd      The communication meta data object holds spatial information about FAB
//!                      boxes that need to be filled and copied from.
//!
//! \param[in]  data_packing A CPO that controls behaviour of preparing buffers and packing the
//!                          source data into the send buffers.
//!
//! \return Returns a CommHandler object that owns context and memory buffers for the whole life
//!         time of the MPI transaction.
template <typename FAB, typename DataPacking,
          typename = std::enable_if_t<IsBaseFab<FAB>::value>,
          typename = std::enable_if_t<IsDataPacking<DataPacking, FAB>::value>>
#ifdef AMREX_USE_MPI
AMREX_NODISCARD CommHandler
ParallelCopy_nowait (NoLocalCopy, FabArray<FAB>& dest, const FabArray<FAB>& src,
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

    if (cmd.m_RcvTags && !(cmd.m_RcvTags->empty())) {
        PrepareRecvBuffers(data_packing, dest, src, handler.recv, *cmd.m_RcvTags);
        PostRecvs(handler.recv, handler.mpi_tag);
    }

    if (cmd.m_SndTags && !(cmd.m_SndTags->empty())) {
        PrepareSendBuffers(data_packing, dest, src, handler.send, *cmd.m_SndTags);
        PackSendBuffers(data_packing, src, handler.send);
        PostSends(handler.send, handler.mpi_tag);
    }
    return handler;
}
#else
CommHandler ParallelCopy_nowait (NoLocalCopy, FabArray<FAB>&, const FabArray<FAB>&,
                                 const FabArrayBase::CommMetaData&, const DataPacking&) {
    return CommHandler{};
}
#endif

//! Initiate recv and send calls for MPI and return after doing local work.
//!
//! DataPacking is a customization point object to control the behaviour of packing and unpacking
//! send or recv data buffers. It is used to perform interpolation or data transformations on either
//! sender or receiver side.
//!
//! This function performs a data packing on sender side and we expect a call to Parallel_finish
//! that performs data unpacking on the receiver side.
//!
//! \param[out] dest     The Multifab that is going to be filled with received data.
//!
//! \param[in]  src      The Multifab that is used to fill the send buffers.
//!
//! \param[in]  cmd      The communication meta data object holds spatial information about FAB
//!                      boxes that need to be filled and copied from.
//!
//! \param[in]  data_packing A CPO that controls behaviour of preparing buffers and packing the
//!                          source data into the send buffers.
//!
//! \return Returns a CommHandler object that owns context and memory buffers for the whole life
//!         time of the MPI transaction.
template <typename FAB, typename DataPacking,
          typename = std::enable_if_t<IsBaseFab<FAB>::value>,
          typename = std::enable_if_t<IsDataPacking<DataPacking, FAB>::value>>
#ifdef AMREX_USE_MPI
AMREX_NODISCARD
#endif
CommHandler
ParallelCopy_nowait (FabArray<FAB>& dest, const FabArray<FAB>& src,
                     const FabArrayBase::CommMetaData& cmd, const DataPacking& data_packing) {
    CommHandler comm = ParallelCopy_nowait(no_local_copy, dest, src, cmd, data_packing);
    // Eagerly do the local work and hope for some overlap with communication
    if (cmd.m_LocTags && !cmd.m_LocTags->empty()) {
        LocalCopy(data_packing, dest, src, *cmd.m_LocTags);
    }
    return comm;
}
////////////////////////////////////////////////////////////////////////////////////
//                                                             [ParallelCopy_finish]

//! Blockingly wait for all communication to be done and fill the local FABs with received data.
//!
//! This function overload performs no local copies, i.e. from this MPI process to itself.
//! It will block the current thread until all MPI recv and send requests are done and
//! calls the DataPacking object to unpack the received buffers.
//!
//! \param[out] dest     The Multifab that is going to be filled with received data.
//!
//! \param[in]  handler  This object holds all data buffers that need to be kept alive as long as
//!                      the data transaction is not done.
//!
//! \param[in]  cmd      The communication meta data object holds spatial information about FAB
//!                      boxes that need to be filled.
//!
//! \param[in]  data_packing A CPO that controls behaviour of unpacking the received buffer to the
//!                          destination FabArray.
//!
//! \returns Nothing.
template <typename FAB, typename DataPacking>
std::enable_if_t<IsBaseFab<FAB>() && IsDataPacking<DataPacking, FAB>()>
#ifdef AMREX_USE_MPI
ParallelCopy_finish (FabArray<FAB>& dest, CommHandler handler,
                     const FabArrayBase::CommMetaData& cmd, const DataPacking& data_packing) {
    // If any FabArray is empty we have nothing to do.
    if (dest.empty()) {
        return;
    }
    // Return if nothing do
    if (ParallelContext::NProcsSub() == 1) {
        return;
    }
    // Unpack receives
    if (cmd.m_RcvTags && !(cmd.m_RcvTags->empty())) {
        ParallelDescriptor::Waitall(handler.recv.request, handler.recv.stats);
#ifdef AMREX_DEBUG
        if (!CheckRcvStats(handler.recv.stats, handler.recv.size, handler.mpi_tag)) {
            amrex::Abort("NonLocalPC::ParallelCopy_finish failed with wrong message size");
        }
#endif
        UnpackRecvBuffers(data_packing, dest, handler.recv);
    }

    // Wait for all sends to be done
    if (cmd.m_SndTags && !(cmd.m_SndTags->empty())) {
        ParallelDescriptor::Waitall(handler.send.request, handler.send.stats);
    }
}
#else
ParallelCopy_finish (FabArray<FAB>&, CommHandler, const FabArrayBase::CommMetaData&, const DataPacking&) {}
#endif

//! Blockingly wait for all communication to be done and fill the local FABs with received data.
//!
//! This function overload performs local copies, i.e. from this MPI process to itself.
//! It will block the current thread until all MPI recv and send requests are done and
//! calls the DataPacking object to unpack the received buffers.
//!
//! \param[out] dest     The Multifab that is going to be filled with received data.
//!
//! \param[in]  src      The Multifab that is used to get the data for the local copies from.
//!
//! \param[in]  handler  This object holds all data buffers that need to be kept alive as long as
//!                      the data transaction is not done.
//!
//! \param[in]  cmd      The communication meta data object holds spatial information about FAB
//!                      boxes that need to be filled.
//!
//! \param[in]  data_packing A CPO that controls behaviour of unpacking the received data.
//!
//! \returns Nothing.
template <typename FAB, typename DataPacking>
std::enable_if_t<IsBaseFab<FAB>() && IsDataPacking<DataPacking, FAB>()>
ParallelCopy_finish (DoLocalCopy, FabArray<FAB>& dest, const FabArray<FAB>& src, CommHandler handler,
                     const FabArrayBase::CommMetaData& cmd, const DataPacking& data_packing) {
    // Eagerly do the local work and hope for some overlap with communication
    if (cmd.m_LocTags && !cmd.m_LocTags->empty()) {
        LocalCopy(data_packing, dest, src, *cmd.m_LocTags);
    }
    ParallelCopy_finish(dest, std::move(handler), cmd, data_packing); // NOLINT
}

//! \brief Call ParallelCopy_nowait followed by ParallelCopy_finish, strong typed version.
//!
//! This function overload uses an already cached CommMetaData. This CommMetaData needs to be
//! compatible with the specified DTOS and projection, otherwise undefined behaviour occurs.
//!
//! \param[out] dest The Multifab that is going to be filled with received data.
//!
//! \param[in]  src  The Multifab that is used to fill the send buffers.
//!
//! \param[in]  cmd  The communication meta data object holds spatial information about FAB
//!                  boxes that need to be filled and copied from.
//!
//! \param[in]  srccomp  The first component in src that will be copied to dest.
//!
//! \param[in]  destcomp The first component in dest that will get written by src.
//!
//! \param[in]  ncomp    The number of successive components that will be copied.
//!
//! \param[in]  dtos An index mapping that maps indices from destination space to source space
//!                  and from source space to destination space.
//!
//! \param[in]  proj A transformation function that might change the data when it is being copied.
//!
//! \return Nothing.
template <typename FAB, typename DTOS = Identity, typename Proj = Identity>
std::enable_if_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
ParallelCopy (FabArray<FAB>& dest, const FabArray<FAB>& src, const FabArrayBase::CommMetaData& cmd,
              SrcComp srccomp, DestComp destcomp, NumComps numcomp, DTOS dtos = DTOS{}, Proj proj = Proj{}) {
    PackComponents components{};
    components.dest_component = destcomp.i;
    components.src_component = srccomp.i;
    components.n_components = numcomp.n;
    ApplyDtosAndProjectionOnReciever<DTOS, Proj> packing{components, std::move(dtos), std::move(proj)};
    CommHandler handler = ParallelCopy_nowait(dest, src, cmd, packing);
    ParallelCopy_finish(dest, std::move(handler), cmd, packing); // NOLINT
}

//! \brief Call ParallelCopy_nowait followed by ParallelCopy_finish.
//!
//! This function overload uses an already cached CommMetaData. This CommMetaData needs to be
//! compatible with the specified DTOS and projection, otherwise undefined behaviour occurs.
//!
//! \param[out] dest The Multifab that is going to be filled with received data.
//!
//! \param[in]  src  The Multifab that is used to fill the send buffers.
//!
//! \param[in]  cmd  The communication meta data object holds spatial information about FAB
//!                  boxes that need to be filled and copied from.
//!
//! \param[in]  srccomp  The first component in src that will be copied to dest.
//!
//! \param[in]  destcomp The first component in dest that will get written by src.
//!
//! \param[in]  ncomp    The number of successive components that will be copied.
//!
//! \param[in]  dtos An index mapping that maps indices from destination space to source space
//!                  and from source space to destination space.
//!
//! \param[in]  proj A transformation function that might change the data when it is being copied.
//!
//! \return Nothing.
template <typename FAB, typename DTOS = Identity, typename Proj = Identity>
std::enable_if_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
ParallelCopy (FabArray<FAB>& dest, const FabArray<FAB>& src, const FabArrayBase::CommMetaData& cmd,
              int srccomp, int destcomp, int numcomp, DTOS dtos = DTOS{}, Proj proj = Proj{}) {
    ParallelCopy(dest, src, cmd, SrcComp(srccomp), DestComp(destcomp), NumComps(numcomp), std::move(dtos), std::move(proj));
}

//! \brief Call ParallelCopy_nowait followed by ParallelCopy_finish, strong typed version.
//!
//! This function constructs a new MultiCommMetaData from the given DTOS, destbox and ngrow.
//!
//! \param[out] dest The Multifab that is going to be filled with received data.
//!
//! \param[in]  destbox The index box in the destination space that will be filled by data from
//!                     src. The source box that describes the dependencies will be computed by
//!                     the specified DTOS.
//!
//! \param[in]  src  The Multifab that is used to fill the send buffers.
//!
//! \param[in]  srccomp  The first component in src that will be copied to dest.
//!
//! \param[in]  destcomp The first component in dest that will get written by src.
//!
//! \param[in]  ncomp    The number of successive components that will be copied.
//!
//! \param[in]  ngrow    The amount of ghost cells that will be taking into consideration.
//!                      Note, even if destbox contains indices outside the domain we need to
//!                      specify an appropriate ngrow that covers the amount of ghost cells that
//!                      we want to copy.
//!
//! \param[in]  dtos An index mapping that maps indices from destination space to source space
//!                  and from source space to destination space.
//!
//! \param[in]  proj A transformation function that might change the data when it is being copied.
//!
//! \return Returns the CommMetaData object that can be cached for future calls to ParallelCopy.
template <typename FAB, typename DTOS = Identity, typename Proj = Identity>
std::enable_if_t<IsBaseFab<FAB>() && IsIndexMapping<DTOS>() && IsFabProjection<Proj, FAB>(),
MultiBlockCommMetaData>
ParallelCopy (FabArray<FAB>& dest, const Box& destbox, const FabArray<FAB>& src, SrcComp srccomp,
              DestComp destcomp, NumComps numcomp, const IntVect& ngrow, DTOS dtos = DTOS{}, Proj proj = Proj{}) {
    MultiBlockCommMetaData cmd(dest, destbox, src, ngrow, dtos);
    ParallelCopy(dest, src, cmd, srccomp, destcomp, numcomp, std::move(dtos), std::move(proj));
    return cmd;
}

//! \brief Call ParallelCopy_nowait followed by ParallelCopy_finish.
//!
//! This function constructs a new MultiCommMetaData from the given DTOS, destbox and ngrow.
//!
//! \param[out] dest The Multifab that is going to be filled with received data.
//!
//! \param[in]  destbox The index box in the destination space that will be filled by data from
//!                     src. The source box that describes the dependencies will be computed by
//!                     the specified DTOS.
//!
//! \param[in]  src  The Multifab that is used to fill the send buffers.
//!
//! \param[in]  srccomp  The first component in src that will be copied to dest.
//!
//! \param[in]  destcomp The first component in dest that will get written by src.
//!
//! \param[in]  ncomp    The number of successive components that will be copied.
//!
//! \param[in]  ngrow    The amount of ghost cells that will be taking into consideration.
//!                      Note, even if destbox contains indices outside the domain we need to
//!                      specify an appropriate ngrow that covers the amount of ghost cells that
//!                      we want to copy.
//!
//! \param[in]  dtos An index mapping that maps indices from destination space to source space
//!                  and from source space to destination space.
//!
//! \param[in]  proj A transformation function that might change the data when it is being copied.
//!
//! \return Returns the CommMetaData object that can be cached for future calls to ParallelCopy.
template <typename FAB, typename DTOS = Identity, typename Proj = Identity>
std::enable_if_t<IsBaseFab<FAB>() && IsIndexMapping<DTOS>() && IsFabProjection<Proj, FAB>(),
MultiBlockCommMetaData>
ParallelCopy (FabArray<FAB>& dest, const Box& destbox, const FabArray<FAB>& src, int srccomp,
              int destcomp, int numcomp, const IntVect& ngrow, DTOS dtos = DTOS{}, Proj proj = Proj{}) {
    return ParallelCopy(dest, destbox, src, SrcComp(srccomp), DestComp(destcomp), NumComps(numcomp), ngrow, std::move(dtos), std::move(proj));
}

// Rotate90 fills the lo-x and lo-y boundary regions by rotating the data
// around (x=0,y=0) by 90 degrees in either direction.  It also fills the
// corner of lo-x and lo-y boundary region by rotating the data by 180
// degrees.

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
Rotate90 (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain);

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
Rotate90 (FabArray<FAB>& mf, Box const& domain);

// Rotate180 fills the lo-x boundary by rotating the data around
// (x=0,y=L_y/2) by 180 degrees.

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
Rotate180 (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain);

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
Rotate180 (FabArray<FAB>& mf, Box const& domain);

// Fill the polar boundaries of the spherical coordinates (theta, phi, r).
// The lo-x boundary is filled with f(-x,y) = f(x,mod(y+pi,2*pi)), and
// the hi-x boundary is filled with f(pi+x,y) = f(pi-x,mod(y+pi,2*pi)).

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
FillPolar (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain);

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
FillPolar (FabArray<FAB>& mf, Box const& domain);

}

#include <AMReX_NonLocalBCImpl.H>

namespace amrex {
    using NonLocalBC::ParallelCopy;
    using NonLocalBC::ParallelCopy_nowait;
    using NonLocalBC::ParallelCopy_finish;
    using NonLocalBC::MultiBlockIndexMapping;
    using NonLocalBC::MultiBlockCommMetaData;
    using NonLocalBC::CommHandler;
}

#endif
