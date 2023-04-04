
#ifndef BL_FABARRAY_H
#define BL_FABARRAY_H
#include <AMReX_Config.H>

#include <AMReX_BLassert.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Box.H>
#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_FabFactory.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Periodicity.H>
#include <AMReX_Print.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_MFIter.H>
#include <AMReX_MakeType.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_LayoutData.H>
#include <AMReX_BaseFabUtility.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_TagParallelFor.H>
#include <AMReX_ParReduce.H>

#include <AMReX_Gpu.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <cstring>
#include <limits>
#include <map>
#include <utility>
#include <vector>
#include <algorithm>
#include <set>
#include <string>

namespace amrex {

template <typename T, typename std::enable_if<!IsBaseFab<T>::value,int>::type = 0>
Long nBytesOwned (T const&) noexcept { return 0; }

template <typename T>
Long nBytesOwned (BaseFab<T> const& fab) noexcept { return fab.nBytesOwned(); }

/*
  A Collection of Fortran Array-like Objects


  The FabArray<FAB> class implements a collection (stored as an array) of
  Fortran array-like objects.  The parameterized type FAB is intended to be
  any class derived from BaseFab<T>.  For example, FAB may be a BaseFab of
  integers, so we could write:

    FabArray<BaseFab<int> > int_fabs;

  Then int_fabs is a FabArray that can hold a collection of BaseFab<int>
  objects.

  FabArray is not just a general container class for Fortran arrays.  It is
  intended to hold "grid" data for use in finite difference calculations in
  which the data is defined on a union of (usually disjoint) rectangular
  regions embedded in a uniform index space.  This region, called the valid
  region, is represented by a BoxArray.  For the purposes of this discussion,
  the Kth Box in the BoxArray represents the interior region of the Kth grid.

  Since the intent is to be used with finite difference calculations a
  FabArray also includes the notion of a boundary region for each grid.  The
  boundary region is specified by the ngrow parameter which tells the FabArray
  to allocate each FAB to be ngrow cells larger in all directions than the
  underlying Box.  The larger region covered by the union of all the FABs is
  called the region of definition.  The underlying notion is that the valid
  region contains the grid interior data and the region of definition includes
  the interior region plus the boundary areas.

  Operations are available to copy data from the valid regions into these
  boundary areas where the two overlap.  The number of components, that is,
  the number of values that can be stored in each cell of a FAB, is either
  given as an argument to the constructor or is inherent in the definition of
  the underlying FAB.  Each FAB in the FabArray will have the same number of
  components.

  In summary, a FabArray is an array of FABs.  The Kth element contains a FAB
  that holds the data for the Kth grid, a Box that defines the valid region
  of the Kth grid.

  A typical use for a FabArray would be to hold the solution vector or
  right-hand-side when solving a linear system of equations on a union of
  rectangular grids.  The copy operations would be used to copy data from the
  valid regions of neighboring grids into the boundary regions after each
  relaxation step of the iterative method.  If a multigrid method is used, a
  FabArray could be used to hold the data at each level in the multigrid
  hierarchy.

  This class is a concrete class not a polymorphic one.

  This class does NOT provide a copy constructor or assignment operator.
*/

//
// alloc: allocate memory or not
//
struct MFInfo {
    bool    alloc = true;
    Arena*  arena = nullptr;
    Vector<std::string> tags;

    MFInfo& SetAlloc (bool a) noexcept { alloc = a; return *this; }

    MFInfo& SetArena (Arena* ar) noexcept { arena = ar; return *this; }

    MFInfo& SetTag () noexcept { return *this; }

    MFInfo& SetTag (const char* t) noexcept {
        tags.emplace_back(t);
        return *this;
    }

    MFInfo& SetTag (const std::string& t) noexcept {
        tags.emplace_back(t);
        return *this;
    }

    template <typename T, typename... Ts>
    MFInfo& SetTag (T&& t, Ts&&... ts) noexcept {
        tags.emplace_back(std::forward<T>(t));
        return SetTag(std::forward<Ts>(ts)...);
    }
};

struct TheFaArenaDeleter {
    using pointer = char*;
    void operator()(pointer p) const noexcept {
        The_FA_Arena()->free(p);
    }
};
using TheFaArenaPointer = std::unique_ptr<char, TheFaArenaDeleter>;

// Data used in non-blocking fill boundary.
template <class FAB>
struct FBData {

    const FabArrayBase::FB*  fb = nullptr;
    int                 scomp;
    int                 ncomp;

    //
    char*               the_recv_data = nullptr;
    char*               the_send_data = nullptr;
    Vector<int>         recv_from;
    Vector<char*>       recv_data;
    Vector<std::size_t> recv_size;
    Vector<MPI_Request> recv_reqs;
    Vector<MPI_Status>  recv_stat;
    //
    Vector<char*>       send_data;
    Vector<MPI_Request> send_reqs;
    int                 tag;

};

// Data used in non-blocking parallel copy.
template <class FAB>
struct PCData {

    const FabArrayBase::CPC*  cpc = nullptr;
    const FabArray<FAB>*      src = nullptr;
    FabArrayBase::CpOp  op;
    int                 tag = -1;
    int                 actual_n_rcvs = -1;
    int                 SC = -1, NC = -1, DC = -1;

    char*               the_recv_data = nullptr;
    char*               the_send_data = nullptr;
    Vector<int>         recv_from;
    Vector<char*>       recv_data;
    Vector<std::size_t> recv_size;
    Vector<MPI_Request> recv_reqs;
    Vector<MPI_Request> send_reqs;

};

template <typename T>
struct MultiArray4
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Array4<T> const& operator[] (int li) const noexcept {
#if AMREX_DEVICE_COMPILE
        return dp[li];
#else
        return hp[li];
#endif
    }

#ifdef AMREX_USE_GPU
    Array4<T> const* AMREX_RESTRICT dp = nullptr;
#endif
    Array4<T> const* AMREX_RESTRICT hp = nullptr;
};

template <class FAB> class FabArray;

template <class DFAB, class SFAB,
          std::enable_if_t<std::conjunction_v<
              IsBaseFab<DFAB>, IsBaseFab<SFAB>,
              std::is_convertible<typename SFAB::value_type,
                                  typename DFAB::value_type>>, int> BAR = 0>
void
Copy (FabArray<DFAB>& dst, FabArray<SFAB> const& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    Copy(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

template <class DFAB, class SFAB,
          std::enable_if_t<std::conjunction_v<
              IsBaseFab<DFAB>, IsBaseFab<SFAB>,
              std::is_convertible<typename SFAB::value_type,
                                  typename DFAB::value_type>>, int> BAR = 0>
void
Copy (FabArray<DFAB>& dst, FabArray<SFAB> const& src, int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
    BL_PROFILE("amrex::Copy()");

    using DT = typename DFAB::value_type;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && dst.isFusingCandidate()) {
        auto const& srcarr = src.const_arrays();
        auto const& dstarr = dst.arrays();
        ParallelFor(dst, nghost, numcomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            dstarr[box_no](i,j,k,dstcomp+n) = DT(srcarr[box_no](i,j,k,srccomp+n));
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            if (bx.ok())
            {
                auto const& srcFab = src.const_array(mfi);
                auto const& dstFab = dst.array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, numcomp, i, j, k, n,
                {
                    dstFab(i,j,k,dstcomp+n) = DT(srcFab(i,j,k,srccomp+n));
                });
            }
        }
    }
}

template <class FAB,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
Add (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    Add(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

template <class FAB,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
Add (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
    BL_PROFILE("amrex::Add()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && dst.isFusingCandidate()) {
        auto const& dstfa = dst.arrays();
        auto const& srcfa = src.const_arrays();
        ParallelFor(dst, nghost, numcomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            dstfa[box_no](i,j,k,n+dstcomp) += srcfa[box_no](i,j,k,n+srccomp);
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            if (bx.ok())
            {
                auto const srcFab = src.array(mfi);
                auto       dstFab = dst.array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, numcomp, i, j, k, n,
                {
                    dstFab(i,j,k,n+dstcomp) += srcFab(i,j,k,n+srccomp);
                });
            }
        }
    }
}

template <class FAB>
class FabArray
    :
    public FabArrayBase
{
public:

    struct FABType {
        typedef FAB value_type;
    };

    /*
    * if FAB is a BaseFab or its child, value_type = FAB::value_type
    * else                              value_type = FAB;
    */
    using value_type = typename std::conditional<IsBaseFab<FAB>::value, FAB, FABType>::type::value_type;

    using fab_type = FAB;

    //
    //! Constructs an empty FabArray<FAB>.
    FabArray () noexcept;

    /**
     * \brief Construct an empty FabArray<FAB> that has a default Arena.  If
     * `define` is called later with a nulltpr as MFInfo's arena, the
     * default Arena `a` will be used.  If the arena in MFInfo is not a
     * nullptr, the MFInfo's arena will be used.
     */
    explicit FabArray (Arena* a) noexcept;

    /**
    * \brief Construct a FabArray<FAB> with a valid region defined by bxs
    * and a region of definition defined by the grow factor ngrow
    * and the number of components nvar.
    */
    FabArray (const BoxArray&            bxs,
              const DistributionMapping& dm,
              int                        nvar,
              int                        ngrow,
#ifdef AMREX_STRICT_MODE
              const MFInfo&              info,
              const FabFactory<FAB>&     factory);
#else
              const MFInfo&              info = MFInfo(),
              const FabFactory<FAB>&     factory = DefaultFabFactory<FAB>());
#endif

    FabArray (const BoxArray&            bxs,
              const DistributionMapping& dm,
              int                        nvar,
              const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
              const MFInfo&              info,
              const FabFactory<FAB>&     factory);
#else
              const MFInfo&              info = MFInfo(),
              const FabFactory<FAB>&     factory = DefaultFabFactory<FAB>());
#endif

    FabArray (const FabArray<FAB>& rhs, MakeType maketype, int scomp, int ncomp);

    //! The destructor -- deletes all FABs in the array.
    ~FabArray ();

    FabArray (FabArray<FAB>&& rhs) noexcept;
    FabArray<FAB>& operator= (FabArray<FAB>&& rhs) noexcept;

    FabArray (const FabArray<FAB>& rhs) = delete;
    FabArray<FAB>& operator= (const FabArray<FAB>& rhs) = delete;

    /**
    * \brief Define this FabArray identically to that performed by
    * the constructor having an analogous function signature.
    * This is only valid if this FabArray was defined using
    * the default constructor.
    */
    void define (const BoxArray& bxs,
                 const DistributionMapping& dm,
                 int                        nvar,
                 int                        ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&              info,
                 const FabFactory<FAB>&     factory);
#else
                 const MFInfo&              info = MFInfo(),
                 const FabFactory<FAB>&     factory = DefaultFabFactory<FAB>());
#endif

    void define (const BoxArray& bxs,
                 const DistributionMapping& dm,
                 int                        nvar,
                 const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&              info,
                 const FabFactory<FAB>&     factory);
#else
                 const MFInfo&              info = MFInfo(),
                 const FabFactory<FAB>&     factory = DefaultFabFactory<FAB>());
#endif

    const FabFactory<FAB>& Factory () const noexcept { return *m_factory; }

    // Provides access to the Arena this FabArray was build with.
    Arena* arena () const noexcept { return m_dallocator.arena(); }

    const Vector<std::string>& tags () const noexcept { return m_tags; }

    bool hasEBFabFactory () const noexcept {
#ifdef AMREX_USE_EB
        const auto *const f = dynamic_cast<EBFArrayBoxFactory const*>(m_factory.get());
        return (f != nullptr);
#else
        return false;
#endif
    }

    bool isAllRegular () const noexcept {
#ifdef AMREX_USE_EB
        const auto *const f = dynamic_cast<EBFArrayBoxFactory const*>(m_factory.get());
        if (f) {
            return f->isAllRegular();
        } else {
            return true;
        }
#else
        return true;
#endif
    }

    /**
    * \brief Return true if the FabArray is well-defined.  That is,
    * the FabArray has a BoxArray and DistributionMapping, the
    * FABs are allocated for each Box in the BoxArray and the
    * sizes of the FABs and the number of components are consistent
    * with the definition of the FabArray.
    */
    bool ok () const;

    /** Has define() been called on this rank?
     *
     * \return true if `define` has been called on this `FabArray`.  Note that all constructors except `FabArray ()`
     * and `FabArray(Arena*a)` call `define`, even if the `MFInfo` argument has `alloc=false`.  One could
     * also use `FabArrayBase::empty()` to find whether `define` is called or not, although they are not exactly
     * the same.
     */
    bool isDefined () const;

    //! Return a constant reference to the FAB associated with mfi.
    const FAB& operator[] (const MFIter& mfi) const noexcept { return *(this->fabPtr(mfi)); }

    //! Return a constant reference to the FAB associated with mfi.
    const FAB& get (const MFIter& mfi) const noexcept { return *(this->fabPtr(mfi)); }

    //! Returns a reference to the FAB associated mfi.
    FAB& operator[] (const MFIter& mfi) noexcept { return *(this->fabPtr(mfi)); }

    //! Returns a reference to the FAB associated mfi.
    FAB& get (const MFIter& mfi) noexcept { return *(this->fabPtr(mfi)); }

    //! Return a constant reference to the FAB associated with the Kth element.
    const FAB& operator[] (int K) const noexcept { return *(this->fabPtr(K)); }

    //! Return a constant reference to the FAB associated with the Kth element.
    const FAB& get (int K) const noexcept { return *(this->fabPtr(K)); }

    //! Return a reference to the FAB associated with the Kth element.
    FAB& operator[] (int K) noexcept { return *(this->fabPtr(K)); }

    //! Return a reference to the FAB associated with the Kth element.
    FAB& get (int K) noexcept { return *(this->fabPtr(K)); }

    //! Return a reference to the FAB associated with local index L
    FAB& atLocalIdx (int L) noexcept { return *m_fabs_v[L]; }
    const FAB& atLocalIdx (int L) const noexcept { return *m_fabs_v[L]; }

    //! Return pointer to FAB
    FAB      * fabPtr (const MFIter& mfi) noexcept;
    FAB const* fabPtr (const MFIter& mfi) const noexcept;
    FAB      * fabPtr (int K) noexcept;  // Here K is global index
    FAB const* fabPtr (int K) const noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void prefetchToHost (const MFIter& mfi) const noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void prefetchToDevice (const MFIter& mfi) const noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> array (const MFIter& mfi) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type> array (const MFIter& mfi) noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> array (int K) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type> array (int K) noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> const_array (const MFIter& mfi) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> const_array (int K) const noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> array (const MFIter& mfi, int start_comp) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type> array (const MFIter& mfi, int start_comp) noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> array (int K, int start_comp) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type> array (int K, int start_comp) noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> const_array (const MFIter& mfi, int start_comp) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> const_array (int K, int start_comp) const noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    MultiArray4<typename FabArray<FAB>::value_type> arrays () noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    MultiArray4<typename FabArray<FAB>::value_type const> arrays () const noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    MultiArray4<typename FabArray<FAB>::value_type const> const_arrays () const noexcept;

    //! Explicitly set the Kth FAB in the FabArray to point to elem.
    void setFab (int boxno, std::unique_ptr<FAB> elem);

    //! Explicitly set the Kth FAB in the FabArray to point to elem.
    template <class F=FAB, std::enable_if_t<std::is_move_constructible<F>::value,int> = 0>
    void setFab (int boxno, FAB&& elem);

    //! Explicitly set the FAB associated with mfi in the FabArray to point to elem.
    void setFab (const MFIter&mfi, std::unique_ptr<FAB> elem);

    //! Explicitly set the FAB associated with mfi in the FabArray to point to elem.
    template <class F=FAB, std::enable_if_t<std::is_move_constructible<F>::value,int> = 0>
    void setFab (const MFIter&mfi, FAB&& elem);

    //! Release ownership of the FAB. This function is not thread safe.
    AMREX_NODISCARD
    FAB* release (int K);

    //! Release ownership of the FAB. This function is not thread safe.
    AMREX_NODISCARD
    FAB* release (const MFIter& mfi);

    //! Releases FAB memory in the FabArray.
    void clear ();

    /**
     * \brief Perform local copy of FabArray data.
     *
     * The two FabArrays must have the same BoxArray and
     * DistributionMapping, although they could have different data types.
     * For example, this could be used to copy from FabArray<BaseFab<float>>
     * to FabArray<BaseFab<double>>.
     *
     * \param src    source FabArray
     * \param scomp  starting component of source
     * \param dcomp  starting component of this FabArray
     * \param ncomp  number of components
     * \param nghost number of ghost cells
     */
    template <typename SFAB, typename DFAB = FAB,
              std::enable_if_t<std::conjunction_v<
                  IsBaseFab<DFAB>, IsBaseFab<SFAB>,
                  std::is_convertible<typename SFAB::value_type,
                                      typename DFAB::value_type>>, int> = 0>
    void LocalCopy (FabArray<SFAB> const& src, int scomp, int dcomp, int ncomp,
                    IntVect const& nghost);

    /**
     * \brief Perform local addition of FabArray data.
     *
     * The two FabArrays must have the same BoxArray and
     * DistributionMapping.
     *
     * \param src    source FabArray
     * \param scomp  starting component of source
     * \param dcomp  starting component of this FabArray
     * \param ncomp  number of components
     * \param nghost number of ghost cells
     */
    template <class F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    void LocalAdd (FabArray<FAB> const& src, int scomp, int dcomp, int ncomp,
                   IntVect const& nghost);

    //! Set all components in the entire region of each FAB to val.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val);

    //! Set all components in the entire region of each FAB to val.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    FabArray<FAB>& operator= (value_type val);

    /**
    * \brief Set the value of num_comp components in the valid region of
    * each FAB in the FabArray, starting at component comp to val.
    * Also set the value of nghost boundary cells.
    */
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val,
                 int        comp,
                 int        ncomp,
                 int        nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val,
                 int        comp,
                 int        ncomp,
                 const IntVect& nghost);

    /**
    * \brief Set the value of num_comp components in the valid region of
    * each FAB in the FabArray, starting at component comp, as well
    * as nghost boundary cells, to val, provided they also intersect
    * with the Box region.
    */
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val,
                 const Box& region,
                 int        comp,
                 int        ncomp,
                 int        nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val,
                 const Box& region,
                 int        comp,
                 int        ncomp,
                 const IntVect& nghost);
    /**
    * \brief Set all components in the valid region of each FAB in the
    * FabArray to val, including nghost boundary cells.
    */
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val, int nghost);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val, const IntVect& nghost);

    /**
    * \brief Set all components in the valid region of each FAB in the
    * FabArray to val, including nghost boundary cells, that also
    * intersect the Box region.
    */
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val, const Box& region, int nghost);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val, const Box& region, const IntVect& nghost);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void abs (int comp, int ncomp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void abs (int comp, int ncomp, const IntVect& nghost);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void plus (value_type val, int comp, int num_comp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void plus (value_type val, const Box& region, int comp, int num_comp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void mult (value_type val, int comp, int num_comp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void mult (value_type val, const Box& region, int comp, int num_comp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void invert (value_type numerator, int comp, int num_comp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void invert (value_type numerator, const Box& region, int comp, int num_comp, int nghost = 0);

    //! Set all values in the boundary region to val.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setBndry (value_type val);

    //! Set ncomp values in the boundary region, starting at start_comp to val.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setBndry (value_type val, int strt_comp, int ncomp);

   //! Set all values outside the Geometry domain to val.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setDomainBndry (value_type val, const Geometry& geom);

    //! Set ncomp values outside the Geometry domain to val, starting at start_comp.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setDomainBndry (value_type val, int strt_comp, int ncomp, const Geometry& geom);

    /**
    * \brief Returns the sum of component "comp"
    *
    * \param comp   component
    * \param nghost number of ghost cells
    * \param local  If true, MPI communciation is skipped.
    */
    template <typename F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    typename F::value_type
    sum (int comp, IntVect const& nghost, bool local = false) const;

    /**
    * \brief This function copies data from fa to this FabArray.  Each FAB
    * in fa is intersected with all FABs in this FabArray and a copy
    * is performed on the region of intersection.  The intersection
    * is restricted to the valid regions.
    */
    void ParallelAdd (const FabArray<FAB>& fa,
                      const Periodicity&   period = Periodicity::NonPeriodic())
       { ParallelCopy(fa,period,FabArray::ADD); }
    void ParallelCopy (const FabArray<FAB>& fa,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       CpOp                 op = FabArrayBase::COPY)
       { ParallelCopy(fa,0,0,nComp(),0,0,period,op); }

    [[deprecated("Use FabArray::ParallelCopy() instead.")]]
    void copy (const FabArray<FAB>& fa,
               const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY)
        { ParallelCopy(fa,period,op); }

    void ParallelAdd_nowait (const FabArray<FAB>& fa,
                             const Periodicity& period = Periodicity::NonPeriodic())
        { ParallelCopy_nowait(fa,period,FabArray::ADD); }
    void ParallelCopy_nowait (const FabArray<FAB>& fa,
                              const Periodicity&   period = Periodicity::NonPeriodic(),
                              CpOp                 op = FabArrayBase::COPY)
       { ParallelCopy_nowait(fa,0,0,nComp(),0,0,period,op); }

    /**
    * \brief This function copies data from src to this FabArray.  Each FAB
    * in src is intersected with all FABs in this FabArray and a copy
    * is performed on the region of intersection.  The intersection
    * is restricted to the num_comp components starting at src_comp
    * in the FabArray src, with the destination components in this
    * FabArray starting at dest_comp.
    */
    void ParallelAdd (const FabArray<FAB>& src,
                      int                  src_comp,
                      int                  dest_comp,
                      int                  num_comp,
                      const Periodicity&   period = Periodicity::NonPeriodic())
       { ParallelCopy(src,src_comp,dest_comp,num_comp, period, FabArrayBase::ADD); }
    void ParallelCopy (const FabArray<FAB>& src,
                       int                  src_comp,
                       int                  dest_comp,
                       int                  num_comp,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       CpOp                 op = FabArrayBase::COPY)
       { ParallelCopy(src,src_comp,dest_comp,num_comp,0,0,period,op); }

    [[deprecated("Use FabArray::ParallelCopy() instead.")]]
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
               const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY)
        { ParallelCopy(src,src_comp,dest_comp,num_comp, period, op); }

    void ParallelAdd_nowait (const FabArray<FAB>& src,
                             int                  src_comp,
                             int                  dest_comp,
                             int                  num_comp,
                             const Periodicity&   period = Periodicity::NonPeriodic())
       { ParallelCopy_nowait(src,src_comp,dest_comp,num_comp, period, FabArrayBase::ADD); }
    void ParallelCopy_nowait (const FabArray<FAB>& src,
                              int                  src_comp,
                              int                  dest_comp,
                              int                  num_comp,
                              const Periodicity&   period = Periodicity::NonPeriodic(),
                              CpOp                 op = FabArrayBase::COPY)
       { ParallelCopy_nowait(src,src_comp,dest_comp,num_comp,0,0,period,op); }

    //! Similar to the above function, except that source and destination are grown by src_nghost and dst_nghost, respectively
    void ParallelAdd (const FabArray<FAB>& src,
                      int                  src_comp,
                      int                  dest_comp,
                      int                  num_comp,
                      int                  src_nghost,
                      int                  dst_nghost,
                      const Periodicity&   period = Periodicity::NonPeriodic())
       { ParallelCopy(src,src_comp,dest_comp,num_comp,IntVect(src_nghost),IntVect(dst_nghost),period,
                      FabArrayBase::ADD); }
    void ParallelAdd (const FabArray<FAB>& src,
                      int                  src_comp,
                      int                  dest_comp,
                      int                  num_comp,
                      const IntVect&       src_nghost,
                      const IntVect&       dst_nghost,
                      const Periodicity&   period = Periodicity::NonPeriodic())
       { ParallelCopy(src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost,period,FabArrayBase::ADD); }
    void ParallelCopy (const FabArray<FAB>& src,
                       int                  src_comp,
                       int                  dest_comp,
                       int                  num_comp,
                       int                  src_nghost,
                       int                  dst_nghost,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       CpOp                 op = FabArrayBase::COPY)
       { ParallelCopy(src,src_comp,dest_comp,num_comp,IntVect(src_nghost),IntVect(dst_nghost),period,op); }
    void ParallelCopy (const FabArray<FAB>& src,
                       int                  scomp,
                       int                  dcomp,
                       int                  ncomp,
                       const IntVect&       snghost,
                       const IntVect&       dnghost,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       CpOp                 op = FabArrayBase::COPY,
                       const FabArrayBase::CPC* a_cpc = nullptr);

    void ParallelAdd_nowait (const FabArray<FAB>& src,
                             int                  src_comp,
                             int                  dest_comp,
                             int                  num_comp,
                             int                  src_nghost,
                             int                  dst_nghost,
                             const Periodicity&   period = Periodicity::NonPeriodic())
       { ParallelCopy_nowait(src,src_comp,dest_comp,num_comp,IntVect(src_nghost),
                               IntVect(dst_nghost),period,FabArrayBase::ADD); }

    void ParallelAdd_nowait (const FabArray<FAB>& src,
                             int                  src_comp,
                             int                  dest_comp,
                             int                  num_comp,
                             const IntVect&       src_nghost,
                             const IntVect&       dst_nghost,
                             const Periodicity&   period = Periodicity::NonPeriodic())
       { ParallelCopy_nowait(src,src_comp,dest_comp,num_comp,src_nghost,
                             dst_nghost,period,FabArrayBase::ADD); }

    void ParallelCopy_nowait (const FabArray<FAB>& src,
                              int                  src_comp,
                              int                  dest_comp,
                              int                  num_comp,
                              int                  src_nghost,
                              int                  dst_nghost,
                              const Periodicity&   period = Periodicity::NonPeriodic(),
                              CpOp                 op = FabArrayBase::COPY)
       { ParallelCopy_nowait(src,src_comp,dest_comp,num_comp,IntVect(src_nghost),
                             IntVect(dst_nghost),period,op); }

    void ParallelCopy_nowait (const FabArray<FAB>& src,
                              int                  scomp,
                              int                  dcomp,
                              int                  ncomp,
                              const IntVect&       snghost,
                              const IntVect&       dnghost,
                              const Periodicity&   period = Periodicity::NonPeriodic(),
                              CpOp                 op = FabArrayBase::COPY,
                              const FabArrayBase::CPC* a_cpc = nullptr,
                              bool                 to_ghost_cells_only = false);

    void ParallelCopy_finish ();

    void ParallelCopyToGhost (const FabArray<FAB>& src,
                              int                  scomp,
                              int                  dcomp,
                              int                  ncomp,
                              const IntVect&       snghost,
                              const IntVect&       dnghost,
                              const Periodicity&   period = Periodicity::NonPeriodic());

    void ParallelCopyToGhost_nowait (const FabArray<FAB>& src,
                                     int                  scomp,
                                     int                  dcomp,
                                     int                  ncomp,
                                     const IntVect&       snghost,
                                     const IntVect&       dnghost,
                                     const Periodicity&   period = Periodicity::NonPeriodic());

    void ParallelCopyToGhost_finish();

    [[deprecated("Use FabArray::ParallelCopy() instead.")]]
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
               int                  src_nghost,
               int                  dst_nghost,
               const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY)
       { ParallelCopy(src,src_comp,dest_comp,num_comp,IntVect(src_nghost),IntVect(dst_nghost),period,op); }

    [[deprecated("Use FabArray::ParallelCopy() instead.")]]
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
               const IntVect&       src_nghost,
               const IntVect&       dst_nghost,
               const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY)
        { ParallelCopy(src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost,period,op); }

    //! Copy from src to this.  this and src have the same BoxArray, but different DistributionMapping
    void Redistribute (const FabArray<FAB>& src,
                       int                  scomp,
                       int                  dcomp,
                       int                  ncomp,
                       const IntVect&       nghost);

    /**
    * \brief Copy the values contained in the intersection of the
    * valid + nghost region of this FabArray with the FAB dest into dest.
    * Note that FAB dest is assumed to be identical on each process.
    */
    void copyTo (FAB& dest, int  nghost = 0) const;

    /**
    * \brief Copy the values contained in the intersection of the
    * num_comp component valid + nghost region of this FabArray, starting at
    * component src_comp, with the FAB dest into dest, starting at
    * component dest_comp in dest.
    * Note that FAB dest is assumed to be identical on each process.
    */
    void copyTo (FAB& dest, int  scomp, int  dcomp, int  ncomp, int  nghost = 0) const;

    //! Shift the boxarray by vector v
    void shift (const IntVect& v);

    bool defined (int K) const noexcept;
    bool defined (const MFIter& mfi) const noexcept;

    /**
    * \brief Copy on intersection within a FabArray.  Data is copied from
    * valid regions to intersecting regions of definition.  The
    * purpose is to fill in the boundary regions of each FAB in
    * the FabArray.  If cross=true, corner cells are not filled.
    * If the length of periodic is provided, periodic boundaries are
    * also filled.  Note that FabArray itself does not contains
    * any periodicity information.
    * FillBoundary expects that its cell-centered version of its BoxArray
    * is non-overlapping.
    */
    template <typename BUF=value_type>
    void FillBoundary (bool cross = false);

    template <typename BUF=value_type>
    void FillBoundary (const Periodicity& period, bool cross = false);

    template <typename BUF=value_type>
    void FillBoundary (const IntVect& nghost, const Periodicity& period, bool cross = false);

    //! Same as FillBoundary(), but only copies ncomp components starting at scomp.
    template <typename BUF=value_type>
    void FillBoundary (int scomp, int ncomp, bool cross = false);

    template <typename BUF=value_type>
    void FillBoundary (int scomp, int ncomp, const Periodicity& period, bool cross = false);

    template <typename BUF=value_type>
    void FillBoundary (int scomp, int ncomp, const IntVect& nghost, const Periodicity& period, bool cross = false);

    template <typename BUF=value_type>
    void FillBoundary_nowait (bool cross = false);

    template <typename BUF=value_type>
    void FillBoundary_nowait (const Periodicity& period, bool cross = false);

    template <typename BUF=value_type>
    void FillBoundary_nowait (const IntVect& nghost, const Periodicity& period, bool cross = false);

    template <typename BUF=value_type>
    void FillBoundary_nowait (int scomp, int ncomp, bool cross = false);

    template <typename BUF=value_type>
    void FillBoundary_nowait (int scomp, int ncomp, const Periodicity& period, bool cross = false);

    template <typename BUF=value_type>
    void FillBoundary_nowait (int scomp, int ncomp, const IntVect& nghost, const Periodicity& period, bool cross = false);

    template <typename BUF=value_type,
              class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void FillBoundary_finish ();

    void FillBoundary_test ();


    /**
     * \brief Fill ghost cells and synchronize nodal data. Ghost regions are
     * filled with data from the intersecting valid regions. The
     * synchronization will override valid regions by the intersecting valid
     * regions with a higher precedence.  The smaller the global box index
     * is, the higher precedence the box has.  With periodic boundaries, for
     * cells in the same box, those near the lower corner have higher
     * precedence than those near the upper corner.
     *
     * \param period periodic length if it's non-zero
     */
    void FillBoundaryAndSync (const Periodicity& period = Periodicity::NonPeriodic());
    /**
     * \brief Fill ghost cells and synchronize nodal data. Ghost regions are
     * filled with data from the intersecting valid regions. The
     * synchronization will override valid regions by the intersecting valid
     * regions with a higher precedence.  The smaller the global box index
     * is, the higher precedence the box has.  With periodic boundaries, for
     * cells in the same box, those near the lower corner have higher
     * precedence than those near the upper corner.
     *
     * \param scomp starting component
     * \param ncomp number of components
     * \param nghost number of ghost cells to fill
     * \param period periodic length if it's non-zero
     */
    void FillBoundaryAndSync (int scomp, int ncomp, const IntVect& nghost,
                              const Periodicity& period);
    void FillBoundaryAndSync_nowait (const Periodicity& period = Periodicity::NonPeriodic());
    void FillBoundaryAndSync_nowait (int scomp, int ncomp, const IntVect& nghost,
                                     const Periodicity& period);
    void FillBoundaryAndSync_finish ();

    /**
     * \brief Synchronize nodal data.  The synchronization will override
     * valid regions by the intersecting valid regions with a higher
     * precedence.  The smaller the global box index is, the higher
     * precedence the box has.  With periodic boundaries, for cells in the
     * same box, those near the lower corner have higher precedence than
     * those near the upper corner.
     *
     * \param period periodic length if it's non-zero
     */
    void OverrideSync (const Periodicity& period = Periodicity::NonPeriodic());
    /**
     * \brief Synchronize nodal data.  The synchronization will override
     * valid regions by the intersecting valid regions with a higher
     * precedence.  The smaller the global box index is, the higher
     * precedence the box has.  With periodic boundaries, for cells in the
     * same box, those near the lower corner have higher precedence than
     * those near the upper corner.
     *
     * \param scomp starting component
     * \param ncomp number of components
     * \param nghost number of ghost cells to fill
     * \param period periodic length if it's non-zero
     */
    void OverrideSync (int scomp, int ncomp, const Periodicity& period);
    void OverrideSync_nowait (const Periodicity& period = Periodicity::NonPeriodic());
    void OverrideSync_nowait (int scomp, int ncomp, const Periodicity& period);
    void OverrideSync_finish ();

    /**
    * \brief Sum values in overlapped cells.  The destination is limited to valid cells.
    */
    void SumBoundary (const Periodicity& period = Periodicity::NonPeriodic());
    void SumBoundary (int scomp, int ncomp, const Periodicity& period = Periodicity::NonPeriodic());
    void SumBoundary_nowait (const Periodicity& period = Periodicity::NonPeriodic());
    void SumBoundary_nowait (int scomp, int ncomp, const Periodicity& period = Periodicity::NonPeriodic());

    /**
    * \brief Sum values in overlapped cells.  The destination is limited to valid + ngrow cells.
    */
    void SumBoundary (int scomp, int ncomp, IntVect const& nghost,
                      const Periodicity& period = Periodicity::NonPeriodic());
    void SumBoundary_nowait (int scomp, int ncomp, IntVect const& nghost,
                             const Periodicity& period = Periodicity::NonPeriodic());

    /**
    * \brief Sum values in overlapped cells.
    *        For computing the overlap, the dst is grown by dst_ngrow, while the src uses src_ngrow.
    */
    void SumBoundary (int scomp, int ncomp, IntVect const& src_nghost, IntVect const& dst_nghost,
                      const Periodicity& period = Periodicity::NonPeriodic());
    void SumBoundary_nowait (int scomp, int ncomp, IntVect const& src_nghost, IntVect const& dst_nghost,
                             const Periodicity& period = Periodicity::NonPeriodic());
    void SumBoundary_finish ();

    /** \brief Fill cells outside periodic domains with their corresponding cells inside
    * the domain.  Ghost cells are treated the same as valid cells.  The BoxArray
    * is allowed to be overlapping.
    */
    void EnforcePeriodicity (const Periodicity& period);
    void EnforcePeriodicity (int scomp, int ncomp, const Periodicity& period);
    void EnforcePeriodicity (int scomp, int ncomp, const IntVect& nghost,
                             const Periodicity& period);

    // covered   : ghost cells covered by valid cells of this FabArray
    //             (including periodically shifted valid cells)
    // notcovered: ghost cells not covered by valid cells
    //             (including ghost cells outside periodic boundaries)
    // physbnd   : boundary cells outside the domain (excluding periodic boundaries)
    // interior  : interior cells (i.e., valid cells)
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void BuildMask (const Box& phys_domain, const Periodicity& period,
                    value_type covered, value_type notcovered,
                    value_type physbnd, value_type interior);

    // The following are private functions.  But we have to make them public for cuda.

    template <typename BUF=value_type,
              class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void FBEP_nowait (int scomp, int ncomp, const IntVect& nghost,
                      const Periodicity& period, bool cross,
                      bool enforce_periodicity_only = false,
                      bool override_sync = false);

    void FB_local_copy_cpu (const FB& TheFB, int scomp, int ncomp);
    void PC_local_cpu (const CPC& thecpc, FabArray<FAB> const& src,
                       int scomp, int dcomp, int ncomp, CpOp op);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val, const CommMetaData& thecmd, int scomp, int ncomp);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    LayoutData<int> RecvLayoutMask (const CommMetaData& thecmd);

#ifdef AMREX_USE_GPU

    void FB_local_copy_gpu (const FB& TheFB, int scomp, int ncomp);
    void PC_local_gpu (const CPC& thecpc, FabArray<FAB> const& src,
                       int scomp, int dcomp, int ncomp, CpOp op);

    void CMD_local_setVal_gpu (value_type x, const CommMetaData& thecmd, int scomp, int ncomp);
    void CMD_remote_setVal_gpu (value_type x, const CommMetaData& thecmd, int scomp, int ncomp);

#if defined(__CUDACC__)

    void FB_local_copy_cuda_graph_1 (const FB& TheFB, int scomp, int ncomp);
    void FB_local_copy_cuda_graph_n (const FB& TheFB, int scomp, int ncomp);

#endif
#endif

#ifdef AMREX_USE_MPI

#ifdef AMREX_USE_GPU
#if defined(__CUDACC__)

    void FB_pack_send_buffer_cuda_graph (const FB& TheFB, int scomp, int ncomp,
                                         Vector<char*>& send_data,
                                         Vector<std::size_t> const& send_size,
                                         Vector<const CopyComTagsContainer*> const& send_cctc);

    void FB_unpack_recv_buffer_cuda_graph (const FB& TheFB, int dcomp, int ncomp,
                                           Vector<char*> const& recv_data,
                                           Vector<std::size_t> const& recv_size,
                                           Vector<const CopyComTagsContainer*> const& recv_cctc,
                                           bool is_thread_safe);

#endif

    template <typename BUF = value_type>
    static void pack_send_buffer_gpu (FabArray<FAB> const& src, int scomp, int ncomp,
                                      Vector<char*> const& send_data,
                                      Vector<std::size_t> const& send_size,
                                      Vector<const CopyComTagsContainer*> const& send_cctc);

    template <typename BUF = value_type>
    static void unpack_recv_buffer_gpu (FabArray<FAB>& dst, int dcomp, int ncomp,
                                        Vector<char*> const& recv_data,
                                        Vector<std::size_t> const& recv_size,
                                        Vector<const CopyComTagsContainer*> const& recv_cctc,
                                        CpOp op, bool is_thread_safe);

#endif

    template <typename BUF = value_type>
    static void pack_send_buffer_cpu (FabArray<FAB> const& src, int scomp, int ncomp,
                                      Vector<char*> const& send_data,
                                      Vector<std::size_t> const& send_size,
                                      Vector<const CopyComTagsContainer*> const& send_cctc);

    template <typename BUF = value_type>
    static void unpack_recv_buffer_cpu (FabArray<FAB>& dst, int dcomp, int ncomp,
                                        Vector<char*> const& recv_data,
                                        Vector<std::size_t> const& recv_size,
                                        Vector<const CopyComTagsContainer*> const& recv_cctc,
                                        CpOp op, bool is_thread_safe);

#endif

    /**
     * \brief Return infinity norm
     *
     * \param comp           starting component
     * \param ncomp          number of components
     * \param nghost         number of ghost cells
     * \param local          If true, MPI communciation is skipped.
     * \param ignore_covered ignore covered cells. Only relevant for cell-centered EB data.
     */
    template <typename F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    typename F::value_type
    norminf (int comp, int ncomp, IntVect const& nghost, bool local = false,
             [[maybe_unused]] bool ignore_covered = false) const;

    /**
     * \brief Return infinity norm in masked region
     *
     * \param mask   only mask=true region is included
     * \param comp   starting component
     * \param ncomp  number of components
     * \param nghost number of ghost cells
     * \param local  If true, MPI communciation is skipped.
     */
    template <typename IFAB, typename F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    typename F::value_type
    norminf (FabArray<IFAB> const& mask, int comp, int ncomp, IntVect const& nghost,
             bool local = false) const;

protected:

    std::unique_ptr<FabFactory<FAB> > m_factory;
    DataAllocator m_dallocator;

    //! has define() been called?
    bool define_function_called = false;

    //
    //! The data.
    std::vector<FAB*> m_fabs_v;

#ifdef AMREX_USE_GPU
    mutable void* m_dp_arrays = nullptr;
#endif
    mutable void* m_hp_arrays = nullptr;
    mutable MultiArray4<value_type> m_arrays;
    mutable MultiArray4<value_type const> m_const_arrays;

    Vector<std::string> m_tags;

    //! for shared memory
    struct ShMem {

        ShMem () noexcept = default;

        ~ShMem () { // NOLINT
#if defined(BL_USE_MPI3)
            if (win != MPI_WIN_NULL) MPI_Win_free(&win);
#endif
#ifdef BL_USE_TEAM
            if (alloc) {
                amrex::update_fab_stats(-n_points, -n_values, sizeof(value_type));
            }
#endif
        }
        ShMem (ShMem&& rhs) noexcept
                 : alloc(rhs.alloc), n_values(rhs.n_values), n_points(rhs.n_points)
#if defined(BL_USE_MPI3)
                 , win(rhs.win)
#endif
        {
            rhs.alloc = false;
#if defined(BL_USE_MPI3)
            rhs.win = MPI_WIN_NULL;
#endif
        }
        ShMem& operator= (ShMem&& rhs) noexcept {
            if (&rhs != this) {
                alloc = rhs.alloc;
                n_values = rhs.n_values;
                n_points = rhs.n_points;
                rhs.alloc = false;
#if defined(BL_USE_MPI3)
                win = rhs.win;
                rhs.win = MPI_WIN_NULL;
#endif
            }
            return *this;
        }
        ShMem (const ShMem&) = delete;
        ShMem& operator= (const ShMem&) = delete;
        bool  alloc{false};
        Long  n_values{0};
        Long  n_points{0};
#if defined(BL_USE_MPI3)
        MPI_Win win = MPI_WIN_NULL;
#endif
    };
    ShMem shmem;

    bool SharedMemory () const noexcept { return shmem.alloc; }

private:
    typedef typename std::vector<FAB*>::iterator    Iterator;

    void AllocFabs (const FabFactory<FAB>& factory, Arena* ar,
                    const Vector<std::string>& tags);

    void setFab_assert (int K, FAB const& fab) const;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void build_arrays () const;

    void clear_arrays ();

public:

#ifdef BL_USE_MPI

    //! Prepost nonblocking receives
    template <typename BUF=value_type>
    void PostRcvs (const MapOfCopyComTagContainers&       RcvTags,
                   char*&                                 the_recv_data,
                   Vector<char*>&                         recv_data,
                   Vector<std::size_t>&                   recv_size,
                   Vector<int>&                           recv_from,
                   Vector<MPI_Request>&                   recv_reqs,
                   int                                    ncomp,
                   int                                    SeqNum) const;

    template <typename BUF=value_type>
    AMREX_NODISCARD TheFaArenaPointer PostRcvs (const MapOfCopyComTagContainers&       RcvTags,
                   Vector<char*>&                         recv_data,
                   Vector<std::size_t>&                   recv_size,
                   Vector<int>&                           recv_from,
                   Vector<MPI_Request>&                   recv_reqs,
                   int                                    ncomp,
                   int                                    SeqNum) const;

    template <typename BUF=value_type>
    void PrepareSendBuffers (const MapOfCopyComTagContainers&     SndTags,
                             char*&                               the_send_data,
                             Vector<char*>&                       send_data,
                             Vector<std::size_t>&                 send_size,
                             Vector<int>&                         send_rank,
                             Vector<MPI_Request>&                 send_reqs,
                             Vector<const CopyComTagsContainer*>& send_cctc,
                             int                                  ncomp) const;

    template <typename BUF=value_type>
    AMREX_NODISCARD TheFaArenaPointer PrepareSendBuffers (const MapOfCopyComTagContainers&     SndTags,
                             Vector<char*>&                       send_data,
                             Vector<std::size_t>&                 send_size,
                             Vector<int>&                         send_rank,
                             Vector<MPI_Request>&                 send_reqs,
                             Vector<const CopyComTagsContainer*>& send_cctc,
                             int                                  ncomp) const;

    static void PostSnds (Vector<char*> const&       send_data,
                          Vector<std::size_t> const& send_size,
                          Vector<int> const&         send_rank,
                          Vector<MPI_Request>&       send_reqs,
                          int                        SeqNum);
#endif

    std::unique_ptr<FBData<FAB>> fbd;
    std::unique_ptr<PCData<FAB>> pcd;

    // Pointer to temporary fab used in non-blocking amrex::OverrideSync
    std::unique_ptr< FabArray<FAB> > os_temp;



    /**
     * \brief y += a*x
     *
     * \param y      FabArray y
     * \param a      scalar a
     * \param x      FabArray x
     * \param xcomp  starting component of x
     * \param ycomp  starting component of y
     * \param ncomp  number of components
     * \param nghost number of ghost cells
     */
    template <class F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    static void Saxpy (FabArray<FAB>& y, value_type a, FabArray<FAB> const& x,
                       int xcomp, int ycomp, int ncomp, IntVect const& nghost);

    /**
     * \brief y = x + a*y
     *
     * \param y      FabArray y
     * \param a      scalar a
     * \param x      FabArray x
     * \param xcomp  starting component of x
     * \param ycomp  starting component of y
     * \param ncomp  number of components
     * \param nghost number of ghost cells
     */
    template <class F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    static void Xpay (FabArray<FAB>& y, value_type a, FabArray<FAB> const& x,
                      int xcomp, int ycomp, int ncomp, IntVect const& nghost);

    /**
     * \brief dst = a*x + b*y
     *
     * \param dst     destination FabArray
     * \param a       scalar a
     * \param x       FabArray x
     * \param xcomp   starting component of x
     * \param b       scalar b
     * \param y       FabArray y
     * \param ycomp   starting component of y
     * \param dstcomp starting component of destination
     * \param numcomp number of components
     * \param nghost  number of ghost cells
     */
    template <class F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    static void LinComb (FabArray<FAB>& dst,
                         value_type a, const FabArray<FAB>& x, int xcomp,
                         value_type b, const FabArray<FAB>& y, int ycomp,
                         int dstcomp, int numcomp, const IntVect& nghost);
};


#include <AMReX_FabArrayCommI.H>

template <class FAB>
bool
FabArray<FAB>::defined (int K) const noexcept
{
    int li = localindex(K);
    if (li >= 0 && li < static_cast<int>(m_fabs_v.size()) && m_fabs_v[li] != 0) {
        return true;
    }
    else {
        return false;
    }
}

template <class FAB>
bool
FabArray<FAB>::defined (const MFIter& mfi) const noexcept
{
    int li = mfi.LocalIndex();
    if (li < static_cast<int>(m_fabs_v.size()) && m_fabs_v[li] != nullptr) {
        return true;
    }
    else {
        return false;
    }
}

template <class FAB>
FAB*
FabArray<FAB>::fabPtr (const MFIter& mfi) noexcept
{
    BL_ASSERT(mfi.LocalIndex() < indexArray.size());
    BL_ASSERT(DistributionMap() == mfi.DistributionMap());
    int li = mfi.LocalIndex();
    return m_fabs_v[li];
}

template <class FAB>
FAB const*
FabArray<FAB>::fabPtr (const MFIter& mfi) const noexcept
{
    BL_ASSERT(mfi.LocalIndex() < indexArray.size());
    BL_ASSERT(DistributionMap() == mfi.DistributionMap());
    int li = mfi.LocalIndex();
    return m_fabs_v[li];
}

template <class FAB>
FAB*
FabArray<FAB>::fabPtr (int K) noexcept
{
    int li = localindex(K);
    BL_ASSERT(li >=0 && li < indexArray.size());
    return m_fabs_v[li];
}

template <class FAB>
FAB const*
FabArray<FAB>::fabPtr (int K) const noexcept
{
    int li = localindex(K);
    BL_ASSERT(li >=0 && li < indexArray.size());
    return m_fabs_v[li];
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::prefetchToHost (const MFIter& mfi) const noexcept
{
#ifdef AMREX_USE_CUDA
    this->fabPtr(mfi)->prefetchToHost();
#else
    amrex::ignore_unused(mfi);
#endif
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::prefetchToDevice (const MFIter& mfi) const noexcept
{
#ifdef AMREX_USE_CUDA
    this->fabPtr(mfi)->prefetchToDevice();
#else
    amrex::ignore_unused(mfi);
#endif
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::array (const MFIter& mfi) const noexcept
{
    return fabPtr(mfi)->const_array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type>
FabArray<FAB>::array (const MFIter& mfi) noexcept
{
    return fabPtr(mfi)->array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::array (int K) const noexcept
{
    return fabPtr(K)->const_array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type>
FabArray<FAB>::array (int K) noexcept
{
    return fabPtr(K)->array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::const_array (const MFIter& mfi) const noexcept
{
    return fabPtr(mfi)->const_array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::const_array (int K) const noexcept
{
    return fabPtr(K)->const_array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::array (const MFIter& mfi, int start_comp) const noexcept
{
    return fabPtr(mfi)->const_array(start_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type>
FabArray<FAB>::array (const MFIter& mfi, int start_comp) noexcept
{
    return fabPtr(mfi)->array(start_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::array (int K, int start_comp) const noexcept
{
    return fabPtr(K)->const_array(start_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type>
FabArray<FAB>::array (int K, int start_comp) noexcept
{
    return fabPtr(K)->array(start_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::const_array (const MFIter& mfi, int start_comp) const noexcept
{
    return fabPtr(mfi)->const_array(start_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::const_array (int K, int start_comp) const noexcept
{
    return fabPtr(K)->const_array(start_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
MultiArray4<typename FabArray<FAB>::value_type>
FabArray<FAB>::arrays () noexcept
{
    build_arrays();
    return m_arrays;
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
MultiArray4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::arrays () const noexcept
{
    build_arrays();
    return m_const_arrays;
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
MultiArray4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::const_arrays () const noexcept
{
    build_arrays();
    return m_const_arrays;
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::build_arrays () const
{
    using A = Array4<value_type>;
    using AC = Array4<value_type const>;
    static_assert(sizeof(A) == sizeof(AC), "sizeof(Array4<T>) != sizeof(Array4<T const>)");
    if (!m_hp_arrays && local_size() > 0) {
        const int n = local_size();
#ifdef AMREX_USE_GPU
        m_hp_arrays = (void*)The_Pinned_Arena()->alloc(n*2*sizeof(A));
        m_dp_arrays = (void*)The_Arena()->alloc(n*2*sizeof(A));
#else
        m_hp_arrays = (void*)std::malloc(n*2*sizeof(A));
#endif
        for (int li = 0; li < n; ++li) {
            new ((A*)m_hp_arrays+li) A(m_fabs_v[li]->array());
            new ((AC*)m_hp_arrays+li+n) AC(m_fabs_v[li]->const_array());
        }
        m_arrays.hp = (A*)m_hp_arrays;
        m_const_arrays.hp = (AC*)m_hp_arrays + n;
#ifdef AMREX_USE_GPU
        m_arrays.dp = (A*)m_dp_arrays;
        m_const_arrays.dp = (AC*)m_dp_arrays + n;
        Gpu::htod_memcpy(m_dp_arrays, m_hp_arrays, n*2*sizeof(A));
#endif
    }
}

template <class FAB>
void
FabArray<FAB>::clear_arrays ()
{
#ifdef AMREX_USE_GPU
    The_Pinned_Arena()->free(m_hp_arrays);
    The_Arena()->free(m_dp_arrays);
    m_dp_arrays = nullptr;
#else
    std::free(m_hp_arrays);
#endif
    m_hp_arrays = nullptr;
    m_arrays.hp = nullptr;
    m_const_arrays.hp = nullptr;
}

template <class FAB>
AMREX_NODISCARD
FAB*
FabArray<FAB>::release (int K)
{
    const int li = localindex(K);
    if (li >= 0 && li < static_cast<int>(m_fabs_v.size()) && m_fabs_v[li] != nullptr) {
        Long nbytes = amrex::nBytesOwned(*m_fabs_v[li]);
        if (nbytes > 0) {
            for (auto const& t : m_tags) {
                updateMemUsage(t, -nbytes, nullptr);
            }
        }
        return std::exchange(m_fabs_v[li], nullptr);
    } else {
        return nullptr;
    }
}

template <class FAB>
AMREX_NODISCARD
FAB*
FabArray<FAB>::release (const MFIter& mfi)
{
    const int li = mfi.LocalIndex();
    if (li >= 0 && li < static_cast<int>(m_fabs_v.size()) && m_fabs_v[li] != nullptr) {
        Long nbytes = amrex::nBytesOwned(*m_fabs_v[li]);
        if (nbytes > 0) {
            for (auto const& t : m_tags) {
                updateMemUsage(t, -nbytes, nullptr);
            }
        }
        return std::exchange(m_fabs_v[li], nullptr);
    } else {
        return nullptr;
    }
}

template <class FAB>
void
FabArray<FAB>::clear ()
{
    if (define_function_called)
    {
        define_function_called = false;
        clearThisBD();  //!< addThisBD is called in define
    }

    Long nbytes = 0L;
    for (auto *x : m_fabs_v) {
        if (x) {
            nbytes += amrex::nBytesOwned(*x);
            m_factory->destroy(x);
        }
    }
    m_fabs_v.clear();
    clear_arrays();
    m_factory.reset();
    m_dallocator.m_arena = nullptr;
    // no need to clear the non-blocking fillboundary stuff

    if (nbytes > 0) {
        for (auto const& t : m_tags) {
            updateMemUsage(t, -nbytes, nullptr);
        }
    }
    m_tags.clear();

    FabArrayBase::clear();
}

template <class FAB>
template <typename SFAB, typename DFAB,
          std::enable_if_t<std::conjunction_v<
              IsBaseFab<DFAB>, IsBaseFab<SFAB>,
              std::is_convertible<typename SFAB::value_type,
                                  typename DFAB::value_type>>, int>>
void
FabArray<FAB>::LocalCopy (FabArray<SFAB> const& src, int scomp, int dcomp, int ncomp,
                          IntVect const& nghost)
{
    amrex::Copy(*this, src, scomp, dcomp, ncomp, nghost);
}

template <class FAB>
template <class F, std::enable_if_t<IsBaseFab<F>::value,int>>
void
FabArray<FAB>::LocalAdd (FabArray<FAB> const& src, int scomp, int dcomp, int ncomp,
                         IntVect const& nghost)
{
    amrex::Add(*this, src, scomp, dcomp, ncomp, nghost);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val, int nghost)
{
    setVal(val,0,n_comp,IntVect(nghost));
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val, const IntVect& nghost)
{
    setVal(val,0,n_comp,nghost);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val, const Box& region, int nghost)
{
    setVal(val,region,0,n_comp,IntVect(nghost));
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val, const Box& region, const IntVect& nghost)
{
    setVal(val,region,0,n_comp,nghost);
}

template <class FAB>
FabArray<FAB>::FabArray () noexcept
    : shmem()
{
    m_FA_stats.recordBuild();
}

template <class FAB>
FabArray<FAB>::FabArray (Arena* a) noexcept
    : m_dallocator(a),
      shmem()
{
    m_FA_stats.recordBuild();
}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray&            bxs,
                         const DistributionMapping& dm,
                         int                        nvar,
                         int                        ngrow,
                         const MFInfo&              info,
                         const FabFactory<FAB>&     factory)
    : FabArray<FAB>(bxs,dm,nvar,IntVect(ngrow),info,factory)
{}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray&            bxs,
                         const DistributionMapping& dm,
                         int                        nvar,
                         const IntVect&             ngrow,
                         const MFInfo&              info,
                         const FabFactory<FAB>&     factory)
    : m_factory(factory.clone()),
      shmem()
{
    m_FA_stats.recordBuild();
    define(bxs,dm,nvar,ngrow,info,*m_factory);
}

template <class FAB>
FabArray<FAB>::FabArray (const FabArray<FAB>& rhs, MakeType maketype, int scomp, int ncomp)
    : m_factory(rhs.Factory().clone()),
      shmem()
{
    m_FA_stats.recordBuild();
    define(rhs.boxArray(), rhs.DistributionMap(), ncomp, rhs.nGrowVect(),
           MFInfo().SetAlloc(false), *m_factory);

    if (maketype == amrex::make_alias)
    {
        for (int i = 0, n = indexArray.size(); i < n; ++i) {
            auto const& rhsfab = *(rhs.m_fabs_v[i]);
            m_fabs_v.push_back(m_factory->create_alias(rhsfab, scomp, ncomp));
        }
    }
    else
    {
        amrex::Abort("FabArray: unknown MakeType");
    }
}

template <class FAB>
FabArray<FAB>::FabArray (FabArray<FAB>&& rhs) noexcept
    : FabArrayBase (static_cast<FabArrayBase&&>(rhs))
    , m_factory    (std::move(rhs.m_factory))
    , m_dallocator (std::move(rhs.m_dallocator))
    , define_function_called(rhs.define_function_called)
    , m_fabs_v     (std::move(rhs.m_fabs_v))
#ifdef AMREX_USE_GPU
    , m_dp_arrays  (std::exchange(rhs.m_dp_arrays, nullptr))
#endif
    , m_hp_arrays  (std::exchange(rhs.m_hp_arrays, nullptr))
    , m_arrays     (rhs.m_arrays)
    , m_const_arrays(rhs.m_const_arrays)
    , m_tags       (std::move(rhs.m_tags))
    , shmem        (std::move(rhs.shmem))
    // no need to worry about the data used in non-blocking FillBoundary.
{
    m_FA_stats.recordBuild();
    rhs.define_function_called = false; // the responsibility of clear BD has been transferred.
    rhs.m_fabs_v.clear(); // clear the data pointers so that rhs.clear does delete them.
    rhs.clear();
}

template <class FAB>
FabArray<FAB>&
FabArray<FAB>::operator= (FabArray<FAB>&& rhs) noexcept
{
    if (&rhs != this)
    {
        clear();

        FabArrayBase::operator=(static_cast<FabArrayBase&&>(rhs));
        m_factory = std::move(rhs.m_factory);
        m_dallocator = std::move(rhs.m_dallocator);
        define_function_called = rhs.define_function_called;
        std::swap(m_fabs_v, rhs.m_fabs_v);
#ifdef AMREX_USE_GPU
        std::swap(m_dp_arrays, rhs.m_dp_arrays);
#endif
        std::swap(m_hp_arrays, rhs.m_hp_arrays);
        m_arrays = rhs.m_arrays;
        m_const_arrays = rhs.m_const_arrays;
        std::swap(m_tags, rhs.m_tags);
        shmem = std::move(rhs.shmem);

        rhs.define_function_called = false;
        rhs.m_fabs_v.clear();
        rhs.m_tags.clear();
        rhs.clear();
    }
    return *this;
}

template <class FAB>
FabArray<FAB>::~FabArray ()
{
    m_FA_stats.recordDelete();
    clear();
}

template <class FAB>
bool
FabArray<FAB>::ok () const
{
    if (!define_function_called) return false;

    int isok = 1;

    for (MFIter fai(*this); fai.isValid() && isok; ++fai)
    {
        if (defined(fai))
        {
            if (get(fai).box() != fabbox(fai.index()))
            {
                isok = 0;
            }
        }
        else
        {
            isok = 0;
        }
    }

    ParallelAllReduce::Min(isok, ParallelContext::CommunicatorSub());

    return isok == 1;
}

template <class FAB>
bool
FabArray<FAB>::isDefined () const
{
    return define_function_called;
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray&            bxs,
                       const DistributionMapping& dm,
                       int                        nvar,
                       int                        ngrow,
                       const MFInfo&              info,
                       const FabFactory<FAB>&     a_factory)
{
    define(bxs,dm,nvar,IntVect(ngrow),info,a_factory);
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray&            bxs,
                       const DistributionMapping& dm,
                       int                        nvar,
                       const IntVect&             ngrow,
                       const MFInfo&              info,
                       const FabFactory<FAB>&     a_factory)
{
    std::unique_ptr<FabFactory<FAB> > factory(a_factory.clone());

    auto *default_arena = m_dallocator.m_arena;
    clear();

    m_factory = std::move(factory);
    m_dallocator.m_arena = info.arena ? info.arena : default_arena;

    define_function_called = true;

    BL_ASSERT(ngrow.allGE(IntVect::TheZeroVector()));
    BL_ASSERT(boxarray.empty());
    FabArrayBase::define(bxs, dm, nvar, ngrow);

    addThisBD();

    if(info.alloc) {
        AllocFabs(*m_factory, m_dallocator.m_arena, info.tags);
#ifdef BL_USE_TEAM
        ParallelDescriptor::MyTeam().MemoryBarrier();
#endif
    }
}

template <class FAB>
void
FabArray<FAB>::AllocFabs (const FabFactory<FAB>& factory, Arena* ar,
                          const Vector<std::string>& tags)
{
    const int n = indexArray.size();
    const int nworkers = ParallelDescriptor::TeamSize();
    shmem.alloc = (nworkers > 1);

    bool alloc = !shmem.alloc;

    FabInfo fab_info;
    fab_info.SetAlloc(alloc).SetShared(shmem.alloc).SetArena(ar);

    m_fabs_v.reserve(n);

    Long nbytes = 0L;
    for (int i = 0; i < n; ++i)
    {
        int K = indexArray[i];
        const Box& tmpbox = fabbox(K);
        m_fabs_v.push_back(factory.create(tmpbox, n_comp, fab_info, K));
        nbytes += amrex::nBytesOwned(*m_fabs_v.back());
    }

    m_tags.clear();
    m_tags.emplace_back("All");
    for (auto const& t : m_region_tag) {
        m_tags.push_back(t);
    }
    for (auto const& t : tags) {
        m_tags.push_back(t);
    }
    for (auto const& t: m_tags) {
        updateMemUsage(t, nbytes, ar);
    }

#ifdef BL_USE_TEAM
    if (shmem.alloc)
    {
        const int teamlead = ParallelDescriptor::MyTeamLead();

        shmem.n_values = 0;
        shmem.n_points = 0;
        Vector<Long> offset(n,0);
        Vector<Long> nextoffset(nworkers,-1);
        for (int i = 0; i < n; ++i) {
            int K = indexArray[i];
            int owner = distributionMap[K] - teamlead;
            Long s = m_fabs_v[i]->size();
            if (ownership[i]) {
                shmem.n_values += s;
                shmem.n_points += m_fabs_v[i]->numPts();
            }
            if (nextoffset[owner] < 0) {
                offset[i] = 0;
                nextoffset[owner] = s;
            } else {
                offset[i] = nextoffset[owner];
                nextoffset[owner] += s;
            }
        }

        size_t bytes = shmem.n_values*sizeof(value_type);

        value_type *mfp;
        Vector<value_type*> dps;

#if defined (BL_USE_MPI3)

        static MPI_Info info = MPI_INFO_NULL;
        if (info == MPI_INFO_NULL) {
            MPI_Info_create(&info);
            MPI_Info_set(info, "alloc_shared_noncontig", "true");
        }

        const MPI_Comm& team_comm = ParallelDescriptor::MyTeam().get();

        BL_MPI_REQUIRE( MPI_Win_allocate_shared(bytes, sizeof(value_type),
                                                info, team_comm, &mfp, &shmem.win) );

        for (int w = 0; w < nworkers; ++w) {
            MPI_Aint sz;
            int disp;
            value_type *dptr = 0;
            BL_MPI_REQUIRE( MPI_Win_shared_query(shmem.win, w, &sz, &disp, &dptr) );
            // BL_ASSERT(disp == sizeof(value_type));
            dps.push_back(dptr);
        }

#else

        amrex::Abort("BaseFab::define: to allocate shared memory, USE_MPI3 must be true");

#endif

        for (int i = 0; i < n; ++i) {
            int K = indexArray[i];
            int owner = distributionMap[K] - teamlead;
            value_type *p = dps[owner] + offset[i];
            m_fabs_v[i]->setPtr(p, m_fabs_v[i]->size());
        }

        for (Long i = 0; i < shmem.n_values; i++, mfp++) {
            new (mfp) value_type;
        }

        amrex::update_fab_stats(shmem.n_points, shmem.n_values, sizeof(value_type));
    }
#endif
}

template <class FAB>
void
FabArray<FAB>::setFab_assert (int K, FAB const& fab) const
{
    amrex::ignore_unused(K,fab);
    AMREX_ASSERT(n_comp == fab.nComp());
    AMREX_ASSERT(!boxarray.empty());
    AMREX_ASSERT(fab.box() == fabbox(K));
    AMREX_ASSERT(distributionMap[K] == ParallelDescriptor::MyProc());
}

template <class FAB>
void
FabArray<FAB>::setFab (int boxno, std::unique_ptr<FAB> elem)
{
    if (n_comp == 0) {
        n_comp = elem->nComp();
    }

    setFab_assert(boxno, *elem);

    if (m_fabs_v.empty()) {
        m_fabs_v.resize(indexArray.size(),nullptr);
    }

    const int li = localindex(boxno);
    if (m_fabs_v[li]) {
        m_factory->destroy(m_fabs_v[li]);
    }
    m_fabs_v[li] = elem.release();
}

template <class FAB>
template <class F, std::enable_if_t<std::is_move_constructible<F>::value,int> >
void
FabArray<FAB>::setFab (int boxno, FAB&& elem)
{
    if (n_comp == 0) {
        n_comp = elem.nComp();
    }

    setFab_assert(boxno, elem);

    if (m_fabs_v.empty()) {
        m_fabs_v.resize(indexArray.size(),nullptr);
    }

    const int li = localindex(boxno);
    if (m_fabs_v[li]) {
        m_factory->destroy(m_fabs_v[li]);
    }
    m_fabs_v[li] = new FAB(std::move(elem));
}

template <class FAB>
void
FabArray<FAB>::setFab (const MFIter& mfi, std::unique_ptr<FAB> elem)
{
    if (n_comp == 0) {
        n_comp = elem->nComp();
    }

    setFab_assert(mfi.index(), *elem);

    if (m_fabs_v.empty()) {
        m_fabs_v.resize(indexArray.size(),nullptr);
    }

    const int li = mfi.LocalIndex();
    if (m_fabs_v[li]) {
        m_factory->destroy(m_fabs_v[li]);
    }
    m_fabs_v[li] = elem.release();
}

template <class FAB>
template <class F, std::enable_if_t<std::is_move_constructible<F>::value,int> >
void
FabArray<FAB>::setFab (const MFIter& mfi, FAB&& elem)
{
    if (n_comp == 0) {
        n_comp = elem.nComp();
    }

    setFab_assert(mfi.index(), elem);

    if (m_fabs_v.empty()) {
        m_fabs_v.resize(indexArray.size(),nullptr);
    }

    const int li = mfi.LocalIndex();
    if (m_fabs_v[li]) {
        m_factory->destroy(m_fabs_v[li]);
    }
    m_fabs_v[li] = new FAB(std::move(elem));
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setBndry (value_type val)
{
    setBndry(val, 0, n_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::setBndry (value_type val,
                         int        strt_comp,
                         int        ncomp)
{
    if (n_grow.max() > 0)
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            bool use_mfparfor = true;
            const int nboxes = local_size();
            if (nboxes == 1) {
                if (boxarray[indexArray[0]].numPts() > Long(65*65*65)) {
                    use_mfparfor = false;
                }
            } else {
                for (int i = 0; i < nboxes; ++i) {
                    const Long npts = boxarray[indexArray[i]].numPts();
                    if (npts >= Long(64*64*64)) {
                        use_mfparfor = false;
                        break;
                    } else if (npts <= Long(17*17*17)) {
                        break;
                    }
                }
            }
            const IntVect nghost = n_grow;
            if (use_mfparfor) {
                auto const& ma = this->arrays();
                ParallelFor(*this, nghost,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    auto const& a = ma[box_no];
                    Box vbx(a);
                    vbx.grow(-nghost);
                    if (!vbx.contains(i,j,k)) {
                        for (int n = 0; n < ncomp; ++n) {
                            a(i,j,k,strt_comp+n) = val;
                        }
                    }
                });
                Gpu::streamSynchronize();
            } else {
                using Tag = Array4BoxTag<value_type>;
                Vector<Tag> tags;
                for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
                    Box const& vbx = mfi.validbox();
                    auto const& a = this->array(mfi);

                    Box b;
#if (AMREX_SPACEDIM == 3)
                    if (nghost[2] > 0) {
                        b = vbx;
                        b.setRange(2, vbx.smallEnd(2)-nghost[2], nghost[2]);
                        b.grow(IntVect(nghost[0],nghost[1],0));
                        tags.emplace_back(Tag{a, b});
                        b.shift(2, vbx.length(2)+nghost[2]);
                        tags.emplace_back(Tag{a, b});
                    }
#endif
#if (AMREX_SPACEDIM >= 2)
                    if (nghost[1] > 0) {
                        b = vbx;
                        b.setRange(1, vbx.smallEnd(1)-nghost[1], nghost[1]);
                        b.grow(0, nghost[0]);
                        tags.emplace_back(Tag{a, b});
                        b.shift(1, vbx.length(1)+nghost[1]);
                        tags.emplace_back(Tag{a, b});
                    }
#endif
                    if (nghost[0] > 0) {
                        b = vbx;
                        b.setRange(0, vbx.smallEnd(0)-nghost[0], nghost[0]);
                        tags.emplace_back(Tag{a, b});
                        b.shift(0, vbx.length(0)+nghost[0]);
                        tags.emplace_back(Tag{a, b});
                    }
                }

                ParallelFor(tags, ncomp,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n, Tag const& tag) noexcept
                {
                    tag.dfab(i,j,k,strt_comp+n) = val;
                });
            }
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            for (MFIter fai(*this); fai.isValid(); ++fai)
            {
                get(fai).template setComplement<RunOn::Host>(val, fai.validbox(), strt_comp, ncomp);
            }
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setDomainBndry (value_type val, const Geometry& geom)
{
    setDomainBndry(val, 0, n_comp, geom);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setDomainBndry (value_type val,
                               int        strt_comp,
                               int        ncomp,
                               const Geometry& geom)
{
    BL_PROFILE("FabArray::setDomainBndry()");

    Box domain_box = amrex::convert(geom.Domain(), boxArray().ixType());
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            int n = domain_box.length(idim);
            domain_box.grow(idim, n);
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        const Box& gbx = fai.fabbox();
        if (! domain_box.contains(gbx))
        {
            get(fai).template setComplement<RunOn::Device>(val, domain_box, strt_comp, ncomp);
        }
    }
}

template <class FAB>
template <class F, std::enable_if_t<IsBaseFab<F>::value,int> FOO>
typename F::value_type
FabArray<FAB>::sum (int comp, IntVect const& nghost, bool local) const
{
    BL_PROFILE("FabArray::sum()");

    using T = typename FAB::value_type;
    auto sm = T(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& ma = this->const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<T>{}, *this, nghost,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                       -> GpuTuple<T>
        {
            return ma[box_no](i,j,k,comp);
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(*this,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.growntilebox(nghost);
            auto const& a = this->const_array(mfi);
            auto tmp = T(0.0);
            AMREX_LOOP_3D(bx, i, j, k,
            {
                tmp += a(i,j,k,comp);
            });
            sm += tmp; // Do it this way so that it does not break regression tests.
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    return sm;
}

template <class FAB>
void
FabArray<FAB>::copyTo (FAB& dest, int  nghost) const
{
    copyTo(dest, 0, 0, dest.nComp(), nghost);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val)
{
    setVal(val,0,n_comp,n_grow);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
FabArray<FAB>&
FabArray<FAB>::operator= (value_type val)
{
    setVal(val);
    return *this;
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val,
                       int        comp,
                       int        ncomp,
                       int        nghost)
{
    setVal(val,comp,ncomp,IntVect(nghost));
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::setVal (value_type val,
                       int        comp,
                       int        ncomp,
                       const IntVect& nghost)
{
    BL_ASSERT(nghost.allGE(IntVect::TheZeroVector()) && nghost.allLE(n_grow));
    BL_ASSERT(comp+ncomp <= n_comp);

    BL_PROFILE("FabArray::setVal()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && this->isFusingCandidate()) {
        auto const& fa = this->arrays();
        ParallelFor(*this, nghost, ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            fa[box_no](i,j,k,n+comp) = val;
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter fai(*this,TilingIfNotGPU()); fai.isValid(); ++fai)
        {
            const Box& bx = fai.growntilebox(nghost);
            auto fab = this->array(fai);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, ncomp, i, j, k, n,
            {
                fab(i,j,k,n+comp) = val;
            });
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val,
                       const Box& region,
                       int        comp,
                       int        ncomp,
                       int        nghost)
{
    setVal(val,region,comp,ncomp,IntVect(nghost));
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::setVal (value_type val,
                       const Box& region,
                       int        comp,
                       int        ncomp,
                       const IntVect& nghost)
{
    BL_ASSERT(nghost.allGE(IntVect::TheZeroVector()) && nghost.allLE(n_grow));
    BL_ASSERT(comp+ncomp <= n_comp);

    BL_PROFILE("FabArray::setVal(val,region,comp,ncomp,nghost)");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && this->isFusingCandidate()) {
        auto const& fa = this->arrays();
        ParallelFor(*this, nghost, ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (region.contains(i,j,k)) {
                fa[box_no](i,j,k,n+comp) = val;
            }
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
        AMREX_ALWAYS_ASSERT(!omp_in_parallel());
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter fai(*this,TilingIfNotGPU()); fai.isValid(); ++fai)
        {
            Box b = fai.growntilebox(nghost) & region;

            if (b.ok()) {
                auto fab = this->array(fai);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( b, ncomp, i, j, k, n,
                {
                    fab(i,j,k,n+comp) = val;
                });
            }
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::abs (int comp, int ncomp, int nghost)
{
    abs(comp, ncomp, IntVect(nghost));
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::abs (int comp, int ncomp, const IntVect& nghost)
{
    BL_ASSERT(nghost.allGE(IntVect::TheZeroVector()) && nghost.allLE(n_grow));
    BL_ASSERT(comp+ncomp <= n_comp);
    BL_PROFILE("FabArray::abs()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && this->isFusingCandidate()) {
        auto const& fa = this->arrays();
        ParallelFor(*this, nghost, ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            fa[box_no](i,j,k,n+comp) = std::abs(fa[box_no](i,j,k,n+comp));
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            auto fab = this->array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, ncomp, i, j, k, n,
            {
                fab(i,j,k,n+comp) = std::abs(fab(i,j,k,n+comp));
            });
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::plus (value_type val, int comp, int num_comp, int nghost)
{
    BL_PROFILE("FabArray::plus()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && this->isFusingCandidate()) {
        auto const& fa = this->arrays();
        ParallelFor(*this, IntVect(nghost), num_comp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            fa[box_no](i,j,k,n+comp) += val;
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            auto fab = this->array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, num_comp, i, j, k, n,
            {
                fab(i,j,k,n+comp) += val;
            });
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::plus (value_type val, const Box& region, int comp, int num_comp, int nghost)
{
    BL_PROFILE("FabArray::plus(val, region, comp, num_comp, nghost)");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && this->isFusingCandidate()) {
        auto const& fa = this->arrays();
        ParallelFor(*this, IntVect(nghost), num_comp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (region.contains(i,j,k)) {
                fa[box_no](i,j,k,n+comp) += val;
            }
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost) & region;
            if (bx.ok()) {
                auto fab = this->array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, num_comp, i, j, k, n,
                {
                    fab(i,j,k,n+comp) += val;
                });
            }
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::mult (value_type val, int comp, int num_comp, int nghost)
{
    BL_PROFILE("FabArray::mult()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && this->isFusingCandidate()) {
        auto const& fa = this->arrays();
        ParallelFor(*this, IntVect(nghost), num_comp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            fa[box_no](i,j,k,n+comp) *= val;
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            auto fab = this->array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, num_comp, i, j, k, n,
            {
                fab(i,j,k,n+comp) *= val;
            });
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::mult (value_type val, const Box& region, int comp, int num_comp, int nghost)
{
    BL_PROFILE("FabArray::mult(val, region, comp, num_comp, nghost)");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && this->isFusingCandidate()) {
        auto const& fa = this->arrays();
        ParallelFor(*this, IntVect(nghost), num_comp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (region.contains(i,j,k)) {
                fa[box_no](i,j,k,n+comp) *= val;
            }
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost) & region;
            if (bx.ok()) {
                auto fab = this->array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, num_comp, i, j, k, n,
                {
                    fab(i,j,k,n+comp) *= val;
                });
            }
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::invert (value_type numerator, int comp, int num_comp, int nghost)
{
    BL_PROFILE("FabArray::invert()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && this->isFusingCandidate()) {
        auto const& fa = this->arrays();
        ParallelFor(*this, IntVect(nghost), num_comp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            fa[box_no](i,j,k,n+comp) = numerator / fa[box_no](i,j,k,n+comp);
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            auto fab = this->array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, num_comp, i, j, k, n,
            {
                fab(i,j,k,n+comp) = numerator / fab(i,j,k,n+comp);
            });
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::invert (value_type numerator, const Box& region, int comp, int num_comp, int nghost)
{
    BL_PROFILE("FabArray::invert(numerator, region, comp, num_comp, nghost)");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && this->isFusingCandidate()) {
        auto const& fa = this->arrays();
        ParallelFor(*this, IntVect(nghost), num_comp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (region.contains(i,j,k)) {
                fa[box_no](i,j,k,n+comp) = numerator / fa[box_no](i,j,k,n+comp);
            }
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost) & region;
            if (bx.ok()) {
                auto fab = this->array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, num_comp, i, j, k, n,
                {
                    fab(i,j,k,n+comp) = numerator / fab(i,j,k,n+comp);
                });
            }
        }
    }
}

template <class FAB>
void
FabArray<FAB>::shift (const IntVect& v)
{
    clearThisBD();  // The new boxarry will have a different ID.
    boxarray.shift(v);
    addThisBD();
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        get(fai).shift(v);
    }
    clear_arrays();
}

template <class FAB>
template <class F, std::enable_if_t<IsBaseFab<F>::value,int> FOO>
void FabArray<FAB>::Saxpy (FabArray<FAB>& y, value_type a, FabArray<FAB> const& x,
                           int xcomp, int ycomp, int ncomp, IntVect const& nghost)
{
    BL_ASSERT(y.boxArray() == x.boxArray());
    BL_ASSERT(y.distributionMap == x.distributionMap);
    BL_ASSERT(y.nGrowVect().allGE(nghost) && x.nGrowVect().allGE(nghost));

    BL_PROFILE("FabArray::Saxpy()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && y.isFusingCandidate()) {
        auto const& yma = y.arrays();
        auto const& xma = x.const_arrays();
        ParallelFor(y, nghost, ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            yma[box_no](i,j,k,ycomp+n) += a * xma[box_no](i,j,k,xcomp+n);
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(y,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);

            if (bx.ok()) {
                auto const& xfab = x.const_array(mfi);
                auto const& yfab = y.array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, ncomp, i, j, k, n,
                {
                    yfab(i,j,k,ycomp+n) += a * xfab(i,j,k,xcomp+n);
                });
            }
        }
    }
}

template <class FAB>
template <class F, std::enable_if_t<IsBaseFab<F>::value,int> FOO>
void
FabArray<FAB>::Xpay (FabArray<FAB>& y, value_type a, FabArray<FAB> const& x,
                     int xcomp, int ycomp, int ncomp, IntVect const& nghost)
{
    BL_ASSERT(y.boxArray() == x.boxArray());
    BL_ASSERT(y.distributionMap == x.distributionMap);
    BL_ASSERT(y.nGrowVect().allGE(nghost) && x.nGrowVect().allGE(nghost));

    BL_PROFILE("FabArray::Xpay()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && y.isFusingCandidate()) {
        auto const& yfa = y.arrays();
        auto const& xfa = x.const_arrays();
        ParallelFor(y, nghost, ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            yfa[box_no](i,j,k,n+ycomp) = xfa[box_no](i,j,k,n+xcomp)
                +                    a * yfa[box_no](i,j,k,n+ycomp);
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(y,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            auto const& xFab = x.const_array(mfi);
            auto const& yFab = y.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, ncomp, i, j, k, n,
            {
                yFab(i,j,k,n+ycomp) = xFab(i,j,k,n+xcomp)
                    +             a * yFab(i,j,k,n+ycomp);
            });
        }
    }
}

template <class FAB>
template <class F, std::enable_if_t<IsBaseFab<F>::value,int> FOO>
void
FabArray<FAB>::LinComb (FabArray<FAB>& dst,
                        value_type a, const FabArray<FAB>& x, int xcomp,
                        value_type b, const FabArray<FAB>& y, int ycomp,
                        int dstcomp, int numcomp, const IntVect& nghost)
{
    BL_ASSERT(dst.boxArray() == x.boxArray());
    BL_ASSERT(dst.distributionMap == x.distributionMap);
    BL_ASSERT(dst.boxArray() == y.boxArray());
    BL_ASSERT(dst.distributionMap == y.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost) && x.nGrowVect().allGE(nghost) && y.nGrowVect().allGE(nghost));

    BL_PROFILE("FabArray::LinComb()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && dst.isFusingCandidate()) {
        auto const& dstma = dst.arrays();
        auto const& xma = x.const_arrays();
        auto const& yma = y.const_arrays();
        ParallelFor(dst, nghost, numcomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            dstma[box_no](i,j,k,dstcomp+n) = a*xma[box_no](i,j,k,xcomp+n)
                +                            b*yma[box_no](i,j,k,ycomp+n);
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            auto const& xfab = x.const_array(mfi);
            auto const& yfab = y.const_array(mfi);
            auto const& dfab = dst.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, numcomp, i, j, k, n,
            {
                dfab(i,j,k,dstcomp+n) = a*xfab(i,j,k,xcomp+n) + b*yfab(i,j,k,ycomp+n);
            });
        }
    }
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary (bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
        FillBoundary_nowait<BUF>(0, nComp(), n_grow, Periodicity::NonPeriodic(), cross);
        FillBoundary_finish<BUF>();
    }
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary (const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
        FillBoundary_nowait<BUF>(0, nComp(), n_grow, period, cross);
        FillBoundary_finish<BUF>();
    }
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary (const IntVect& nghost, const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nghost.allLE(nGrowVect()),
                                     "FillBoundary: asked to fill more ghost cells than we have");
    if ( nghost.max() > 0 ) {
        FillBoundary_nowait<BUF>(0, nComp(), nghost, period, cross);
        FillBoundary_finish<BUF>();
    }
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary (int scomp, int ncomp, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
        FillBoundary_nowait<BUF>(scomp, ncomp, n_grow, Periodicity::NonPeriodic(), cross);
        FillBoundary_finish<BUF>();
    }
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary (int scomp, int ncomp, const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
        FillBoundary_nowait<BUF>(scomp, ncomp, n_grow, period, cross);
        FillBoundary_finish<BUF>();
    }
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary (int scomp, int ncomp, const IntVect& nghost,
                             const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nghost.allLE(nGrowVect()),
                                     "FillBoundary: asked to fill more ghost cells than we have");
    if ( nghost.max() > 0 ) {
        FillBoundary_nowait<BUF>(scomp, ncomp, nghost, period, cross);
        FillBoundary_finish<BUF>();
    }
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary_nowait (bool cross)
{
    FillBoundary_nowait<BUF>(0, nComp(), nGrowVect(), Periodicity::NonPeriodic(), cross);
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary_nowait (const Periodicity& period, bool cross)
{
    FillBoundary_nowait<BUF>(0, nComp(), nGrowVect(), period, cross);
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary_nowait (const IntVect& nghost, const Periodicity& period, bool cross)
{
    FillBoundary_nowait<BUF>(0, nComp(), nghost, period, cross);
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary_nowait (int scomp, int ncomp, bool cross)
{
    FillBoundary_nowait<BUF>(scomp, ncomp, nGrowVect(), Periodicity::NonPeriodic(), cross);
}

template <class FAB>
void
FabArray<FAB>::FillBoundaryAndSync (const Periodicity& period)
{
    BL_PROFILE("FabArray::FillBoundaryAndSync()");
    if (n_grow.max() > 0 || !is_cell_centered()) {
        FillBoundaryAndSync_nowait(0, nComp(), n_grow, period);
        FillBoundaryAndSync_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundaryAndSync (int scomp, int ncomp, const IntVect& nghost,
                                    const Periodicity& period)
{
    BL_PROFILE("FabArray::FillBoundaryAndSync()");
    if (nghost.max() > 0 || !is_cell_centered()) {
        FillBoundaryAndSync_nowait(scomp, ncomp, nghost, period);
        FillBoundaryAndSync_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundaryAndSync_nowait (const Periodicity& period)
{
    FillBoundaryAndSync_nowait(0, nComp(), nGrowVect(), period);
}

template <class FAB>
void
FabArray<FAB>::FillBoundaryAndSync_nowait (int scomp, int ncomp, const IntVect& nghost,
                                           const Periodicity& period)
{
    BL_PROFILE("FillBoundaryAndSync_nowait()");
    FBEP_nowait(scomp, ncomp, nghost, period, false, false, true);
}

template <class FAB>
void
FabArray<FAB>::FillBoundaryAndSync_finish ()
{
    BL_PROFILE("FillBoundaryAndSync_finish()");
    FillBoundary_finish();
}

template <class FAB>
void
FabArray<FAB>::OverrideSync (const Periodicity& period)
{
    BL_PROFILE("FAbArray::OverrideSync()");
    if (!is_cell_centered()) {
        OverrideSync_nowait(0, nComp(), period);
        OverrideSync_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::OverrideSync (int scomp, int ncomp, const Periodicity& period)
{
    BL_PROFILE("FAbArray::OverrideSync()");
    if (!is_cell_centered()) {
        OverrideSync_nowait(scomp, ncomp, period);
        OverrideSync_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::OverrideSync_nowait (const Periodicity& period)
{
    OverrideSync_nowait(0, nComp(), period);
}

template <class FAB>
void
FabArray<FAB>::OverrideSync_nowait (int scomp, int ncomp, const Periodicity& period)
{
    BL_PROFILE("OverrideSync_nowait()");
    FBEP_nowait(scomp, ncomp, IntVect(0), period, false, false, true);
}

template <class FAB>
void
FabArray<FAB>::OverrideSync_finish ()
{
    BL_PROFILE("OverrideSync_finish()");
    FillBoundary_finish();
}

template <class FAB>
void
FabArray<FAB>::SumBoundary (const Periodicity& period)
{
    SumBoundary(0, n_comp, IntVect(0), period);
}

template <class FAB>
void
FabArray<FAB>::SumBoundary (int scomp, int ncomp, const Periodicity& period)
{
    SumBoundary(scomp, ncomp, IntVect(0), period);
}

template <class FAB>
void
FabArray<FAB>::SumBoundary (int scomp, int ncomp, IntVect const& nghost, const Periodicity& period)
{
    SumBoundary(scomp, ncomp, this->nGrowVect(), nghost, period);
}

template <class FAB>
void
FabArray<FAB>::SumBoundary (int scomp, int ncomp, IntVect const& src_nghost, IntVect const& dst_nghost, const Periodicity& period)
{
    BL_PROFILE("FabArray<FAB>::SumBoundary()");

    SumBoundary_nowait(scomp, ncomp, src_nghost, dst_nghost, period);
    SumBoundary_finish();
}

template <class FAB>
void
FabArray<FAB>::SumBoundary_nowait (const Periodicity& period)
{
    SumBoundary_nowait(0, n_comp, IntVect(0), period);
}

template <class FAB>
void
FabArray<FAB>::SumBoundary_nowait (int scomp, int ncomp, const Periodicity& period)
{
    SumBoundary_nowait(scomp, ncomp, IntVect(0), period);
}

template <class FAB>
void
FabArray<FAB>::SumBoundary_nowait (int scomp, int ncomp, IntVect const& nghost, const Periodicity& period)
{
    SumBoundary_nowait(scomp, ncomp, this->nGrowVect(), nghost, period);
}

template <class FAB>
void
FabArray<FAB>::SumBoundary_nowait (int scomp, int ncomp, IntVect const& src_nghost, IntVect const& dst_nghost, const Periodicity& period)
{
    BL_PROFILE("FabArray<FAB>::SumBoundary_nowait()");

    if ( n_grow == IntVect::TheZeroVector() && boxArray().ixType().cellCentered()) return;

    AMREX_ASSERT(src_nghost <= n_grow);

    auto* tmp = new FabArray<FAB>( boxArray(), DistributionMap(), ncomp, src_nghost, MFInfo(), Factory() );
    amrex::Copy(*tmp, *this, scomp, 0, ncomp, src_nghost);
    this->setVal(typename FAB::value_type(0), scomp, ncomp, dst_nghost);
    this->ParallelCopy_nowait(*tmp,0,scomp,ncomp,src_nghost,dst_nghost,period,FabArrayBase::ADD);

    // All local. Operation complete.
    if (!this->pcd) { delete tmp; }
}

template <class FAB>
void
FabArray<FAB>::SumBoundary_finish ()
{
    BL_PROFILE("FabArray<FAB>::SumBoundary_finish()");

    // If pcd doesn't exist, ParallelCopy was all local and operation was fully completed in "SumBoundary_nowait".
    if ( (n_grow == IntVect::TheZeroVector() && boxArray().ixType().cellCentered()) || !(this->pcd) ) return;

    auto* tmp = const_cast<FabArray<FAB>*> (this->pcd->src);
    this->ParallelCopy_finish();
    delete tmp;
}

template <class FAB>
void
FabArray<FAB>::EnforcePeriodicity (const Periodicity& period)
{
    BL_PROFILE("FabArray::EnforcePeriodicity");
    if (period.isAnyPeriodic()) {
        FBEP_nowait(0, nComp(), nGrowVect(), period, false, true);
        FillBoundary_finish(); // unsafe unless isAnyPeriodic()
    }
}

template <class FAB>
void
FabArray<FAB>::EnforcePeriodicity (int scomp, int ncomp, const Periodicity& period)
{
    BL_PROFILE("FabArray::EnforcePeriodicity");
    if (period.isAnyPeriodic()) {
        FBEP_nowait(scomp, ncomp, nGrowVect(), period, false, true);
        FillBoundary_finish(); // unsafe unless isAnyPeriodic()
    }
}

template <class FAB>
void
FabArray<FAB>::EnforcePeriodicity (int scomp, int ncomp, const IntVect& nghost,
                                   const Periodicity& period)
{
    BL_PROFILE("FabArray::EnforcePeriodicity");
    if (period.isAnyPeriodic()) {
        FBEP_nowait(scomp, ncomp, nghost, period, false, true);
        FillBoundary_finish(); // unsafe unless isAnyPeriodic()
    }
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary_nowait (int scomp, int ncomp, const Periodicity& period, bool cross)
{
    FBEP_nowait<BUF>(scomp, ncomp, nGrowVect(), period, cross);
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FillBoundary_nowait (int scomp, int ncomp, const IntVect& nghost,
                                    const Periodicity& period, bool cross)
{
    FBEP_nowait<BUF>(scomp, ncomp, nghost, period, cross);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::BuildMask (const Box& phys_domain, const Periodicity& period,
                          value_type covered, value_type notcovered,
                          value_type physbnd, value_type interior)
{
    BL_PROFILE("FabArray::BuildMask()");

    int ncomp = this->nComp();
    const IntVect& ngrow = this->nGrowVect();

    Box domain = amrex::convert(phys_domain, boxArray().ixType());
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (period.isPeriodic(i)) {
            domain.grow(i, ngrow[i]);
        }
    }

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && this->isFusingCandidate()) {
        auto const& fa = this->arrays();
        ParallelFor(*this, ngrow, ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            auto const& fab = fa[box_no];
            Box vbx(fab);
            vbx.grow(-ngrow);
            if (vbx.contains(i,j,k)) {
                fab(i,j,k,n) = interior;
            } else if (domain.contains(i,j,k)) {
                fab(i,j,k,n) = notcovered;
            } else {
                fab(i,j,k,n) = physbnd;
            }
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            auto const& fab = this->array(mfi);
            Box const& fbx = mfi.growntilebox();
            Box const& gbx = fbx & domain;
            Box const& vbx = mfi.validbox();
            AMREX_HOST_DEVICE_FOR_4D(fbx, ncomp, i, j, k, n,
            {
                if (vbx.contains(i,j,k)) {
                    fab(i,j,k,n) = interior;
                } else if (gbx.contains(i,j,k)) {
                    fab(i,j,k,n) = notcovered;
                } else {
                    fab(i,j,k,n) = physbnd;
                }
            });
        }
    }

    const FabArrayBase::FB& TheFB = this->getFB(ngrow,period);
    setVal(covered, TheFB, 0, ncomp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val, const CommMetaData& thecmd, int scomp, int ncomp)
{
    BL_PROFILE("FabArray::setVal(val, thecmd, scomp, ncomp)");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        CMD_local_setVal_gpu(val, thecmd, scomp, ncomp);
        CMD_remote_setVal_gpu(val, thecmd, scomp, ncomp);
    }
    else
#endif
    {
        const CopyComTagsContainer&      LocTags = *(thecmd.m_LocTags);
        const MapOfCopyComTagContainers& RcvTags = *(thecmd.m_RcvTags);
        auto N_locs = static_cast<int>(LocTags.size());
#ifdef AMREX_USE_OMP
#pragma omp parallel for if (thecmd.m_threadsafe_loc)
#endif
        for (int i = 0; i < N_locs; ++i) {
            const CopyComTag& tag = LocTags[i];
            (*this)[tag.dstIndex].template setVal<RunOn::Host>(val, tag.dbox, scomp, ncomp);
        }

        for (const auto & RcvTag : RcvTags) {
            auto N = static_cast<int>(RcvTag.second.size());
#ifdef AMREX_USE_OMP
#pragma omp parallel for if (thecmd.m_threadsafe_rcv)
#endif
            for (int i = 0; i < N; ++i) {
                const CopyComTag& tag = RcvTag.second[i];
                (*this)[tag.dstIndex].template setVal<RunOn::Host>(val, tag.dbox, scomp, ncomp);
            }
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
LayoutData<int>
FabArray<FAB>::RecvLayoutMask (const CommMetaData& thecmd)
{
    BL_PROFILE("FabArray::RecvLayoutMask()");

    LayoutData<int> r(this->boxArray(), this->DistributionMap());
#ifdef AMREX_USE_OMP
#pragma omp parallel if (thecmd.m_threadsafe_rcv)
#endif
    for (MFIter mfi(r); mfi.isValid(); ++mfi) {
        r[mfi] = 0;
    }

    const CopyComTagsContainer&      LocTags = *(thecmd.m_LocTags);
    const MapOfCopyComTagContainers& RcvTags = *(thecmd.m_RcvTags);

    auto N_locs = static_cast<int>(LocTags.size());
    for (int i = 0; i < N_locs; ++i) {
        const CopyComTag& tag = LocTags[i];
        r[tag.dstIndex] = 1;
    }

    for (const auto & RcvTag : RcvTags) {
        auto N = static_cast<int>(RcvTag.second.size());
        for (int i = 0; i < N; ++i) {
            const CopyComTag& tag = RcvTag.second[i];
            r[tag.dstIndex] = 1;
        }
    }
    return r;
}

template <class FAB>
template <typename F, std::enable_if_t<IsBaseFab<F>::value,int> FOO>
typename F::value_type
FabArray<FAB>::norminf (int comp, int ncomp, IntVect const& nghost, bool local,
                        [[maybe_unused]] bool ignore_covered) const
{
    BL_PROFILE("FabArray::norminf()");

    using RT = typename F::value_type;

    auto nm0 = RT(0.0);

#ifdef AMREX_USE_EB
    if ( this->is_cell_centered() && this->hasEBFabFactory() && ignore_covered )
    {
        const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(this->Factory());
        auto const& flags = ebfactory.getMultiEBCellFlagFab();
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& flagsma = flags.const_arrays();
            auto const& ma = this->const_arrays();
            nm0 = ParReduce(TypeList<ReduceOpMax>{}, TypeList<RT>{}, *this, nghost,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<RT>
            {
                if (flagsma[box_no](i,j,k).isCovered()) {
                    return RT(0.0);
                } else {
                    auto tmp = RT(0.0);
                    auto const& a = ma[box_no];
                    for (int n = 0; n < ncomp; ++n) {
                        tmp = amrex::max(tmp, std::abs(a(i,j,k,comp+n)));
                    }
                    return tmp;
                }
            });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:nm0)
#endif
            for (MFIter mfi(*this,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.growntilebox(nghost);
                if (flags[mfi].getType(bx) != FabType::covered) {
                    auto const& flag = flags.const_array(mfi);
                    auto const& a = this->const_array(mfi);
                    AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
                    {
                        if (!flag(i,j,k).isCovered()) {
                            nm0 = std::max(nm0, std::abs(a(i,j,k,comp+n)));
                        }
                    });
                }
            }
        }
    }
    else
#endif
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& ma = this->const_arrays();
            nm0 = ParReduce(TypeList<ReduceOpMax>{}, TypeList<RT>{}, *this, nghost, ncomp,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept -> GpuTuple<RT>
            {
                return std::abs(ma[box_no](i,j,k,comp+n));
            });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:nm0)
#endif
            for (MFIter mfi(*this,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.growntilebox(nghost);
                auto const& a = this->const_array(mfi);
                AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
                {
                    nm0 = std::max(nm0, std::abs(a(i,j,k,comp+n)));
                });
            }
        }
    }

    if (!local) {
        ParallelAllReduce::Max(nm0, ParallelContext::CommunicatorSub());
    }

    return nm0;
}

template <class FAB>
template <typename IFAB, typename F, std::enable_if_t<IsBaseFab<F>::value,int> FOO>
typename F::value_type
FabArray<FAB>::norminf (FabArray<IFAB> const& mask, int comp, int ncomp,
                        IntVect const& nghost, bool local) const
{
    BL_PROFILE("FabArray::norminf(mask)");

    using RT = typename F::value_type;

    auto nm0 = RT(0.0);

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& ma = this->const_arrays();
        auto const& maskma = mask.const_arrays();
        nm0 = ParReduce(TypeList<ReduceOpMax>{}, TypeList<RT>{}, *this, IntVect(nghost),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<RT>
        {
            if (maskma[box_no](i,j,k)) {
                auto tmp = RT(0.0);
                auto const& a = ma[box_no];
                for (int n = 0; n < ncomp; ++n) {
                    tmp = amrex::max(tmp, std::abs(a(i,j,k,comp+n)));
                }
                return tmp;
            } else {
                return RT(0.0);
            }
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:nm0)
#endif
        for (MFIter mfi(*this,true); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.growntilebox(nghost);
            auto const& a = this->const_array(mfi);
            auto const& mskfab = mask.const_array(mfi);
            AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
            {
                if (mskfab(i,j,k)) {
                    nm0 = std::max(nm0, std::abs(a(i,j,k,comp+n)));
                }
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Max(nm0, ParallelContext::CommunicatorSub());
    }

    return nm0;
}

}

#endif /*BL_FABARRAY_H*/
