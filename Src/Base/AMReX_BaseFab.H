#ifndef AMREX_BASEFAB_H_
#define AMREX_BASEFAB_H_

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <climits>
#include <array>
#include <type_traits>
#include <memory>
#include <atomic>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Extension.H>
#include <AMReX_BLassert.H>
#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_Loop.H>
#include <AMReX_BoxList.H>
#include <AMReX_BArena.H>
#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_MakeType.H>
#include <AMReX_Utility.H>
#include <AMReX_Reduce.H>
//#include <AMReX_Scan.H>
#include <AMReX_Gpu.H>
#include <AMReX_Math.H>

namespace amrex
{

extern std::atomic<Long> atomic_total_bytes_allocated_in_fabs;
extern std::atomic<Long> atomic_total_bytes_allocated_in_fabs_hwm;
extern std::atomic<Long> atomic_total_cells_allocated_in_fabs;
extern std::atomic<Long> atomic_total_cells_allocated_in_fabs_hwm;
extern Long private_total_bytes_allocated_in_fabs;     //!< total bytes at any given time
extern Long private_total_bytes_allocated_in_fabs_hwm; //!< high-water-mark over a given interval
extern Long private_total_cells_allocated_in_fabs;     //!< total cells at any given time
extern Long private_total_cells_allocated_in_fabs_hwm; //!< high-water-mark over a given interval
#ifdef _OPENMP
#pragma omp threadprivate(private_total_bytes_allocated_in_fabs)
#pragma omp threadprivate(private_total_bytes_allocated_in_fabs_hwm)
#pragma omp threadprivate(private_total_cells_allocated_in_fabs)
#pragma omp threadprivate(private_total_cells_allocated_in_fabs_hwm)
#endif

Long TotalBytesAllocatedInFabs () noexcept;
Long TotalBytesAllocatedInFabsHWM () noexcept;
Long TotalCellsAllocatedInFabs () noexcept;
Long TotalCellsAllocatedInFabsHWM () noexcept;
void ResetTotalBytesAllocatedInFabsHWM () noexcept;
void update_fab_stats (Long n, Long s, std::size_t szt) noexcept;

void BaseFab_Initialize ();
void BaseFab_Finalize ();

struct SrcComp {
    AMREX_GPU_HOST_DEVICE
    explicit SrcComp (int ai) noexcept : i(ai) {}
    int i;
};

struct DestComp {
    AMREX_GPU_HOST_DEVICE
    explicit DestComp (int ai) noexcept : i(ai) {}
    int i;
};

struct NumComps {
    AMREX_GPU_HOST_DEVICE
    explicit NumComps (int an) noexcept : n(an) {}
    int n;
};

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Array4<T>
makeArray4 (T* p, Box const& bx, int ncomp) noexcept
{
    return Array4<T>{p, amrex::begin(bx), amrex::end(bx), ncomp};
}

/**
*  \brief A Fortran Array-like Object
*  BaseFab emulates the Fortran array concept.
*  Useful operations can be performed upon
*  BaseFabs in C++, and they provide a convenient interface to
*  Fortran when it is necessary to retreat into that language.

*  BaseFab is a template class.  Through use of the
*  template, a BaseFab may be based upon any class.  So far at least,
*  most applications have been based upon simple types like integers,
*  real*4s, or real*8s.  Most applications do not use BaseFabs
*  directly, but utilize specialized classes derived from BaseFab.

*  Classes derived from BaseFab include FArrayBox, IArrayBox, TagBox,
*  Mask, EBFArrayBox, EBCellFlag and CutFab.

*  BaseFab objects depend on the dimensionality of space
*  (indirectly through the DOMAIN Box member).  It is
*  typical to define the macro SPACEDIM to be 1, 2, or 3 to indicate
*  the dimension of space.  See the discussion of class Box for more
*  information.  A BaseFab contains a Box DOMAIN, which indicates the
*  integer indexing space over which the array is defined.  A BaseFab
*  also has NVAR components.  By components, we mean that for each
*  point in the rectangular indexing space, there are NVAR values
*  associated with that point.  A Fortran array corresponding to a
*  BaseFab would have (SPACEDIM+1) dimensions.

*  By design, the array layout in a BaseFab mirrors that of a
*  Fortran array.  The first index (x direction for example) varies
*  most rapidly, the next index (y direction), if any, varies next
*  fastest. The component index varies last, after all the spatial
*  indices.

*  It is sometimes convenient to be able to treat a sub-array within an
*  existing BaseFab as a BaseFab in its own right.  This is often
*  referred to as aliasing the BaseFab.  Note that when aliasing is
*  used, the BaseFabs domain will not, in general, be the same as the
*  parent BaseFabs domain, nor will the number of components.
*  BaseFab is a dimension dependent class, so SPACEDIM must be
*  defined as either 1, 2, or 3 when compiling.

*  This is NOT a polymorphic class.

*  It does NOT provide a copy constructor or assignment operator.

*  T MUST have a default constructor and an assignment operator.
*/

template <typename T>
typename std::enable_if<std::is_arithmetic<T>::value>::type
placementNew (T* const /*ptr*/, Long /*n*/)
{}

template <typename T>
typename std::enable_if< AMREX_IS_TRIVIALLY_DEFAULT_CONSTRUCTIBLE(T) && !std::is_arithmetic<T>::value >::type
placementNew (T* const ptr, Long n)
{
    for (Long i = 0; i < n; ++i) {
        new (ptr+i) T;
    }
}

template <typename T>
typename std::enable_if<!AMREX_IS_TRIVIALLY_DEFAULT_CONSTRUCTIBLE(T)>::type
placementNew (T* const ptr, Long n)
{
    AMREX_HOST_DEVICE_FOR_1D ( n, i,
    {
        new (ptr+i) T;
    });
}

template <typename T>
typename std::enable_if<std::is_trivially_destructible<T>::value>::type
placementDelete (T* const /*ptr*/, Long /*n*/)
{}

template <typename T>
typename std::enable_if<!std::is_trivially_destructible<T>::value>::type
placementDelete (T* const ptr, Long n)
{
    AMREX_HOST_DEVICE_FOR_1D (n, i,
    {
        (ptr+i)->~T();
    });
}

struct DataAllocator {
    Arena* m_arena = nullptr;

    DataAllocator ()  noexcept : m_arena(nullptr) {}

    DataAllocator (Arena* ar) noexcept : m_arena(ar) {}

    void* alloc (std::size_t sz) const noexcept {
        return arena()->alloc(sz);
    }

    void free (void* pt) const noexcept {
        arena()->free(pt);
    }

    Arena* arena () const noexcept {
        return (m_arena) ? m_arena : The_Arena();
    }
};

struct DataDeleter {
    Arena* m_arena = nullptr;
    DataDeleter () noexcept : m_arena(nullptr) {}
    DataDeleter (Arena* ar) noexcept : m_arena(ar) {}
    Arena* arena () const noexcept { return (m_arena) ? m_arena : The_Arena(); }
    void operator() (void* pt) const noexcept { arena()->free(pt); }
};

template <class T>
class BaseFab
    : public DataAllocator
{
public:

    template <class U> friend class BaseFab;

    typedef T         value_type;

    //! Construct an empty BaseFab, which must be resized (see BaseFab::resize) before use.
    BaseFab () noexcept;

    explicit BaseFab (Arena* ar) noexcept;

    BaseFab (const Box& bx, int n, Arena* ar);

    //!  Make BaseFab with desired domain (box) and number of components.
    explicit BaseFab (const Box& bx, int n = 1, bool alloc = true,
		      bool shared = false, Arena* ar = nullptr);

    BaseFab (const BaseFab<T>& rhs, MakeType make_type, int scomp, int ncomp);

    /**
     * \brief Create an NON-OWNING BaseFab.  Thus BaseFab is not
     * responsible for memory management.  And it's caller's responsibility that
     * p points to a chunk of memory large enough.
     */
    BaseFab (const Box& bx, int ncomp, T* p);
    BaseFab (const Box& bx, int ncomp, T const* p);

    explicit BaseFab (Array4<T> const& a) noexcept;

    explicit BaseFab (Array4<T> const& a, IndexType t) noexcept;

    explicit BaseFab (Array4<T const> const& a) noexcept;

    explicit BaseFab (Array4<T const> const& a, IndexType t) noexcept;

    //! The destructor deletes the array memory.
    ~BaseFab () noexcept;

    BaseFab (const BaseFab<T>& rhs) = delete;
    BaseFab<T>& operator= (const BaseFab<T>& rhs) = delete;
    BaseFab<T>& operator= (BaseFab<T>&& rhs) = delete;

    BaseFab (BaseFab<T>&& rhs) noexcept;

#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab& operator= (T const&) noexcept;

    static void Initialize();
    static void Finalize();

    /**
    * \brief This function resizes a BaseFab so it covers the Box b
    * with N components.

    * The default action is that under resizing, the memory allocated for the
    * BaseFab only grows and never shrinks.  This function is
    * particularly useful when a BaseFab is used as a temporary
    * space which must be a different size whenever it is used.
    * Resizing is typically faster than re-allocating a
    * BaseFab because memory allocation can often be avoided.
    */
    void resize (const Box& b, int N = 1);

    template <class U=T, typename std::enable_if<std::is_trivially_destructible<U>::value,int>::type = 0>
    Elixir elixir () noexcept;

    /**
     * \brief The function returns the BaseFab to the invalid state.
     * The memory is freed.
     */
    void clear () noexcept;

    // Release ownership of memory
    std::unique_ptr<T,DataDeleter> release () noexcept;

    //! Returns how many bytes used
    std::size_t nBytes () const noexcept { return this->truesize*sizeof(T); }

    std::size_t nBytesOwned () const noexcept {
        return (this->ptr_owner) ? nBytes() : 0;
    }

    //! Returns bytes used in the Box for those components
    std::size_t nBytes (const Box& bx, int ncomps) const noexcept
        { return bx.numPts() * sizeof(T) * ncomps; }

    //! Returns the number of components
    int nComp () const noexcept { return this->nvar; }

    //! for calls to fortran.
    const int* nCompPtr() const noexcept {
        return &(this->nvar);
    }

    //! Returns the number of points
    Long numPts () const noexcept { return this->domain.numPts(); }

    //! Returns the total number of points of all components
    Long size () const noexcept { return this->nvar*this->domain.numPts(); }

    //! Returns the domain (box) where the array is defined
    const Box& box () const noexcept { return this->domain; }

    /**
    * \brief Returns a pointer to an array of SPACEDIM integers
    * giving the length of the domain in each direction
    */
    IntVect length () const noexcept { return this->domain.length(); }

    /**
    * \brief Returns the lower corner of the domain
    * See class Box for analogue.
    */
    const IntVect& smallEnd () const noexcept { return this->domain.smallEnd(); }

    //!  Returns the upper corner of the domain.  See class Box for analogue.
    const IntVect& bigEnd () const noexcept { return this->domain.bigEnd(); }

    /**
    * \brief Returns the lower corner of the domain.

    *Instead of returning them in the form of INTVECTs, as in smallEnd and
    * bigEnd, it returns the values as a pointer to an array of
    * constant integers.  This is useful when interfacing to
    * Fortran subroutines.
    */
    const int* loVect () const noexcept { return this->domain.loVect(); }

    /**
    * \brief Returns the upper corner of the domain.

    *Instead of returning them in the form of INTVECTs, as in smallEnd and
    * bigEnd, it returns the values as a pointer to an array of
    * constant integers.  This is useful when interfacing to
    * Fortran subroutines.
    */
    const int* hiVect () const noexcept { return this->domain.hiVect(); }

    /**
    * \brief Returns true if the domain of fab is totally contained within
    * the domain of this BaseFab.
    */
    bool contains (const BaseFab<T>& fab) const noexcept
    {
        return box().contains(fab.box()) && this->nvar <= fab.nvar;
    }

    /**
    * \brief Returns true if bx is totally contained
    * within the domain of this BaseFab.
    */
    bool contains (const Box& bx) const noexcept { return box().contains(bx); }

    /**
    * \brief Returns a pointer to an object of type T that is the
    * value of the Nth component associated with the cell at the
    * low end of the domain.  This is commonly used to get a pointer
    * to data in the array which is then handed off to a Fortran
    * subroutine.  Remember that data is stored in Fortran array
    * order, with the component index coming last.   In other words,
    * dataPtr returns a pointer to all the Nth components.
    */
    T* dataPtr (int n = 0) noexcept { AMREX_ASSERT(!(this->dptr == 0)); return &(this->dptr[n*this->domain.numPts()]); }

    //! Same as above except works on const FABs.
    const T* dataPtr (int n = 0) const noexcept { AMREX_ASSERT(!(this->dptr == 0)); return &(this->dptr[n*this->domain.numPts()]); }

    T* dataPtr (const IntVect& iv, int n = 0) noexcept;

    const T* dataPtr (const IntVect& iv, int n = 0) const noexcept;

    void setPtr (T* p, Long sz) noexcept { AMREX_ASSERT(this->dptr == 0 && this->truesize == 0); this->dptr = p; this->truesize = sz; }

    void prefetchToHost () const noexcept;
    void prefetchToDevice () const noexcept;

    AMREX_FORCE_INLINE
    Array4<T const> array () const noexcept
    {
        return makeArray4<T const>(this->dptr, this->domain, this->nvar);
    }

    AMREX_FORCE_INLINE
    Array4<T const> array (int start_comp) const noexcept
    {
        return Array4<T const>(makeArray4<T const>(this->dptr, this->domain, this->nvar),start_comp);
    }

    AMREX_FORCE_INLINE
    Array4<T> array () noexcept
    {
        return makeArray4<T>(this->dptr, this->domain, this->nvar);
    }

    AMREX_FORCE_INLINE
    Array4<T> array (int start_comp) noexcept
    {
        return Array4<T>(makeArray4<T>(this->dptr, this->domain, this->nvar),start_comp);
    }

    AMREX_FORCE_INLINE
    Array4<T const> const_array () const noexcept
    {
        return makeArray4<T const>(this->dptr, this->domain, this->nvar);
    }

    AMREX_FORCE_INLINE
    Array4<T const> const_array (int start_comp) const noexcept
    {
        return Array4<T const>(makeArray4<T const>(this->dptr, this->domain, this->nvar),start_comp);
    }

    //! Returns true if the data for the FAB has been allocated.
    bool isAllocated () const noexcept { return this->dptr != 0; }

    /**
    * \brief Returns a reference to the Nth component value
    * defined at position p in the domain.  This operator may be
    * inefficient if the C++ compiler is unable to optimize the
    * C++ code.
    */
    T& operator() (const IntVect& p, int N) noexcept;

    //! Same as above, except returns component 0.
    T& operator() (const IntVect& p) noexcept;

    //! Same as above except works on const FABs.
    const T& operator() (const IntVect& p, int N) const noexcept;

    //! Same as above, except returns component 0.
    const T& operator() (const IntVect& p) const noexcept;

    /**
    * \brief This function puts numcomp component values, starting at
    * component N, from position pos in the domain into array data,
    * that must be allocated by the user.
    */
    void getVal (T* data, const IntVect& pos, int N, int numcomp) const noexcept;
    //! Same as above, except that starts at component 0 and copies all comps.
    void getVal (T* data, const IntVect& pos) const noexcept;
    /**
    * \brief The setVal functions set sub-regions in the BaseFab to a
    * constant value.  This most general form specifies the sub-box,
    * the starting component number, and the number of components
    * to be set.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setVal (T const& x, const Box& bx, int nstart, int ncomp) noexcept;
    //! Same as above, except the number of modified components is one. N is the component to be modified.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setVal (T const& x, const Box& bx, int N = 0) noexcept;
    //! Same as above, except the sub-box defaults to the entire domain.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setVal (T const& x, int N) noexcept;

#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setValIfNot (T const& x, const Box& bx, const BaseFab<int>& mask, int nstart, int ncomp) noexcept;

    /**
    * \brief This function is analogous to the fourth form of
    * setVal above, except that instead of setting values on the
    * Box b, values are set on the complement of b in the domain.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setComplement (T const& x, const Box& b, int ns, int num) noexcept;

    /**
    * \brief The copy functions copy the contents of one BaseFab into
    * another.  The destination BaseFab is always the object which
    * invokes the function.  This, the most general form of copy,
    * specifies the contents of any sub-box srcbox in BaseFab src
    * may be copied into a (possibly different) destbox in the
    * destination BaseFab.  Note that although the srcbox and the
    * destbox may be disjoint, they must be the same size and shape.
    * If the sizes differ, the copy is undefined and a runtime error
    * results.  This copy function is the only one of the copy
    * functions to allow a copy between differing boxes. The user
    * also specifies how many components are copied, starting at
    * component srccomp in src and stored starting at component
    * destcomp. The results are UNDEFINED if the src and dest are the
    * same and the srcbox and destbox overlap.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& copy (const BaseFab<T>& src, const Box& srcbox, int srccomp,
                      const Box& destbox, int destcomp, int numcomp) noexcept;

    /**
    * \brief As above, except the destination Box and the source Box
    * are taken to be the entire domain of the destination.   A copy
    * of the intersecting region is performed.
    * class.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& copy (const BaseFab<T>& src, int srccomp, int destcomp,
                      int numcomp = 1) noexcept;
    /**
    * \brief As above, except that the destination Box is specified,
    * but the source Box is taken to the equal to the source
    * Box, and all components of the destination BaseFab are
    * copied.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& copy (const BaseFab<T>& src, const Box& destbox) noexcept;

    //! Copy from the srcbox of this Fab to raw memory and return the number of bytes copied
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    std::size_t copyToMem (const Box& srcbox, int srccomp,
                           int numcomp, void* dst) const noexcept;

    //! Copy from raw memory to the dstbox of this Fab and return the number of bytes copied
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    std::size_t copyFromMem (const Box& dstbox, int dstcomp,
                             int numcomp, const void* src) noexcept;

    //! Add from raw memory to the dstbox of this Fab and return the number of bytes copied
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    std::size_t addFromMem (const Box& dstbox, int dstcomp,
                            int numcomp, const void* src) noexcept;

    /**
    * \brief Perform shifts upon the domain of the BaseFab. They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    BaseFab<T>& shift (const IntVect& v) noexcept;
    /**
    * \brief Perform shifts upon the domain of the BaseFab.  They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    BaseFab<T>& shift (int idir, int n_cell) noexcept;
    /**
    * \brief Perform shifts upon the domain of the BaseFab.  They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    BaseFab<T>& shiftHalf (int dir, int num_halfs) noexcept;
    /**
    * \brief Perform shifts upon the domain of the BaseFab. They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    BaseFab<T>& shiftHalf (const IntVect& num_halfs) noexcept;

#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Real norminfmask (const Box& subbox, const BaseFab<int>& mask, int scomp=0, int ncomp=1) const noexcept;

    /**
    * \brief Compute the Lp-norm of this FAB using components (scomp : scomp+ncomp-1).
    *   p < 0  -> ERROR
    *   p = 0  -> infinity norm (max norm)
    *   p = 1  -> sum of ABS(FAB)
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Real norm (int p, int scomp = 0, int ncomp = 1) const noexcept;

    //! Same as above except only on given subbox.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    Real norm (const Box& subbox, int p, int scomp = 0, int ncomp = 1) const noexcept;
    //!Compute absolute value for all components of this FAB.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void abs () noexcept;
    //! Same as above except only for components (comp: comp+numcomp-1)
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void abs (int comp, int numcomp=1) noexcept;
    /**
    * \brief Calculate abs() on subbox for given component range.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void abs (const Box& subbox, int comp = 0, int numcomp=1) noexcept;
    /**
    * \return Minimum value of given component.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T min (int comp = 0) const noexcept;
    /**
    * \return Minimum value of given component in given subbox.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T min (const Box& subbox, int comp = 0) const noexcept;
    /**
    * \return Maximum value of given component.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T max (int comp = 0) const noexcept;
    /**
    * \return Maximum value of given component in given subbox.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T max (const Box& subbox, int comp = 0) const noexcept;
    /**
    * \return Maximum of the absolute value of given component.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T maxabs (int comp = 0) const noexcept;
    /**
    * \return Maximum of the absolute value of given component in given subbox.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T maxabs (const Box& subbox, int comp = 0) const noexcept;

    /*(
    * \return location of given value
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    IntVect indexFromValue (const Box& subbox, int comp, T const& value) const noexcept;

    /**
    * \return location of minimum value in given component.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    IntVect minIndex (int comp = 0) const noexcept;
    /**
    * \return location of minimum value in given component in
    * given subbox.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    IntVect minIndex (const Box& subbox, int comp = 0) const noexcept;
    /**
    * \return return mininum value and location to allow
    * efficient looping over multiple boxes.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void  minIndex (const Box& subbox, Real& min_val, IntVect& min_idx, int comp = 0) const noexcept;

    /**
    * \return location of maximum value in given component.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    IntVect maxIndex (int comp = 0) const noexcept;
    /**
    * \return location of maximum value in given component in given
    * subbox.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    IntVect maxIndex (const Box& subbox, int comp = 0) const noexcept;
    /**
    * \return return maximum value and location to allow
    * efficient looping over multiple boxes.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void  maxIndex (const Box& subbox, Real& max_value, IntVect& max_index, int comp = 0) const noexcept;

    /**
    * \brief Compute mask array with value of 1 in cells where
    * BaseFab has value less than val, 0 otherwise.
    * mask is resized by this function.
    * The number of cells marked with 1 returned.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    int maskLT (BaseFab<int>& mask, T const& val, int comp = 0) const noexcept;
    //! Same as above except mark cells with value less than or equal to val.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    int maskLE (BaseFab<int>& mask, T const& val, int comp = 0) const noexcept;

    //! Same as above except mark cells with value equal to val.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    int maskEQ (BaseFab<int>& mask, T const& val, int comp = 0) const noexcept;
    //! Same as above except mark cells with value greater than val.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    int maskGT (BaseFab<int>& mask, T const& val, int comp = 0) const noexcept;
    //! Same as above except mark cells with value greater than or equal to val.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    int maskGE (BaseFab<int>& mask, T const& val, int comp = 0) const noexcept;

    //! Returns sum of given component of FAB state vector.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T sum (int comp, int numcomp = 1) const noexcept;
    //! Compute sum of given component of FAB state vector in given subbox.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T sum (const Box& subbox, int comp, int numcomp = 1) const noexcept;

    //! Most general version, specify subbox and which components.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& invert (T const& v, const Box& subbox, int comp=0, int numcomp=1) noexcept;
    //! As above except on entire domain.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& invert (T const& v, int comp, int numcomp=1) noexcept;

    //! Negate BaseFab, most general.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& negate (const Box& subbox, int comp=0, int numcomp=1) noexcept;
    //! As above, except on entire domain.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& negate (int comp, int numcomp=1) noexcept;

    //! Scalar addition (a[i] <- a[i] + r), most general.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& plus (T const& r, const Box& b, int comp=0, int numcomp=1) noexcept;

    //! As above, except on entire domain.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& plus (T const& r, int comp, int numcomp=1) noexcept;

    /**
    * \brief Add src components (srccomp:srccomp+numcomp-1) to
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& plus (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp=1) noexcept;
    /**
    * \brief Same as above except addition is restricted to intersection
    * of subbox and src FAB. NOTE: subbox must be contained in this
    * FAB.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& plus (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp, int numcomp=1) noexcept;
    /**
    * \brief Add srcbox region of src FAB to destbox region of this FAB.
    * The srcbox and destbox must be same size.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& plus (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                      int srccomp, int destcomp, int numcomp=1) noexcept;

    //! Atomic FAB addition (a[i] <- a[i] + b[i]).
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& atomicAdd (const BaseFab<T>& src) noexcept;

    /**
    * \brief Atomically add src components (srccomp:srccomp+numcomp-1) to
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& atomicAdd (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp=1) noexcept;
    /**
    * \brief Same as above except addition is restricted to intersection
    * of subbox and src FAB. NOTE: subbox must be contained in this
    * FAB.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& atomicAdd (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp,
                           int numcomp=1) noexcept;
    /**
    * \brief Atomically add srcbox region of src FAB to destbox region of this FAB.
    * The srcbox and destbox must be same size.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& atomicAdd (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                           int srccomp, int destcomp, int numcomp=1) noexcept;

    //! FAB SAXPY (y[i] <- y[i] + a * x[i]), in place.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& saxpy (T a, const BaseFab<T>& x, const Box& srcbox, const Box& destbox,
                       int srccomp, int destcomp, int numcomp=1) noexcept;
    //! FAB SAXPY (y[i] <- y[i] + a * x[i]), in place.  All components.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& saxpy (T a, const BaseFab<T>& x) noexcept;

    //! FAB XPAY (y[i] <- x[i] + a * y[i])
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& xpay (T a, const BaseFab<T>& x, const Box& srcbox, const Box& destbox,
                      int srccomp, int destcomp, int numcomp=1) noexcept;

    //! y[i] <- y[i] + x1[i] * x2[i])
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& addproduct (const Box& destbox, int destcomp, int numcomp,
                            const BaseFab<T>& src1, int comp1,
                            const BaseFab<T>& src2, int comp2) noexcept;

    /**
    * \brief Subtract src components (srccomp:srccomp+numcomp-1) to
    * this FABs components (destcomp:destcomp+numcomp-1) where
    * the two FABs intersect.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& minus (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp=1) noexcept;
    /**
    * \brief Same as above except subtraction is restricted to intersection
    * of subbox and src FAB.  NOTE: subbox must be contained in
    * this FAB.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& minus (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp,
                       int numcomp=1) noexcept;
    /**
    * \brief Subtract srcbox region of src FAB from destbox region
    * of this FAB. srcbox and destbox must be same size.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& minus (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                       int srccomp, int destcomp, int numcomp=1) noexcept;

    //! Scalar multiplication, except control which components are multiplied.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& mult (T const& r, int comp, int numcomp=1) noexcept;
    /**
    * \brief As above, except specify sub-box.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& mult (T const& r, const Box& b, int comp=0, int numcomp=1) noexcept;

    /**
    * \brief Multiply src components (srccomp:srccomp+numcomp-1) with
    * this FABs components (destcomp:destcomp+numcomp-1) where
    * the two FABs intersect.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& mult (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp=1) noexcept;

    /**
    * \brief Same as above except multiplication is restricted to
    * intersection of subbox and src FAB.  NOTE: subbox must be
    * contained in this FAB.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& mult (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp,
                      int numcomp=1) noexcept;

    /**
    * \brief Multiply srcbox region of src FAB with destbox region
    * of this FAB. The srcbox and destbox must be same size.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& mult (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                      int srccomp, int destcomp, int numcomp=1) noexcept;

    //! As above except specify which components.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& divide (T const& r, int comp, int numcomp=1) noexcept;

    //! As above except specify sub-box.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& divide (T const& r, const Box& b, int comp=0, int numcomp=1) noexcept;

    /**
    * \brief This FAB is numerator, src FAB is denominator
    * divide src components (srccomp:srccomp+numcomp-1) into
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& divide (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp=1) noexcept;
    /**
    * \brief Same as above except division is restricted to
    * intersection of subbox and src FAB.  NOTE: subbox must be
    * contained in this FAB.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& divide (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp,
                        int numcomp=1) noexcept;
    /**
    * \brief destbox region of this FAB is numerator. srcbox regions of
    * src FAB is denominator. srcbox and destbox must be same size.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& divide (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                        int srccomp, int destcomp, int numcomp=1) noexcept;
    /**
    * \brief Divide wherever "src" is "true" or "non-zero".
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& protected_divide (const BaseFab<T>& src) noexcept;

    /**
    * \brief Divide wherever "src" is "true" or "non-zero".
    * This FAB is numerator, src FAB is denominator
    * divide src components (srccomp:srccomp+numcomp-1) into
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& protected_divide (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp=1) noexcept;

    /**
    * \brief Divide wherever "src" is "true" or "non-zero".
    * Same as above except division is restricted to
    * intersection of subbox and src FAB.  NOTE: subbox must be
    * contained in this FAB.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& protected_divide (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp,
                                  int numcomp=1) noexcept;

    /**
    * Divide wherever "src" is "true" or "non-zero".
    * destbox region of this FAB is numerator. srcbox regions of
    * src FAB is denominator. srcbox and destbox must be same size.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& protected_divide (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                                  int srccomp, int destcomp, int numcomp=1) noexcept;

    /**
    * \brief Linear interpolation / extrapolation.
    * Result is (t2-t)/(t2-t1)*f1 + (t-t1)/(t2-t1)*f2
    * Data is taken from b1 region of f1, b2 region of f2
    * and stored in b region of this FAB.
    * Boxes b, b1 and b2 must be the same size.
    * Data is taken from component comp1 of f1, comp2 of f2,
    * and stored in component comp of this FAB.
    * This FAB is returned as a reference for chaining.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& linInterp (const BaseFab<T>& f1, const Box& b1, int comp1,
                           const BaseFab<T>& f2, const Box& b2, int comp2,
                           Real t1, Real t2, Real t,
                           const Box& b, int comp, int numcomp = 1) noexcept;

    //! Version of linInterp() in which b, b1, & b2 are the same.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& linInterp (const BaseFab<T>& f1, int comp1,
                           const BaseFab<T>& f2, int comp2,
                           Real t1, Real t2, Real t,
                           const Box& b, int comp, int numcomp = 1) noexcept;

    /**
    * \brief Linear combination.  Result is alpha*f1 + beta*f2.
    * Data is taken from b1 region of f1, b2 region of f2
    * and stored in b region of this FAB.
    * Boxes b, b1 and b2 must be the same size.
    * Data is taken from component comp1 of f1, comp2 of f2,
    * and stored in component comp of this FAB.
    * This FAB is returned as a reference for chaining.
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& linComb (const BaseFab<T>& f1, const Box& b1, int comp1,
                         const BaseFab<T>& f2, const Box& b2, int comp2,
                         Real alpha, Real beta, const Box& b,
                         int comp, int numcomp = 1) noexcept;

    //! Dot product of x (i.e.,this) and y
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T dot (const Box& xbx, int xcomp, const BaseFab<T>& y, const Box& ybx, int ycomp,
	   int numcomp = 1) const noexcept;

#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T dotmask (const BaseFab<int>& mask, const Box& xbx, int xcomp,
               const BaseFab<T>& y, const Box& ybx, int ycomp,
               int numcomp) const noexcept;

    //! Change the Box type without change the length
    void SetBoxType (const IndexType& typ) noexcept { this->domain.setType(typ); }

    //
    // New interfaces
    //

    //! Set value on the whole domain and all components
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setVal (T const& val) noexcept;
    //
    //! Do nothing if bx is empty.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setVal (T const& val, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept;

#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setValIf (T const& val, const BaseFab<int>& mask) noexcept;
    //
    //! Do nothing if bx is empty.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setValIf (T const& val, Box const& bx, const BaseFab<int>& mask, DestComp dcomp, NumComps ncomp) noexcept;

#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setValIfNot (T const& val, const BaseFab<int>& mask) noexcept;
    //
    //! Do nothing if bx is empty.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setValIfNot (T const& val, Box const& bx, const BaseFab<int>& mask, DestComp dcomp, NumComps ncomp) noexcept;

    //! setVal on the complement of bx in the fab's domain
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    void setComplement (T const& val, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept;

    /**
    * copy is performed on the intersection of dest and src fabs.
    * All components of dest fab are copied. src fab must have enough
    * components (more is OK).
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& copy (const BaseFab<T>& src) noexcept;
    //
    //! Do nothing if bx does not intersect with src fab.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& copy (const BaseFab<T>& src, Box bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept;

    //! Scalar addition on the whole domain and all components
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& plus (T const& val) noexcept;
    //
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& operator+= (T const& val) noexcept;
    //
    //! Do nothing if bx is empty.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& plus (T const& val, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept;
    /**
    * Fab addition is performed on the intersection of dest and src fabs.
    * All components of dest fab are copied. src fab must have enough
    * components (more is OK).
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& plus (const BaseFab<T>& src) noexcept;
    //
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& operator+= (const BaseFab<T>& src) noexcept;
    //
    //! Do nothing if bx does not intersect with src fab.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& plus (const BaseFab<T>& src, Box bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept;

    //! Scalar subtraction on the whole domain and all components
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& minus (T const& val) noexcept;
    //
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& operator-= (T const& val) noexcept;
    //
    //! Do nothing if bx is empty.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& minus (T const& val, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept;
    /**
    * Fab subtraction is performed on the intersection of dest and src fabs.
    * All components of dest fab are copied. src fab must have enough
    * components (more is OK).
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& minus (const BaseFab<T>& src) noexcept;
    //
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& operator-= (const BaseFab<T>& src) noexcept;
    //
    //! Do nothing if bx does not intersect with src fab.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& minus (const BaseFab<T>& src, Box bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept;

    //! Scalar multiplication on the whole domain and all components
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& mult (T const& val) noexcept;
    //
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& operator*= (T const& val) noexcept;
    //
    //! Do nothing if bx is empty.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& mult (T const& val, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept;
    /**
    * Fab multiplication is performed on the intersection of dest and src fabs.
    * All components of dest fab are copied. src fab must have enough
    * components (more is OK).
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& mult (const BaseFab<T>& src) noexcept;
    //
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& operator*= (const BaseFab<T>& src) noexcept;
    //
    //! Do nothing if bx does not intersect with src fab.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& mult (const BaseFab<T>& src, Box bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept;

    //! Scalar division on the whole domain and all components
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& divide (T const& val) noexcept;
    //
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& operator/= (T const& val) noexcept;
    //
    //! Do nothing if bx is empty.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& divide (T const& val, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept;
    /**
    * Fab division is performed on the intersection of dest and src fabs.
    * All components of dest fab are copied. src fab must have enough
    * components (more is OK).
    */
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& divide (const BaseFab<T>& src) noexcept;
    //
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& operator/= (const BaseFab<T>& src) noexcept;
    //
    //! Do nothing if bx does not intersect with src fab.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& divide (const BaseFab<T>& src, Box bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept;

    //! on the whole domain and all components
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& negate () noexcept;
    //
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& negate (const Box& bx, DestComp dcomp, NumComps ncomp) noexcept;

    //! Fab <- Fab/r on the whole domain and all components
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& invert (T const& r) noexcept;
    //
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    BaseFab<T>& invert (T const& r, const Box& bx, DestComp dcomp, NumComps ncomp) noexcept;

    //! Sum
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T sum (const Box& bx, DestComp dcomp, NumComps ncomp) const noexcept;

    //! Dot product of two Fabs
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T dot (const BaseFab<T>& src, const Box& bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) const noexcept;

    //! Int wrapper for dot.
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T dot (const Box& bx, int destcomp, int numcomp) const noexcept;

    //! Dot product
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T dot (const Box& bx, DestComp dcomp, NumComps ncomp) const noexcept;

    //! Dot product of two Fabs with mask
#if defined(AMREX_USE_GPU)
    template <RunOn run_on>
#else
    template <RunOn run_on=RunOn::Host>
#endif
    T dotmask (const BaseFab<T>& src, const Box& bx, const BaseFab<int>& mask,
               SrcComp scomp, DestComp dcomp, NumComps ncomp) const noexcept;

protected:
    //! Allocates memory for the BaseFab<T>.
    void define ();

    T*   dptr = nullptr;        //!< The data pointer.
    Box  domain;                //!< My index space.
    int  nvar = 0;              //!< Number components.
    Long truesize = 0L;         //!< nvar*numpts that was allocated on heap.
    bool ptr_owner = false;     //!< Owner of T*?
    bool shared_memory = false; //!< Is the memory allocated in shared memory?
};

template <class T>
AMREX_FORCE_INLINE
T*
BaseFab<T>::dataPtr (const IntVect& p, int n) noexcept
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr + (this->domain.index(p)+n*this->domain.numPts());
}

template <class T>
AMREX_FORCE_INLINE
const T*
BaseFab<T>::dataPtr (const IntVect& p, int n) const noexcept
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr + (this->domain.index(p)+n*this->domain.numPts());
}

template <class T>
void
BaseFab<T>::prefetchToHost () const noexcept
{
#ifdef AMREX_USE_GPU
    if (this->arena() == The_Arena() or this->arena() == The_Managed_Arena()) {
#if defined(AMREX_USE_DPCPP)
        // xxxxx DPCPP todo: prefetchToHost
        // std::size_t s = sizeof(T)*this->nvar*this->domain.numPts();
        // auto& q = Gpu::Device::streamQueue();
        // q.submit([&] (sycl::handler& h) { h.prefetch(this->dptr, s); });
#elif defined(AMREX_USE_CUDA)
        std::size_t s = sizeof(T)*this->nvar*this->domain.numPts();
        AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(this->dptr, s,
                                                  cudaCpuDeviceId,
                                                  Gpu::gpuStream()));
#elif defined(AMREX_USE_HIP)
        // xxxxx HIP FIX HERE after managed memory is supported
#endif
    }
#endif
}

template <class T>
void
BaseFab<T>::prefetchToDevice () const noexcept
{
#ifdef AMREX_USE_GPU
    if (this->arena() == The_Arena() or this->arena() == The_Managed_Arena()) {
#if defined(AMREX_USE_DPCPP)
        std::size_t s = sizeof(T)*this->nvar*this->domain.numPts();
        auto& q = Gpu::Device::streamQueue();
        q.submit([&] (sycl::handler& h) { h.prefetch(this->dptr, s); });
#elif defined(AMREX_USE_CUDA)
        std::size_t s = sizeof(T)*this->nvar*this->domain.numPts();
        AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(this->dptr, s,
                                                  Gpu::Device::deviceId(),
                                                  Gpu::gpuStream()));
#elif defined(AMREX_USE_HIP)
        // xxxxx HIP FIX HERE after managed memory is supported
#endif
    }
#endif
}

template <class T>
AMREX_FORCE_INLINE
T&
BaseFab<T>::operator() (const IntVect& p, int n) noexcept
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr[this->domain.index(p)+n*this->domain.numPts()];
}

template <class T>
AMREX_FORCE_INLINE
T&
BaseFab<T>::operator() (const IntVect& p) noexcept
{
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr[this->domain.index(p)];
}

template <class T>
AMREX_FORCE_INLINE
const T&
BaseFab<T>::operator() (const IntVect& p, int n) const noexcept
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr[this->domain.index(p)+n*this->domain.numPts()];
}

template <class T>
AMREX_FORCE_INLINE
const T&
BaseFab<T>::operator() (const IntVect& p) const noexcept
{
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr[this->domain.index(p)];
}

template <class T>
void
BaseFab<T>::getVal  (T*             data,
                     const IntVect& pos,
                     int            n,
                     int            numcomp) const noexcept
{
    const int loc      = this->domain.index(pos);
    const Long sz    = this->domain.numPts();

    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(n >= 0 && n + numcomp <= this->nvar);

    for (int k = 0; k < numcomp; k++)
        data[k] = this->dptr[loc+(n+k)*sz];
}

template <class T>
void
BaseFab<T>::getVal (T*             data,
                    const IntVect& pos) const noexcept
{
    getVal(data,pos,0,this->nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::shift (const IntVect& v) noexcept
{
    this->domain += v;
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::shift (int idir, int n_cell) noexcept
{
    this->domain.shift(idir,n_cell);
    return *this;
}

template <class T>
BaseFab<T> &
BaseFab<T>::shiftHalf (const IntVect& v) noexcept
{
    this->domain.shiftHalf(v);
    return *this;
}

template <class T>
BaseFab<T> &
BaseFab<T>::shiftHalf (int idir, int n_cell) noexcept
{
    this->domain.shiftHalf(idir,n_cell);
    return *this;
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setVal (T const& x, const Box& bx, int n) noexcept
{
    this->setVal<run_on>(x, bx, DestComp{n}, NumComps{1});
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setVal (T const& x, int n) noexcept
{
    this->setVal<run_on>(x, this->domain, DestComp{n}, NumComps{1});
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setVal (T const& x, const Box& bx, int dcomp, int ncomp) noexcept
{
    this->setVal<run_on>(x, bx, DestComp{dcomp}, NumComps{ncomp});
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setValIfNot (T const& val, const Box& bx, const BaseFab<int>& mask, int ns, int num) noexcept
{
    this->setValIfNot<run_on>(val, bx, mask, DestComp{ns}, NumComps{num});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src, const Box& srcbox, int srccomp,
                  const Box& destbox, int destcomp, int numcomp) noexcept
{
    AMREX_ASSERT(destbox.ok());
    AMREX_ASSERT(srcbox.sameSize(destbox));
    AMREX_ASSERT(src.box().contains(srcbox));
    AMREX_ASSERT(this->domain.contains(destbox));
    AMREX_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    AMREX_ASSERT(destcomp >= 0 && destcomp+numcomp <= this->nvar);

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};

    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        d(i,j,k,n+destcomp) = s(i+offset.x,j+offset.y,k+offset.z,n+srccomp);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src, const Box& destbox) noexcept
{
    return this->copy<run_on>(src, destbox, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    return copy<run_on>(src, this->domain, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T>
void
BaseFab<T>::define ()
{
    AMREX_ASSERT(this->dptr == 0);
    AMREX_ASSERT(this->domain.numPts() > 0);
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
    AMREX_ASSERT(std::numeric_limits<Long>::max()/this->nvar > this->domain.numPts());
#else
    AMREX_ASSERT(LONG_MAX/this->nvar > this->domain.numPts());
#endif

    AMREX_ASSERT(this->nvar >= 0);
    if (this->nvar == 0) return;

    this->truesize  = this->nvar*this->domain.numPts();
    this->ptr_owner = true;
    this->dptr = static_cast<T*>(this->alloc(this->truesize*sizeof(T)));

    placementNew(this->dptr, this->truesize);

    amrex::update_fab_stats(this->domain.numPts(), this->truesize, sizeof(T));
}

template <class T>
BaseFab<T>::BaseFab () noexcept
{}

template <class T>
BaseFab<T>::BaseFab (Arena* ar) noexcept
    : DataAllocator{ar}
{}

template <class T>
BaseFab<T>::BaseFab (const Box& bx, int n, Arena* ar)
    : DataAllocator{ar}, domain(bx), nvar(n)
{
    define();
}

template <class T>
BaseFab<T>::BaseFab (const Box& bx, int n, bool alloc, bool shared, Arena* ar)
    : DataAllocator{ar}, domain(bx), nvar(n), shared_memory(shared)
{
    if (!this->shared_memory && alloc) define();
}

template <class T>
BaseFab<T>::BaseFab (const BaseFab<T>& rhs, MakeType make_type, int scomp, int ncomp)
    : DataAllocator{rhs.arena()}, 
      dptr(const_cast<T*>(rhs.dataPtr(scomp))),
      domain(rhs.domain), nvar(ncomp),
      truesize(ncomp*rhs.domain.numPts())
{
    AMREX_ASSERT(scomp+ncomp <= rhs.nComp());
    if (make_type == amrex::make_deep_copy)
    {
        this->dptr = nullptr;
        define();
        this->copy<RunOn::Device>(rhs, this->domain, scomp, this->domain, 0, ncomp);
    } else if (make_type == amrex::make_alias) {
        ; // nothing to do
    } else {
        amrex::Abort("BaseFab: unknown MakeType");
    }
}

template<class T>
BaseFab<T>::BaseFab (const Box& bx, int ncomp, T* p)
    : dptr(p), domain(bx), nvar(ncomp), truesize(bx.numPts()*ncomp)
{
}

template<class T>
BaseFab<T>::BaseFab (const Box& bx, int ncomp, T const* p)
    : dptr(const_cast<T*>(p)), domain(bx), nvar(ncomp), truesize(bx.numPts()*ncomp)
{
}

template<class T>
BaseFab<T>::BaseFab (Array4<T> const& a) noexcept
    : dptr(a.p),
      domain(IntVect(AMREX_D_DECL(a.begin.x,a.begin.y,a.begin.z)),
             IntVect(AMREX_D_DECL(a.end.x-1,a.end.y-1,a.end.z-1))),
      nvar(a.ncomp), truesize(a.ncomp*a.nstride)
{}

template<class T>
BaseFab<T>::BaseFab (Array4<T> const& a, IndexType t) noexcept
    : dptr(a.p),
      domain(IntVect(AMREX_D_DECL(a.begin.x,a.begin.y,a.begin.z)),
             IntVect(AMREX_D_DECL(a.end.x-1,a.end.y-1,a.end.z-1)), t),
      nvar(a.ncomp), truesize(a.ncomp*a.nstride)
{}

template<class T>
BaseFab<T>::BaseFab (Array4<T const> const& a) noexcept
    : dptr(const_cast<T*>(a.p)),
      domain(IntVect(AMREX_D_DECL(a.begin.x,a.begin.y,a.begin.z)),
             IntVect(AMREX_D_DECL(a.end.x-1,a.end.y-1,a.end.z-1))),
      nvar(a.ncomp), truesize(a.ncomp*a.nstride)
{}

template<class T>
BaseFab<T>::BaseFab (Array4<T const> const& a, IndexType t) noexcept
    : dptr(const_cast<T*>(a.p)),
      domain(IntVect(AMREX_D_DECL(a.begin.x,a.begin.y,a.begin.z)),
             IntVect(AMREX_D_DECL(a.end.x-1,a.end.y-1,a.end.z-1)), t),
      nvar(a.ncomp), truesize(a.ncomp*a.nstride)
{}

template <class T>
BaseFab<T>::~BaseFab () noexcept
{
    clear();
}

template <class T>
BaseFab<T>::BaseFab (BaseFab<T>&& rhs) noexcept
    : DataAllocator{rhs.arena()},
      dptr(rhs.dptr), domain(rhs.domain),
      nvar(rhs.nvar), truesize(rhs.truesize),
      ptr_owner(rhs.ptr_owner), shared_memory(rhs.shared_memory)
{
    rhs.dptr = nullptr;
    rhs.ptr_owner = false;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::operator= (T const& t) noexcept
{
    setVal<run_on>(t);
    return *this;
}

template <class T>
void
BaseFab<T>::resize (const Box& b, int n)
{
    this->nvar   = n;
    this->domain = b;

    if (this->dptr == 0 || !this->ptr_owner)
    {
	if (this->shared_memory)
	    amrex::Abort("BaseFab::resize: BaseFab in shared memory cannot increase size");

        this->dptr = nullptr;
        define();
    }
    else if (this->nvar*this->domain.numPts() > this->truesize)
    {
	if (this->shared_memory)
	    amrex::Abort("BaseFab::resize: BaseFab in shared memory cannot increase size");

        clear();

        define();
    }
}

template <class T>
template <class U, typename std::enable_if<std::is_trivially_destructible<U>::value,int>::type>
Elixir
BaseFab<T>::elixir () noexcept
{
    bool o;
    if (Gpu::inLaunchRegion()) {
        o = this->ptr_owner;
        this->ptr_owner = false;
        if (o && this->dptr) {
            if (this->nvar > 1) {
                amrex::update_fab_stats(-this->truesize/this->nvar, -this->truesize, sizeof(T));
            } else {
                amrex::update_fab_stats(0, -this->truesize, sizeof(T));
            }
        }
    } else {
        o = false;
    }
    return Elixir((o ? this->dptr : nullptr), this->arena());
}

template <class T>
void
BaseFab<T>::clear () noexcept
{
    if (this->dptr)
    {
        //
        // Call T::~T() on the to-be-destroyed memory.
        //
	if (this->ptr_owner)
	{
	    if (this->shared_memory)
	    {
		amrex::Abort("BaseFab::clear: BaseFab cannot be owner of shared memory");
	    }

            placementDelete(this->dptr, this->truesize);

	    this->free(this->dptr);

	    if (this->nvar > 1) {
		amrex::update_fab_stats(-this->truesize/this->nvar, -this->truesize, sizeof(T));
	    } else {
		amrex::update_fab_stats(0, -this->truesize, sizeof(T));
	    }
	}

	this->dptr = 0;
	this->truesize = 0;
    }
}

template <class T>
std::unique_ptr<T,DataDeleter>
BaseFab<T>::release () noexcept
{
    std::unique_ptr<T,DataDeleter> r(nullptr, DataDeleter{this->arena()});
    if (this->dptr && this->ptr_owner) {
        r.reset(this->dptr);
        this->ptr_owner = false;
        if (this->nvar > 1) {
            amrex::update_fab_stats(-this->truesize/this->nvar, -this->truesize, sizeof(T));
        } else {
            amrex::update_fab_stats(0, -this->truesize, sizeof(T));
        }
    }
    return r;
}

template <class T>
template <RunOn run_on>
std::size_t
BaseFab<T>::copyToMem (const Box& srcbox,
                       int        srccomp,
                       int        numcomp,
                       void*      dst) const noexcept
{
    BL_ASSERT(box().contains(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= nComp());

    if (srcbox.ok())
    {
        Array4<T> d(static_cast<T*>(dst),amrex::begin(srcbox),amrex::end(srcbox),numcomp);
        Array4<T const> const& s = this->const_array();
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, srcbox, numcomp, i, j, k, n,
        {
            d(i,j,k,n) = s(i,j,k,n+srccomp);
        });
        return sizeof(T)*d.size();
    }
    else
    {
        return 0;
    }
}

template <class T>
template <RunOn run_on>
std::size_t
BaseFab<T>::copyFromMem (const Box&  dstbox,
                         int         dstcomp,
                         int         numcomp,
                         const void* src) noexcept
{
    BL_ASSERT(box().contains(dstbox));
    BL_ASSERT(dstcomp >= 0 && dstcomp+numcomp <= nComp());

    if (dstbox.ok())
    {
        Array4<T const> s(static_cast<T const*>(src), amrex::begin(dstbox),
                          amrex::end(dstbox), numcomp);
        Array4<T> const& d = this->array();
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, dstbox, numcomp, i, j, k, n,
        {
            d(i,j,k,n+dstcomp) = s(i,j,k,n);
        });
        return sizeof(T)*s.size();
    }
    else
    {
        return 0;
    }
}

template <class T>
template <RunOn run_on>
std::size_t
BaseFab<T>::addFromMem (const Box&  dstbox,
                        int         dstcomp,
                        int         numcomp,
                        const void* src) noexcept
{
    BL_ASSERT(box().contains(dstbox));
    BL_ASSERT(dstcomp >= 0 && dstcomp+numcomp <= nComp());

    if (dstbox.ok())
    {
        Array4<T const> s(static_cast<T const*>(src), amrex::begin(dstbox),
                          amrex::end(dstbox), numcomp);
        Array4<T> const& d = this->array();
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, dstbox, numcomp, i, j, k, n,
        {
            d(i,j,k,n+dstcomp) += s(i,j,k,n);
        });
        return sizeof(T)*s.size();
    }
    else
    {
        return 0;
    }
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setComplement (T const& x, const Box& b, int ns, int num) noexcept
{
    this->setComplement<run_on>(x, b, DestComp{ns}, NumComps{num});
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::abs () noexcept
{
    this->abs<run_on>(this->domain,0,this->nvar);
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::abs (int comp, int numcomp) noexcept
{
    this->abs<run_on>(this->domain,comp,numcomp);
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::abs (const Box& subbox, int comp, int numcomp) noexcept
{
    Array4<T> const& a = this->array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG (run_on, subbox, numcomp, i, j, k, n,
    {
        a(i,j,k,n+comp) = amrex::Math::abs(a(i,j,k,n+comp));
    });
}

template <class T>
template <RunOn run_on>
Real
BaseFab<T>::norminfmask (const Box& subbox, const BaseFab<int>& mask,
                         int scomp, int ncomp) const noexcept
{
    BL_ASSERT(this->domain.contains(subbox));
    BL_ASSERT(scomp >= 0 && scomp + ncomp <= this->nvar);

    Array4<T const> const& a = this->const_array();
    Array4<int const> const& m = mask.const_array();
    Real r = 0.0;
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = ReduceData<Real>::Type;
        reduce_op.eval(subbox, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real t = 0.0;
            if (m(i,j,k)) {
                for (int n = 0; n < ncomp; ++n) {
                    t = amrex::max(t,static_cast<Real>(amrex::Math::abs(a(i,j,k,n+scomp))));
                }
            }
            return {t};
        });
        ReduceTuple hv = reduce_data.value();
        r = amrex::get<0>(hv);
    } else
#endif
    {
        amrex::LoopOnCpu(subbox, ncomp, [=,&r] (int i, int j, int k, int n) noexcept
        {
            if (m(i,j,k)) {
                Real t = static_cast<Real>(amrex::Math::abs(a(i,j,k,n+scomp)));
                r = amrex::max(r,t);
            }
        });
    }
    return r;
}

template <class T>
template <RunOn run_on>
Real
BaseFab<T>::norm (int p, int comp, int numcomp) const noexcept
{
    return norm<run_on>(this->domain,p,comp,numcomp);
}

template <class T>
template <RunOn run_on>
Real
BaseFab<T>::norm (const Box& subbox, int p, int comp, int numcomp) const noexcept
{
    BL_ASSERT(this->domain.contains(subbox));
    BL_ASSERT(comp >= 0 && comp + numcomp <= this->nvar);

    Array4<T const> const& a = this->const_array();
    Real nrm = 0.;
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        if (p == 0) {
            ReduceOps<ReduceOpMax> reduce_op;
            ReduceData<Real> reduce_data(reduce_op);
            using ReduceTuple = ReduceData<Real>::Type;
            reduce_op.eval(subbox, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                Real t = 0.0;
                for (int n = 0; n < numcomp; ++n) {
                    t = amrex::max(t, static_cast<Real>(amrex::Math::abs(a(i,j,k,n+comp))));
                }
                return {t};
            });
            ReduceTuple hv = reduce_data.value();
            nrm = amrex::get<0>(hv);
        } else if (p == 1) {
            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<Real> reduce_data(reduce_op);
            using ReduceTuple = ReduceData<Real>::Type;
            reduce_op.eval(subbox, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                Real t = 0.0;
                for (int n = 0; n < numcomp; ++n) {
                    t += static_cast<Real>(amrex::Math::abs(a(i,j,k,n+comp)));
                }
                return {t};
            });
            ReduceTuple hv = reduce_data.value();
            nrm = amrex::get<0>(hv);
        } else if (p == 2) {
            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<Real> reduce_data(reduce_op);
            using ReduceTuple = ReduceData<Real>::Type;
            reduce_op.eval(subbox, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                Real t = 0.0;
                for (int n = 0; n < numcomp; ++n) {
                    t += static_cast<Real>(a(i,j,k,n+comp)*a(i,j,k,n+comp));
                }
                return {t};
            });
            ReduceTuple hv = reduce_data.value();
            nrm = amrex::get<0>(hv);
        } else {
            amrex::Error("BaseFab<T>::norm: wrong p");
        }
    } else
#endif
    {
        if (p == 0) {
            amrex::LoopOnCpu(subbox, numcomp, [=,&nrm] (int i, int j, int k, int n) noexcept
            {
                Real t = static_cast<Real>(amrex::Math::abs(a(i,j,k,n+comp)));
                nrm = amrex::max(nrm,t);
            });
        } else if (p == 1) {
            amrex::LoopOnCpu(subbox, numcomp, [=,&nrm] (int i, int j, int k, int n) noexcept
            {
                nrm += amrex::Math::abs(a(i,j,k,n+comp));
            });
        } else if (p == 2) {
            amrex::LoopOnCpu(subbox, numcomp, [=,&nrm] (int i, int j, int k, int n) noexcept
            {
                nrm += a(i,j,k,n+comp)*a(i,j,k,n+comp);
            });
        } else {
            amrex::Error("BaseFab<T>::norm: wrong p");
        }
    }

    return nrm;
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::min (int comp) const noexcept
{
    return this->min<run_on>(this->domain,comp);
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::min (const Box& subbox, int comp) const noexcept
{
    Array4<T const> const& a = this->const_array(comp);
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpMin> reduce_op;
        ReduceData<T> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(subbox, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            return { a(i,j,k) };
        });
        ReduceTuple hv = reduce_data.value();
        return amrex::get<0>(hv);
    } else
#endif
    {
        T r = std::numeric_limits<T>::max();
        amrex::LoopOnCpu(subbox, [=,&r] (int i, int j, int k) noexcept
        {
            r = amrex::min(r, a(i,j,k));
        });
        return r;
    }
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::max (int comp) const noexcept
{
    return this->max<run_on>(this->domain,comp);
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::max (const Box& subbox, int comp) const noexcept
{
    Array4<T const> const& a = this->const_array(comp);
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<T> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(subbox, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            return { a(i,j,k) };
        });
        ReduceTuple hv = reduce_data.value();
        return amrex::get<0>(hv);
    } else
#endif
    {
        T r = std::numeric_limits<T>::lowest();
        amrex::LoopOnCpu(subbox, [=,&r] (int i, int j, int k) noexcept
        {
            r = amrex::max(r, a(i,j,k));
        });
        return r;
    }
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::maxabs (int comp) const noexcept
{
    return this->maxabs<run_on>(this->domain,comp);
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::maxabs (const Box& subbox, int comp) const noexcept
{
    Array4<T const> const& a = this->const_array(comp);
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<T> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(subbox, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            return { amrex::Math::abs(a(i,j,k)) };
        });
        ReduceTuple hv = reduce_data.value();
        return amrex::get<0>(hv);
    } else
#endif
    {
        T r = 0;
        amrex::Loop(subbox, [=,&r] (int i, int j, int k) noexcept
        {
            r = amrex::max(r, amrex::Math::abs(a(i,j,k)));
        });
        return r;
    }
}


template <class T>
template <RunOn run_on>
IntVect
BaseFab<T>::indexFromValue (Box const& subbox, int comp, T const& value) const noexcept
{
    Array4<T const> const& a = this->const_array(comp);
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        Array<int,1+AMREX_SPACEDIM> ha{0,AMREX_D_DECL(std::numeric_limits<int>::lowest(),
                                                      std::numeric_limits<int>::lowest(),
                                                      std::numeric_limits<int>::lowest())};
        Gpu::DeviceVector<int> dv(1+AMREX_SPACEDIM);
        int* p = dv.data();
        Gpu::htod_memcpy_async(p, ha.data(), sizeof(int)*ha.size());
        amrex::ParallelFor(subbox, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int* flag = p;
            if ((*flag == 0) and (a(i,j,k) == value)) {
                if (Gpu::Atomic::Exch(flag,1) == 0) {
                    AMREX_D_TERM(p[1] = i;,
                                 p[2] = j;,
                                 p[3] = k;);
                }
            }
        });
        Gpu::dtoh_memcpy_async(ha.data(), p, sizeof(int)*ha.size());
        Gpu::streamSynchronize();
        return IntVect(AMREX_D_DECL(ha[1],ha[2],ha[2]));
    } else
#endif
    {
        AMREX_LOOP_3D(subbox, i, j, k,
        {
            if (a(i,j,k) == value) return IntVect(AMREX_D_DECL(i,j,k));
        });
        return IntVect::TheMinVector();
    }
}

template <class T>
template <RunOn run_on>
IntVect
BaseFab<T>::minIndex (int comp) const noexcept
{
    return this->minIndex<run_on>(this->domain,comp);
}

template <class T>
template <RunOn run_on>
IntVect
BaseFab<T>::minIndex (const Box& subbox, int comp) const noexcept
{
    T min_val = this->min<run_on>(subbox, comp);
    return this->indexFromValue<run_on>(subbox, comp, min_val);
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::minIndex (const Box& subbox, Real& min_val, IntVect& min_idx, int comp) const noexcept
{
    min_val = this->min<run_on>(subbox, comp);
    min_idx = this->indexFromValue<run_on>(subbox, comp, min_val);
}

template <class T>
template <RunOn run_on>
IntVect
BaseFab<T>::maxIndex (int comp) const noexcept
{
    return this->maxIndex<run_on>(this->domain,comp);
}

template <class T>
template <RunOn run_on>
IntVect
BaseFab<T>::maxIndex (const Box& subbox, int comp) const noexcept
{
    T max_val = this->max<run_on>(subbox, comp);
    return this->indexFromValue<run_on>(subbox, comp, max_val);
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::maxIndex (const Box& subbox, Real& max_val, IntVect& max_idx, int comp) const noexcept
{
    max_val = this->max<run_on>(subbox, comp);
    max_idx = this->indexFromValue<run_on>(subbox, comp, max_val);
}

template <class T>
template <RunOn run_on>
int
BaseFab<T>::maskLT (BaseFab<int>& mask, T const& val, int comp) const noexcept
{
    mask.resize(this->domain,1);
    int cnt = 0;
    Array4<int> const& m = mask.array();
    Array4<T const> const& a = this->const_array(comp);
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(this->domain, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            int t;
            if (a(i,j,k) < val) {
                m(i,j,k) = 1;
                t = 1;
            } else {
                m(i,j,k) = 0;
                t = 0;
            }
            return {t};
        });
        ReduceTuple hv = reduce_data.value();
        cnt = amrex::get<0>(hv);
    } else
#endif
    {
        AMREX_LOOP_3D(this->domain, i, j, k,
        {
            if (a(i,j,k) < val) {
                m(i,j,k) = 1;
                ++cnt;
            } else {
                m(i,j,k) = 0;
            }
        });
    }

    return cnt;
}

template <class T>
template <RunOn run_on>
int
BaseFab<T>::maskLE (BaseFab<int>& mask, T const& val, int comp) const noexcept
{
    mask.resize(this->domain,1);
    int cnt = 0;
    Array4<int> const& m = mask.array();
    Array4<T const> const& a = this->const_array(comp);
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(this->domain, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            int t;
            if (a(i,j,k) <= val) {
                m(i,j,k) = 1;
                t = 1;
            } else {
                m(i,j,k) = 0;
                t = 0;
            }
            return {t};
        });
        ReduceTuple hv = reduce_data.value();
        cnt = amrex::get<0>(hv);
    } else
#endif
    {
        AMREX_LOOP_3D(this->domain, i, j, k,
        {
            if (a(i,j,k) <= val) {
                m(i,j,k) = 1;
                ++cnt;
            } else {
                m(i,j,k) = 0;
            }
        });
    }

    return cnt;
}

template <class T>
template <RunOn run_on>
int
BaseFab<T>::maskEQ (BaseFab<int>& mask, T const& val, int comp) const noexcept
{
    mask.resize(this->domain,1);
    int cnt = 0;
    Array4<int> const& m = mask.array();
    Array4<T const> const& a = this->const_array(comp);
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(this->domain, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            int t;
            if (a(i,j,k) == val) {
                m(i,j,k) = 1;
                t = 1;
            } else {
                m(i,j,k) = 0;
                t = 0;
            }
            return {t};
        });
        ReduceTuple hv = reduce_data.value();
        cnt = amrex::get<0>(hv);
    } else
#endif
    {
        AMREX_LOOP_3D(this->domain, i, j, k,
        {
            if (a(i,j,k) == val) {
                m(i,j,k) = 1;
                ++cnt;
            } else {
                m(i,j,k) = 0;
            }
        });
    }

    return cnt;
}

template <class T>
template <RunOn run_on>
int
BaseFab<T>::maskGT (BaseFab<int>& mask, T const& val, int comp) const noexcept
{
    mask.resize(this->domain,1);
    int cnt = 0;
    Array4<int> const& m = mask.array();
    Array4<T const> const& a = this->const_array(comp);
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(this->domain, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            int t;
            if (a(i,j,k) > val) {
                m(i,j,k) = 1;
                t = 1;
            } else {
                m(i,j,k) = 0;
                t = 0;
            }
            return {t};
        });
        ReduceTuple hv = reduce_data.value();
        cnt = amrex::get<0>(hv);
    } else
#endif
    {
        AMREX_LOOP_3D(this->domain, i, j, k,
        {
            if (a(i,j,k) > val) {
                m(i,j,k) = 1;
                ++cnt;
            } else {
                m(i,j,k) = 0;
            }
        });
    }

    return cnt;
}

template <class T>
template <RunOn run_on>
int
BaseFab<T>::maskGE (BaseFab<int>& mask, T const& val, int comp) const noexcept
{
    mask.resize(this->domain,1);
    int cnt = 0;
    Array4<int> const& m = mask.array();
    Array4<T const> const& a = this->const_array(comp);
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(this->domain, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            int t;
            if (a(i,j,k) >= val) {
                m(i,j,k) = 1;
                t = 1;
            } else {
                m(i,j,k) = 0;
                t = 0;
            }
            return {t};
        });
        ReduceTuple hv = reduce_data.value();
        cnt = amrex::get<0>(hv);
    } else
#endif
    {
        AMREX_LOOP_3D(this->domain, i, j, k,
        {
            if (a(i,j,k) >= val) {
                m(i,j,k) = 1;
                ++cnt;
            } else {
                m(i,j,k) = 0;
            }
        });
    }

    return cnt;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::atomicAdd (const BaseFab<T>& x) noexcept
{
    Box ovlp(this->domain);
    ovlp &= x.domain;
    return ovlp.ok() ? this->atomicAdd<run_on>(x,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::saxpy (T a, const BaseFab<T>& x, const Box& srcbox, const Box& destbox,
                   int srccomp, int destcomp, int numcomp) noexcept
{
    BL_ASSERT(srcbox.ok());
    BL_ASSERT(x.box().contains(srcbox));
    BL_ASSERT(destbox.ok());
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT( srccomp >= 0 &&  srccomp+numcomp <= x.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <=   nComp());

    Array4<T> const& d = this->array();
    Array4<T const> const& s = x.const_array();
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        d(i,j,k,n+destcomp) += a * s(i+offset.x,j+offset.y,k+offset.z,n+srccomp);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::saxpy (T a, const BaseFab<T>& x) noexcept
{
    Box ovlp(this->domain);
    ovlp &= x.domain;
    return ovlp.ok() ? saxpy<run_on>(a,x,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::xpay (T a, const BaseFab<T>& x,
                  const Box& srcbox, const Box& destbox,
                  int srccomp, int destcomp, int numcomp) noexcept
{
    BL_ASSERT(srcbox.ok());
    BL_ASSERT(x.box().contains(srcbox));
    BL_ASSERT(destbox.ok());
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT( srccomp >= 0 &&  srccomp+numcomp <= x.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <=   nComp());

    Array4<T> const& d = this->array();
    Array4<T const> const& s = x.const_array();
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        d(i,j,k,n+destcomp) = s(i+offset.x,j+offset.y,k+offset.z,n+srccomp) + a*d(i,j,k,n+destcomp);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::addproduct (const Box& destbox, int destcomp, int numcomp,
                        const BaseFab<T>& src1, int comp1,
                        const BaseFab<T>& src2, int comp2) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(   comp1 >= 0 &&    comp1+numcomp <= src1.nComp());
    BL_ASSERT(   comp2 >= 0 &&    comp2+numcomp <= src2.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <=      nComp());

    Array4<T> const& d = this->array();
    Array4<T const> const& s1 = src1.const_array();
    Array4<T const> const& s2 = src2.const_array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        d(i,j,k,n+destcomp) += s1(i,j,k,n+comp1) * s2(i,j,k,n+comp2);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::linComb (const BaseFab<T>& f1, const Box& b1, int comp1,
                     const BaseFab<T>& f2, const Box& b2, int comp2,
                     Real alpha, Real beta, const Box& b,
                     int comp, int numcomp) noexcept
{
    BL_ASSERT(b1.ok());
    BL_ASSERT(f1.box().contains(b1));
    BL_ASSERT(b2.ok());
    BL_ASSERT(f2.box().contains(b2));
    BL_ASSERT(b.ok());
    BL_ASSERT(box().contains(b));
    BL_ASSERT(b.sameSize(b1));
    BL_ASSERT(b.sameSize(b2));
    BL_ASSERT(comp1 >= 0 && comp1+numcomp <= f1.nComp());
    BL_ASSERT(comp2 >= 0 && comp2+numcomp <= f2.nComp());
    BL_ASSERT(comp  >= 0 && comp +numcomp <=    nComp());

    Array4<T> const& d = this->array();
    Array4<T const> const& s1 = f1.const_array();
    Array4<T const> const& s2 = f2.const_array();
    const auto dlo = amrex::lbound(b);
    const auto slo1 = amrex::lbound(b1);
    const auto slo2 = amrex::lbound(b2);
    const Dim3 off1{slo1.x-dlo.x,slo1.y-dlo.y,slo1.z-dlo.z};
    const Dim3 off2{slo2.x-dlo.x,slo2.y-dlo.y,slo2.z-dlo.z};    

    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, b, numcomp, i, j, k, n,
    {
        d(i,j,k,n+comp) = alpha*s1(i+off1.x,j+off1.y,k+off1.z,n+comp1)
                         + beta*s2(i+off2.x,j+off2.y,k+off2.z,n+comp2);
    });
    return *this;
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::dot (const Box& xbx, int xcomp,
                 const BaseFab<T>& y, const Box& ybx, int ycomp,
                 int numcomp) const noexcept
{
    BL_ASSERT(xbx.ok());
    BL_ASSERT(box().contains(xbx));
    BL_ASSERT(y.box().contains(ybx));
    BL_ASSERT(xbx.sameSize(ybx));
    BL_ASSERT(xcomp >= 0 && xcomp+numcomp <=   nComp());
    BL_ASSERT(ycomp >= 0 && ycomp+numcomp <= y.nComp());

    T r = 0;

    const auto xlo = amrex::lbound(xbx);
    const auto ylo = amrex::lbound(ybx);
    const Dim3 offset{ylo.x-xlo.x,ylo.y-xlo.y,ylo.z-xlo.z};
    Array4<T const> const& xa = this->const_array();
    Array4<T const> const& ya = y.const_array();

#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<T> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(xbx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            T t = 0;
            for (int n = 0; n < numcomp; ++n) {
                t += xa(i,j,k,n+xcomp) * ya(i+offset.x,j+offset.y,k+offset.z,n+ycomp);
            }
            return {t};
        });
        ReduceTuple hv = reduce_data.value();
        r = amrex::get<0>(hv);
    } else
#endif
    {
        AMREX_LOOP_4D(xbx, numcomp, i, j, k, n,
        {
            r += xa(i,j,k,n+xcomp) * ya(i+offset.x,j+offset.y,k+offset.z,n+ycomp);
        });
    }

    return r;
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::dotmask (const BaseFab<int>& mask, const Box& xbx, int xcomp,
                     const BaseFab<T>& y, const Box& ybx, int ycomp,
                     int numcomp) const noexcept
{
    BL_ASSERT(xbx.ok());
    BL_ASSERT(box().contains(xbx));
    BL_ASSERT(y.box().contains(ybx));
    BL_ASSERT(xbx.sameSize(ybx));
    BL_ASSERT(xcomp >= 0 && xcomp+numcomp <=   nComp());
    BL_ASSERT(ycomp >= 0 && ycomp+numcomp <= y.nComp());

    T r = 0;

    const auto xlo = amrex::lbound(xbx);
    const auto ylo = amrex::lbound(ybx);
    const Dim3 offset{ylo.x-xlo.x,ylo.y-xlo.y,ylo.z-xlo.z};

    Array4<T const> const& xa = this->const_array();
    Array4<T const> const& ya = y.const_array();
    Array4<int const> const& ma = mask.const_array();

#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<T> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(xbx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            int m = static_cast<int>(static_cast<bool>(ma(i,j,k)));
            T t = 0;
            for (int n = 0; n < numcomp; ++n) {
                t += xa(i,j,k,n+xcomp) * ya(i+offset.x,j+offset.y,k+offset.z,n+ycomp) * m;
            }
            return {t};
        });
        ReduceTuple hv = reduce_data.value();
        r = amrex::get<0>(hv);
    } else
#endif
    {
        AMREX_LOOP_4D(xbx, numcomp, i, j, k, n,
        {
            int m = static_cast<int>(static_cast<bool>(ma(i,j,k)));
            r += xa(i,j,k,n+xcomp) * ya(i+offset.x,j+offset.y,k+offset.z,n+ycomp) * m;
        });
    }

    return r;
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::sum (int comp, int numcomp) const noexcept
{
    return this->sum<run_on>(this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::sum (const Box& subbox, int comp, int numcomp) const noexcept
{
    return this->sum<run_on>(subbox, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::negate (int comp, int numcomp) noexcept
{
    return this->negate<run_on>(this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::negate (const Box& b, int comp, int numcomp) noexcept
{
    return this->negate<run_on>(b, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::invert (T const& r, int comp, int numcomp) noexcept
{
    return this->invert<run_on>(r, this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::invert (T const& r, const Box& b, int comp, int numcomp) noexcept
{
    return this->invert<run_on>(r, b, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::plus (T const& r, int comp, int numcomp) noexcept
{
    return this->plus<run_on>(r, this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::plus (T const& r, const Box& b, int comp, int numcomp) noexcept
{
    return this->plus<run_on>(r, b, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->plus<run_on>(src, this->domain, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::atomicAdd (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->atomicAdd<run_on>(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp,
                  int numcomp) noexcept
{
    return this->plus<run_on>(src, subbox, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::atomicAdd (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp,
                       int numcomp) noexcept
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    ovlp &= subbox;
    return ovlp.ok() ? this->atomicAdd<run_on>(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                  int srccomp, int destcomp, int numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        d(i,j,k,n+destcomp) += s(i+offset.x,j+offset.y,k+offset.z,n+srccomp);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::atomicAdd (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                       int srccomp, int destcomp, int numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        T* p = d.ptr(i,j,k,n+destcomp);
        HostDevice::Atomic::Add(p, s(i+offset.x,j+offset.y,k+offset.z,n+srccomp));
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->minus<run_on>(src, this->domain, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->minus<run_on>(src, subbox, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                   int srccomp, int destcomp, int numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        d(i,j,k,n+destcomp) -= s(i+offset.x,j+offset.y,k+offset.z,n+srccomp);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::mult (T const& r, int comp, int numcomp) noexcept
{
    return this->mult<run_on>(r, this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::mult (T const& r, const Box& b, int comp, int numcomp) noexcept
{
    return this->mult<run_on>(r, b, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->mult<run_on>(src, this->domain, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->mult<run_on>(src, subbox, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                  int srccomp, int destcomp, int numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        d(i,j,k,n+destcomp) *= s(i+offset.x,j+offset.y,k+offset.z,n+srccomp);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::divide (T const& r, int comp, int numcomp) noexcept
{
    return this->divide<run_on>(r, this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::divide (T const& r, const Box& b, int comp, int numcomp) noexcept
{
    return this->divide<run_on>(r, b, DestComp{comp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->divide<run_on>(src, this->domain, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->divide<run_on>(src, subbox, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                    int srccomp, int destcomp, int numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        d(i,j,k,n+destcomp) /= s(i+offset.x,j+offset.y,k+offset.z,n+srccomp);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::protected_divide (const BaseFab<T>& src) noexcept
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->protected_divide<run_on>(src,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::protected_divide (const BaseFab<T>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->protected_divide<run_on>(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::protected_divide (const BaseFab<T>& src, const Box& subbox, int srccomp, int destcomp,
                              int numcomp) noexcept
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    ovlp &= subbox;
    return ovlp.ok() ? this->protected_divide<run_on>(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::protected_divide (const BaseFab<T>& src, const Box& srcbox, const Box& destbox,
                              int srccomp, int destcomp, int numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const Dim3 offset{slo.x-dlo.x,slo.y-dlo.y,slo.z-dlo.z};
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, destbox, numcomp, i, j, k, n,
    {
        if (s(i+offset.x,j+offset.y,k+offset.z,n+srccomp) != 0) {
            d(i,j,k,n+destcomp) /= s(i+offset.x,j+offset.y,k+offset.z,n+srccomp);
        }
    });

    return *this;
}

/**
* Linear Interpolation / Extrapolation
* Result is (t2-t)/(t2-t1)*f1 + (t-t1)/(t2-t1)*f2
* Data is taken from b1 region of f1, b2 region of f2
* and stored in b region of this FAB.
* Boxes b, b1 and b2 must be the same size.
* Data is taken from component comp1 of f1, comp2 of f2,
* and stored in component comp of this FAB.
* This fab is returned as a reference for chaining.
*/

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::linInterp (const BaseFab<T>& f1, const Box& b1, int comp1,
                       const BaseFab<T>& f2, const Box& b2, int comp2,
                       Real t1, Real t2, Real t,
                       const Box& b, int comp, int numcomp) noexcept
{
    Real alpha = (t2-t)/(t2-t1);
    Real beta = (t-t1)/(t2-t1);
    return linComb<run_on>(f1,b1,comp1,f2,b2,comp2,alpha,beta,b,comp,numcomp);
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::linInterp (const BaseFab<T>& f1, int comp1,
                       const BaseFab<T>& f2, int comp2,
                       Real t1, Real t2, Real t,
                       const Box& b, int comp, int numcomp) noexcept
{
    Real tol = 1.0e-16;
    if(amrex::Math::abs(t2-t1) > tol)
    {
        Real alpha = (t2-t)/(t2-t1);
        Real beta = (t-t1)/(t2-t1);

        return linComb<run_on>(f1,b,comp1,f2,b,comp2,alpha,beta,b,comp,numcomp);
    }
    else
    {
        copy<run_on>(f1,b,comp1,b,comp,numcomp);
    }
    return *this;
}

//
// New interfaces
//

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setVal (T const& val) noexcept
{
    this->setVal<run_on>(val, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setVal (T const& x, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);
    Array4<T> const& a = this->array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG (run_on, bx, ncomp.n, i, j, k, n,
    {
        a(i,j,k,n+dcomp.i) = x;
    });
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setValIf (T const& val, const BaseFab<int>& mask) noexcept
{
    this->setValIf<run_on>(val, this->domain, mask, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setValIf (T const& val, Box const& bx, const BaseFab<int>& mask, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);
    Array4<T> const& a = this->array();
    Array4<int const> const& m = mask.const_array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG (run_on, bx, ncomp.n, i, j, k, n,
    {
        if (m(i,j,k)) { a(i,j,k,n+dcomp.i) = val; }
    });
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setValIfNot (T const& val, const BaseFab<int>& mask) noexcept
{
    this->setValIfNot<run_on>(val, this->domain, mask, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setValIfNot (T const& val, Box const& bx, const BaseFab<int>& mask, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);
    Array4<T> const& a = this->array();
    Array4<int const> const& m = mask.const_array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG (run_on, bx, ncomp.n, i, j, k, n,
    {
        if (!m(i,j,k)) { a(i,j,k,n+dcomp.i) = val; }
    });
}

template <class T>
template <RunOn run_on>
void
BaseFab<T>::setComplement (T const& x, const Box& bx, DestComp dcomp, NumComps ncomp) noexcept
{
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        Array4<T> const& a = this->array();
        amrex::ParallelFor(this->domain, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (not bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                for (int n = dcomp.i; n < ncomp.n+dcomp.i; ++n) {
                    a(i,j,k,n) = x;
                }
            }
        });
    } else
#endif
    {
        const BoxList b_lst = amrex::boxDiff(this->domain,bx);
        for (auto const& b : b_lst) {
            this->setVal<RunOn::Host>(x, b, dcomp, ncomp);
        }
    }
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src) noexcept
{
    this->copy<run_on>(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src, Box bx,
                  SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    bx &= src.domain;

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, bx, ncomp.n, i, j, k, n,
    {
        d(i,j,k,n+dcomp.i) = s(i,j,k,n+scomp.i);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::plus (T const& val) noexcept
{
    return this->plus<run_on>(val, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::operator+= (T const& val) noexcept
{
    return this->plus<run_on>(val, this->domain,  DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::plus (T const& val, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    Array4<T> const& a = this->array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, bx, ncomp.n, i, j, k, n,
    {
        a(i,j,k,n+dcomp.i) += val;
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src) noexcept
{
    return this->plus<run_on>(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::operator+= (const BaseFab<T>& src) noexcept
{
    return this->plus<run_on>(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src, Box bx,
                  SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    bx &= src.domain;

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, bx, ncomp.n, i, j, k, n,
    {
        d(i,j,k,n+dcomp.i) += s(i,j,k,n+scomp.i);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::minus (T const& val) noexcept
{
    return this->minus<run_on>(val, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::operator-= (T const& val) noexcept
{
    return this->minus<run_on>(val, this->domain,  DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::minus (T const& val, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    Array4<T> const& a = this->array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, bx, ncomp.n, i, j, k, n,
    {
        a(i,j,k,n+dcomp.i) -= val;
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src) noexcept
{
    return this->minus<run_on>(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::operator-= (const BaseFab<T>& src) noexcept
{
    return this->minus<run_on>(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src, Box bx,
                   SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    bx &= src.domain;

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, bx, ncomp.n, i, j, k, n,
    {
        d(i,j,k,n+dcomp.i) -= s(i,j,k,n+scomp.i);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::mult (T const& val) noexcept
{
    return this->mult<run_on>(val, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::operator*= (T const& val) noexcept
{
    return this->mult<run_on>(val, this->domain,  DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::mult (T const& val, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    Array4<T> const& a = this->array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, bx, ncomp.n, i, j, k, n,
    {
        a(i,j,k,n+dcomp.i) *= val;
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src) noexcept
{
    return this->mult<run_on>(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::operator*= (const BaseFab<T>& src) noexcept
{
    return this->mult<run_on>(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src, Box bx,
                  SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    bx &= src.domain;

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, bx, ncomp.n, i, j, k, n,
    {
        d(i,j,k,n+dcomp.i) *= s(i,j,k,n+scomp.i);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::divide (T const& val) noexcept
{
    return this->divide<run_on>(val, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::operator/= (T const& val) noexcept
{
    return this->divide<run_on>(val, this->domain,  DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::divide (T const& val, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    Array4<T> const& a = this->array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, bx, ncomp.n, i, j, k, n,
    {
        a(i,j,k,n+dcomp.i) /= val;
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src) noexcept
{
    return this->divide<run_on>(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::operator/= (const BaseFab<T>& src) noexcept
{
    return this->divide<run_on>(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src, Box bx,
                    SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    bx &= src.domain;

    Array4<T> const& d = this->array();
    Array4<T const> const& s = src.const_array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, bx, ncomp.n, i, j, k, n,
    {
        d(i,j,k,n+dcomp.i) /= s(i,j,k,n+scomp.i);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::negate () noexcept
{
    return this->negate<run_on>(this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::negate (const Box& bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    Array4<T> const& a = this->array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, bx, ncomp.n, i, j, k, n,
    {
        a(i,j,k,n+dcomp.i) = -a(i,j,k,n+dcomp.i);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::invert (T const& r) noexcept
{
    return this->invert<run_on>(r, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T>
template <RunOn run_on>
BaseFab<T>&
BaseFab<T>::invert (T const& r, const Box& bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    Array4<T> const& a = this->array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(run_on, bx, ncomp.n, i, j, k, n,
    {
        a(i,j,k,n+dcomp.i) = r / a(i,j,k,n+dcomp.i);
    });

    return *this;
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::sum (const Box& bx, DestComp dcomp, NumComps ncomp) const noexcept
{
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    T r = 0;
    Array4<T const> const& a = this->const_array();
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<T> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            T t = 0;
            for (int n = 0; n < ncomp.n; ++n) {
                t += a(i,j,k,n+dcomp.i);
            }
            return { t };
        });
        ReduceTuple hv = reduce_data.value();
        r = amrex::get<0>(hv);
    } else
#endif
    {
        amrex::LoopOnCpu(bx, ncomp.n, [=,&r] (int i, int j, int k, int n) noexcept
        {
            r += a(i,j,k,n+dcomp.i);
        });
    }

    return r;
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::dot (const BaseFab<T>& src, const Box& bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) const noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    T r = 0;
    Array4<T const> const& d = this->const_array();
    Array4<T const> const& s = src.const_array();
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<T> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            T t = 0;
            for (int n = 0; n < ncomp.n; ++n) {
                t += d(i,j,k,n+dcomp.i) * s(i,j,k,n+scomp.i);
            }
            return { t };
        });
        ReduceTuple hv = reduce_data.value();
        r = amrex::get<0>(hv);
    } else
#endif
    {
        amrex::LoopOnCpu(bx, ncomp.n, [=,&r] (int i, int j, int k, int n) noexcept
        {
            r += d(i,j,k,n+dcomp.i) * s(i,j,k,n+scomp.i);
        });
    }

    return r;
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::dot (const Box& bx, int destcomp, int numcomp) const noexcept
{
    return dot<run_on>(bx, DestComp{destcomp}, NumComps{numcomp});
}


template <class T>
template <RunOn run_on>
T
BaseFab<T>::dot (const Box& bx, DestComp dcomp, NumComps ncomp) const noexcept
{
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    T r = 0;
    Array4<T const> const& a = this->const_array();
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<T> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            T t = 0;
            for (int n = 0; n < ncomp.n; ++n) {
                t += a(i,j,k,n+dcomp.i)*a(i,j,k,n+dcomp.i);
            }
            return { t };
        });
        ReduceTuple hv = reduce_data.value();
        r = amrex::get<0>(hv);
    } else
#endif
    {
        amrex::LoopOnCpu(bx, ncomp.n, [=,&r] (int i, int j, int k, int n) noexcept
        {
            r += a(i,j,k,n+dcomp.i)*a(i,j,k,n+dcomp.i);
        });
    }

    return r;
}

template <class T>
template <RunOn run_on>
T
BaseFab<T>::dotmask (const BaseFab<T>& src, const Box& bx, const BaseFab<int>& mask,
                     SrcComp scomp, DestComp dcomp, NumComps ncomp) const noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(this->domain.sameType(mask.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    T r = 0;
    Array4<T const> const& d = this->const_array();
    Array4<T const> const& s = src.const_array();
    Array4<int const> const& m = mask.const_array();
#ifdef AMREX_USE_GPU
    if (run_on == RunOn::Device and Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<T> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            T t = 0;
            T mi = static_cast<T>(static_cast<int>(static_cast<bool>(m(i,j,k))));
            for (int n = 0; n < ncomp.n; ++n) {
                t += d(i,j,k,n+dcomp.i)*s(i,j,k,n+scomp.i)*mi;
            }
            return { t };
        });
        ReduceTuple hv = reduce_data.value();
        r = amrex::get<0>(hv);
    } else
#endif
    {
        amrex::LoopOnCpu(bx, ncomp.n, [=,&r] (int i, int j, int k, int n) noexcept
        {
            int mi = static_cast<int>(static_cast<bool>(m(i,j,k)));
            r += d(i,j,k,n+dcomp.i)*s(i,j,k,n+scomp.i)*mi;
        });
    }

    return r;
}

}

#endif /*BL_BASEFAB_H*/
