#ifndef AMREX_BASEFAB_H_
#define AMREX_BASEFAB_H_

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <climits>
#include <array>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Extension.H>
#include <AMReX_BLassert.H>
#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_BoxList.H>
#include <AMReX_BArena.H>
#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_MakeType.H>
#include <AMReX_FabAllocator.H>

#include <AMReX_Gpu.H>

#ifdef USE_PERILLA
#include <LocalConnection.H>
#include <RemoteConnection.H>
#endif

namespace amrex
{
    extern long private_total_bytes_allocated_in_fabs;     //!< total bytes at any given time
    extern long private_total_bytes_allocated_in_fabs_hwm; //!< high-water-mark over a given interval
    extern long private_total_cells_allocated_in_fabs;     //!< total cells at any given time
    extern long private_total_cells_allocated_in_fabs_hwm; //!< high-water-mark over a given interval
#ifdef _OPENMP
#pragma omp threadprivate(private_total_bytes_allocated_in_fabs)
#pragma omp threadprivate(private_total_bytes_allocated_in_fabs_hwm)
#pragma omp threadprivate(private_total_cells_allocated_in_fabs)
#pragma omp threadprivate(private_total_cells_allocated_in_fabs_hwm)
#endif

    long TotalBytesAllocatedInFabs () noexcept;
    long TotalBytesAllocatedInFabsHWM () noexcept;
    long TotalCellsAllocatedInFabs () noexcept;
    long TotalCellsAllocatedInFabsHWM () noexcept;
    void ResetTotalBytesAllocatedInFabsHWM () noexcept;
    void update_fab_stats (long n, long s, std::size_t szt) noexcept;

    void BaseFab_Initialize ();
    void BaseFab_Finalize ();

template <typename T>
struct FabView
{
    T* AMREX_RESTRICT p;
    int jstride;
    int kstride;
    long nstride;
#ifdef AMREX_DEBUG
    int iend, jend, kend, nend;
#endif
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE T& operator() (int i, int j, int k, int n=0) const noexcept {
#ifdef AMREX_DEBUG
        AMREX_ASSERT(i<iend && j<jend && k<kend && n<nend);
#endif
        return *(p + (i+j*jstride+k*kstride+n*nstride));
    }
};

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
makeArray4 (T* p, Box const& bx) noexcept
{
    return Array4<T>{p, amrex::begin(bx), amrex::end(bx)};
}

template <typename T>
struct BaseFabData
    : public FabAllocator
{
    AMREX_GPU_HOST_DEVICE
    BaseFabData (T* a_dptr, Box const& a_domain, int a_nvar,
                 long a_truesize, bool a_ptr_owner, bool a_shared_memory) noexcept
        : dptr(a_dptr), domain(a_domain), nvar(a_nvar),
          truesize(a_truesize), ptr_owner(a_ptr_owner), shared_memory(a_shared_memory)
        {}
    BaseFabData () = default;
    T*      dptr;     //!< The data pointer.
    Box     domain;   //!< My index space.
    int     nvar;     //!< Number components.
    long    truesize; //!< nvar*numpts that was allocated on heap.
    bool    ptr_owner;//!< Owner of T*?
    bool    shared_memory;  //!< Is the memory allocated in shared memory?
    void setOwner (bool b) noexcept { ptr_owner = b; }
};

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

template <class T>
struct DataAllocator {
    void* alloc (std::size_t sz) const noexcept { return The_Arena()->alloc(sz); }
    void free (void* pt) const noexcept { The_Arena()->free(pt); }
    Arena* arena () const noexcept { return The_Arena(); }
};

#ifdef AMREX_USE_GPU
template <class T>
struct CpuDataAllocator {
    void* alloc (std::size_t sz) const noexcept { return The_Cpu_Arena()->alloc(sz); }
    void free (void* pt) const noexcept { The_Cpu_Arena()->free(pt); }
    Arena* arena () const noexcept { return The_Cpu_Arena(); }
};
#else
template <class T> using CpuDataAllocator = DataAllocator<T>;
#endif

template <class T, class Allocator = DataAllocator<T> >
class BaseFab : public BaseFabData<T>
{
public:

    template <class U, class A> friend class BaseFab;

    typedef T         value_type;
    typedef Allocator allocator_type;
    //! Construct an empty BaseFab, which must be resized (see BaseFab::resize) before use.
    BaseFab ();

    //!  Make BaseFab with desired domain (box) and number of components.
    explicit BaseFab (const Box& bx,
                      int        n = 1,
		      bool       alloc = true,
		      bool       shared = false);

    BaseFab (const BaseFab<T,Allocator>& rhs, MakeType make_type, int scomp, int ncomp);

#ifdef AMREX_USE_GPU
    BaseFab (const BaseFab<T,Allocator>& rhs, MakeType make_type);
#endif

    /**
     * \brief Create an NON-OWNING BaseFab.  Thus BaseFab is not
     * responsible for memory management.  And it's caller's responsibility that
     * p points to a chunk of memory large enough.
     */
    BaseFab (const Box& bx, int ncomp, T* p);

    //! The destructor deletes the array memory.
    ~BaseFab ();

    BaseFab (const BaseFab<T,Allocator>& rhs) = delete;
    BaseFab<T,Allocator>& operator= (const BaseFab<T,Allocator>& rhs) = delete;
    BaseFab<T,Allocator>& operator= (BaseFab<T,Allocator>&& rhs) = delete;

    BaseFab (BaseFab<T,Allocator>&& rhs) noexcept;

    AMREX_GPU_HOST_DEVICE
    BaseFab& operator= (T) noexcept;

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
    void resize (const Box& b,
                 int        N = 1);

    template <class U=T, class = typename std::enable_if<std::is_trivially_destructible<U>::value>::type >
    Elixir elixir () noexcept;

    /**
     * \brief The function returns the BaseFab to the invalid state.
     * The memory is freed.
     */
    void clear () noexcept;

    //! Returns how many bytes used
    AMREX_GPU_HOST_DEVICE
    std::size_t nBytes () const noexcept { return this->truesize*sizeof(T); }

    static bool preAllocatable () noexcept { return true; }

    static bool isCopyOMPSafe () noexcept { return true; }

    //! Returns bytes used in the Box for those components
    AMREX_GPU_HOST_DEVICE
    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const noexcept
        { return bx.numPts() * sizeof(T) * ncomps; }

    //! Returns the number of components
    AMREX_GPU_HOST_DEVICE
    int nComp () const noexcept { return this->nvar; }

    //! for calls to fortran.
    AMREX_GPU_HOST_DEVICE
    const int* nCompPtr() const noexcept {
        return &(this->nvar);
    }

    //! Returns the number of points
    AMREX_GPU_HOST_DEVICE
    long numPts () const noexcept { return this->domain.numPts(); }

    //! Returns the total number of points of all components
    AMREX_GPU_HOST_DEVICE
    long size () const noexcept { return this->nvar*this->domain.numPts(); }

    //! Returns the domain (box) where the array is defined
    AMREX_GPU_HOST_DEVICE
    const Box& box () const noexcept { return this->domain; }

    /**
    * \brief Returns a pointer to an array of SPACEDIM integers
    * giving the length of the domain in each direction
    */
    AMREX_GPU_HOST_DEVICE
    IntVect length () const noexcept { return this->domain.length(); }

    /**
    * \brief Returns the lower corner of the domain
    * See class Box for analogue.
    */
    AMREX_GPU_HOST_DEVICE
    const IntVect& smallEnd () const noexcept { return this->domain.smallEnd(); }

    //!  Returns the upper corner of the domain.  See class Box for analogue.
    AMREX_GPU_HOST_DEVICE
    const IntVect& bigEnd () const noexcept { return this->domain.bigEnd(); }

    /**
    * \brief Returns the lower corner of the domain.

    *Instead of returning them in the form of INTVECTs, as in smallEnd and
    * bigEnd, it returns the values as a pointer to an array of
    * constant integers.  This is useful when interfacing to
    * Fortran subroutines.
    */
    AMREX_GPU_HOST_DEVICE
    const int* loVect () const noexcept { return this->domain.loVect(); }

    /**
    * \brief Returns the upper corner of the domain.

    *Instead of returning them in the form of INTVECTs, as in smallEnd and
    * bigEnd, it returns the values as a pointer to an array of
    * constant integers.  This is useful when interfacing to
    * Fortran subroutines.
    */
    AMREX_GPU_HOST_DEVICE
    const int* hiVect () const noexcept { return this->domain.hiVect(); }

    /**
    * \brief Returns true if the domain of fab is totally contained within
    * the domain of this BaseFab.
    */
    AMREX_GPU_HOST_DEVICE
    bool contains (const BaseFab<T,Allocator>& fab) const noexcept
    {
        return box().contains(fab.box()) && this->nvar <= fab.nvar;
    }

    /**
    * \brief Returns true if bx is totally contained
    * within the domain of this BaseFab.
    */
    AMREX_GPU_HOST_DEVICE
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
    AMREX_GPU_HOST_DEVICE
    T* dataPtr (int n = 0) noexcept { AMREX_ASSERT(!(this->dptr == 0)); return &(this->dptr[n*this->domain.numPts()]); }

    //! Same as above except works on const FABs.
    AMREX_GPU_HOST_DEVICE
    const T* dataPtr (int n = 0) const noexcept { AMREX_ASSERT(!(this->dptr == 0)); return &(this->dptr[n*this->domain.numPts()]); }

    AMREX_GPU_HOST_DEVICE
    T* dataPtr (const IntVect& iv, int n = 0) noexcept;

    AMREX_GPU_HOST_DEVICE
    const T* dataPtr (const IntVect& iv, int n = 0) const noexcept;

    void setPtr (T* p, long sz) noexcept { AMREX_ASSERT(this->dptr == 0 && this->truesize == 0); this->dptr = p; this->truesize = sz; }

    void prefetchToHost () const noexcept;
    void prefetchToDevice () const noexcept;

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    Array4<T const> array () const noexcept
    {
        return makeArray4<T const>(this->dptr, this->domain);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    Array4<T> array () noexcept
    {
        return makeArray4<T>(this->dptr, this->domain);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    Array4<T const> const_array () const noexcept
    {
        return makeArray4<T const>(this->dptr, this->domain);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T const> view (int n = 0) const noexcept
    {
        return view(this->domain.loVect3d(),n);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T> view (int n = 0) noexcept
    {
        return view(this->domain.loVect3d(),n);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T const> view (const IntVect& iv, int n = 0) const noexcept
    {
#if (AMREX_SPACEDIM == 1)
        return view(GpuArray<int,3>{iv[0],    0,    0}, n);
#elif (AMREX_SPACEDIM == 2)
        return view(GpuArray<int,3>{iv[0],iv[1],    0}, n);
#else
        return view(GpuArray<int,3>{iv[0],iv[1],iv[2]}, n);
#endif
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T> view (const IntVect& iv, int n = 0) noexcept
    {
#if (AMREX_SPACEDIM == 1)
        return view(GpuArray<int,3>{iv[0],    0,    0}, n);
#elif (AMREX_SPACEDIM == 2)
        return view(GpuArray<int,3>{iv[0],iv[1],    0}, n);
#else
        return view(GpuArray<int,3>{iv[0],iv[1],iv[2]}, n);
#endif
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T const> view (const Box& subbox, int n = 0) const noexcept
    {
        return view(subbox.loVect3d(),n);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T> view (const Box& subbox, int n = 0) noexcept
    {
        return view(subbox.loVect3d(),n);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T const> view (const GpuArray<int,3>& slo, int n = 0) const noexcept
    {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        AMREX_ASSERT(this->domain.numPts() < static_cast<long>(std::numeric_limits<int>::max()));
#else
        AMREX_ASSERT(this->domain.numPts() < static_cast<long>(INT_MAX));
#endif
        const auto& len = this->domain.length3d();
        int jstride = len[0];
        int kstride = jstride*len[1];
        long nstride = kstride*len[2];
        const auto& dlo = this->domain.loVect3d();
        T const* p = this->dptr + ((slo[0]-dlo[0])
                                +  (slo[1]-dlo[1])*jstride
                                +  (slo[2]-dlo[2])*kstride
                                +                n*nstride);
#ifdef AMREX_DEBUG
        return FabView<T const>{p, jstride, kstride, nstride,
                len[0]-(slo[0]-dlo[0]),
                len[1]-(slo[1]-dlo[1]),
                len[2]-(slo[2]-dlo[2]),
                this->nvar-n};
#else
        return FabView<T const>{p, jstride, kstride, nstride};
#endif
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T> view (const GpuArray<int,3>& slo, int n = 0) noexcept
    {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        AMREX_ASSERT(this->domain.numPts() < static_cast<long>(std::numeric_limits<int>::max()));
#else
        AMREX_ASSERT(this->domain.numPts() < static_cast<long>(INT_MAX));
#endif
        const auto& len = this->domain.length3d();
        int jstride = len[0];
        int kstride = jstride*len[1];
        long nstride = kstride*len[2];
        const auto& dlo = this->domain.loVect3d();
        T* p = this->dptr + ((slo[0]-dlo[0])
                          +  (slo[1]-dlo[1])*jstride
                          +  (slo[2]-dlo[2])*kstride
                          +                n*nstride);
#ifdef AMREX_DEBUG
        return FabView<T>{p, jstride, kstride, nstride,
                len[0]-(slo[0]-dlo[0]),
                len[1]-(slo[1]-dlo[1]),
                len[2]-(slo[2]-dlo[2]),
                this->nvar-n};
#else
        return FabView<T>{p, jstride, kstride, nstride};
#endif
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T const> view (const Dim3& slo, int n = 0) const noexcept
    {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        AMREX_ASSERT(this->domain.numPts() < static_cast<long>(std::numeric_limits<int>::max()));
#else
        AMREX_ASSERT(this->domain.numPts() < static_cast<long>(INT_MAX));
#endif
        const auto len = amrex::length(this->domain);
        const auto dlo = amrex::lbound(this->domain);
        int jstride = len.x;
        int kstride = jstride*len.y;
        long nstride = kstride*len.z;
        T const* p = this->dptr + ((slo.x-dlo.x)
                                +  (slo.y-dlo.y)*jstride
                                +  (slo.z-dlo.z)*kstride
                                +              n*nstride);
#ifdef AMREX_DEBUG
        return FabView<T const>{p, jstride, kstride, nstride,
                len.x-(slo.x-dlo.x),
                len.y-(slo.y-dlo.y),
                len.z-(slo.z-dlo.z),
                this->nvar-n};
#else
        return FabView<T const>{p, jstride, kstride, nstride};
#endif
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T> view (const Dim3& slo, int n = 0) noexcept
    {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        AMREX_ASSERT(this->domain.numPts() < static_cast<long>(std::numeric_limits<int>::max()));
#else
        AMREX_ASSERT(this->domain.numPts() < static_cast<long>(INT_MAX));
#endif
        const auto len = amrex::length(this->domain);
        const auto dlo = amrex::lbound(this->domain);
        int jstride = len.x;
        int kstride = jstride*len.y;
        long nstride = kstride*len.z;
        T* p = this->dptr + ((slo.x-dlo.x)
                          +  (slo.y-dlo.y)*jstride
                          +  (slo.z-dlo.z)*kstride
                          +              n*nstride);
#ifdef AMREX_DEBUG
        return FabView<T>{p, jstride, kstride, nstride,
                len.x-(slo.x-dlo.x),
                len.y-(slo.y-dlo.y),
                len.z-(slo.z-dlo.z),
                this->nvar-n};
#else
        return FabView<T>{p, jstride, kstride, nstride};
#endif
    }

    //! Returns true if the data for the FAB has been allocated.
    AMREX_GPU_HOST_DEVICE
    bool isAllocated () const noexcept { return this->dptr != 0; }

    /**
    * \brief Returns a reference to the Nth component value
    * defined at position p in the domain.  This operator may be
    * inefficient if the C++ compiler is unable to optimize the
    * C++ code.
    */
    AMREX_GPU_HOST_DEVICE
    T& operator() (const IntVect& p, int N) noexcept;

    //! Same as above, except returns component 0.
    AMREX_GPU_HOST_DEVICE
    T& operator() (const IntVect& p) noexcept;

    //! Same as above except works on const FABs.
    AMREX_GPU_HOST_DEVICE
    const T& operator() (const IntVect& p, int N) const noexcept;

    //! Same as above, except returns component 0.
    AMREX_GPU_HOST_DEVICE
    const T& operator() (const IntVect& p) const noexcept;

    /**
    * \brief This function puts numcomp component values, starting at
    * component N, from position pos in the domain into array data,
    * that must be allocated by the user.
    */
    AMREX_GPU_HOST_DEVICE
    void getVal (T*             data,
                 const IntVect& pos,
                 int            N,
                 int            numcomp) const noexcept;
    //! Same as above, except that starts at component 0 and copies all comps.
    AMREX_GPU_HOST_DEVICE
    void getVal (T*             data,
                 const IntVect& pos) const noexcept;
    /**
    * \brief The setVal functions set sub-regions in the BaseFab to a
    * constant value.  This most general form specifies the sub-box,
    * the starting component number, and the number of components
    * to be set.
    */
    AMREX_GPU_HOST_DEVICE
    void setVal (T x, const Box& bx, int nstart, int ncomp) noexcept;
    //! Same as above, except the number of modified components is one. N is the component to be modified.
    AMREX_GPU_HOST_DEVICE
    void setVal (T x, const Box& bx, int N = 0) noexcept;
    //! Same as above, except the sub-box defaults to the entire domain.
    AMREX_GPU_HOST_DEVICE
    void setVal (T x, int N) noexcept;

    AMREX_GPU_HOST_DEVICE
    void setValIfNot (T x, const Box& bx, const BaseFab<int>& mask, int nstart, int ncomp) noexcept;

    /**
    * \brief This function is analogous to the fourth form of
    * setVal above, except that instead of setting values on the
    * Box b, values are set on the complement of b in the domain.
    */
    void setComplement (T x, const Box& b, int ns, int num) noexcept;

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
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& copy (const BaseFab<T,Allocator>& src,
                                const Box&        srcbox,
                                int               srccomp,
                                const Box&        destbox,
                                int               destcomp,
                                int               numcomp) noexcept;

    /**
    * \brief As above, except the destination Box and the source Box
    * are taken to be the entire domain of the destination.   A copy
    * of the intersecting region is performed.
    * class.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& copy (const BaseFab<T,Allocator>& src,
                                int               srccomp,
                                int               destcomp,
                                int               numcomp = 1) noexcept;
    /**
    * \brief As above, except that the destination Box is specified,
    * but the source Box is taken to the equal to the source
    * Box, and all components of the destination BaseFab are
    * copied.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& copy (const BaseFab<T,Allocator>& src,
                                const Box&        destbox) noexcept;

    //! Copy from the srcbox of this Fab to raw memory and return the number of bytes copied
    AMREX_GPU_HOST_DEVICE
    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const noexcept;

    AMREX_GPU_HOST_DEVICE
    //! Copy from raw memory to the dstbox of this Fab and return the number of bytes copied
    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src) noexcept;

    AMREX_GPU_HOST_DEVICE
    //! Add from raw memory to the dstbox of this Fab and return the number of bytes copied
    std::size_t addFromMem (const Box&  dstbox,
                            int         dstcomp,
                            int         numcomp,
                            const void* src) noexcept;

    /**
    * \brief Perform shifts upon the domain of the BaseFab. They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& shift (const IntVect& v) noexcept;
    /**
    * \brief Perform shifts upon the domain of the BaseFab.  They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& shift (int idir, int n_cell) noexcept;
    /**
    * \brief Perform shifts upon the domain of the BaseFab.  They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& shiftHalf (int dir, int num_halfs) noexcept;
    /**
    * \brief Perform shifts upon the domain of the BaseFab. They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& shiftHalf (const IntVect& num_halfs) noexcept;

    AMREX_GPU_HOST_DEVICE
    Real norminfmask (const Box& subbox, const BaseFab<int>& mask, int scomp=0, int ncomp=1) const noexcept;

    /**
    * \brief Compute the Lp-norm of this FAB using components (scomp : scomp+ncomp-1).
    *   p < 0  -> ERROR
    *   p = 0  -> infinity norm (max norm)
    *   p = 1  -> sum of ABS(FAB)
    */
    AMREX_GPU_HOST_DEVICE
    Real norm (int p,
               int scomp = 0,
               int ncomp = 1) const noexcept;

    //! Same as above except only on given subbox.
    AMREX_GPU_HOST_DEVICE
    Real norm (const Box& subbox,
               int        p,
               int        scomp = 0,
               int        ncomp = 1) const noexcept;
    //!Compute absolute value for all components of this FAB.
    AMREX_GPU_HOST_DEVICE
    void abs () noexcept;
    //! Same as above except only for components (comp: comp+numcomp-1)
    AMREX_GPU_HOST_DEVICE
    void abs (int comp,
              int numcomp=1) noexcept;
    /**
    * \brief Calculate abs() on subbox for given component range.
    */
    AMREX_GPU_HOST_DEVICE
    void abs (const Box& subbox,
              int        comp = 0,
              int        numcomp=1) noexcept;
    /**
    * \return Minimum value of given component.
    */
    AMREX_GPU_HOST_DEVICE
    T min (int comp = 0) const noexcept;
    /**
    * \return Minimum value of given component in given subbox.
    */
    AMREX_GPU_HOST_DEVICE
    T min (const Box& subbox,
           int        comp = 0) const noexcept;
    /**
    * \return Maximum value of given component.
    */
    AMREX_GPU_HOST_DEVICE
    T max (int comp = 0) const noexcept;
    /**
    * \return Maximum value of given component in given subbox.
    */
    AMREX_GPU_HOST_DEVICE
    T max (const Box& subbox,
           int        comp = 0) const noexcept;
    /**
    * \return Maximum of the absolute value of given component.
    */
    AMREX_GPU_HOST_DEVICE
    T maxabs (int comp = 0) const noexcept;
    /**
    * \return Maximum of the absolute value of given component in given subbox.
    */
    AMREX_GPU_HOST_DEVICE
    T maxabs (const Box& subbox,
              int        comp = 0) const noexcept;

    /**
    * \return location of a cell containing the specified value
    * given subbox. Returns IntVect outside box if value not present.
    */
    AMREX_GPU_HOST_DEVICE
    IntVect indexFromValue (Real       value,
                            const Box& subbox,
                            int        comp = 0) const noexcept;

    /**
    * \return location of minimum value in given component.
    */
    AMREX_GPU_HOST_DEVICE
    IntVect minIndex (int comp = 0) const noexcept;
    /**
    * \return location of minimum value in given component in
    * given subbox.
    */
    AMREX_GPU_HOST_DEVICE
    IntVect minIndex (const Box& subbox,
                      int        comp = 0) const noexcept;
    /**
    * \return return mininum value and location to allow
    * efficient looping over multiple boxes.
    */
    AMREX_GPU_HOST_DEVICE
    void  minIndex (const Box&     subbox,
                          Real&    min_val,
                          IntVect& min_idx,
                          int      comp = 0) const noexcept;

    /**
    * \return location of maximum value in given component.
    */
    AMREX_GPU_HOST_DEVICE
    IntVect maxIndex (int comp = 0) const noexcept;
    /**
    * \return location of maximum value in given component in given
    * subbox.
    */
    AMREX_GPU_HOST_DEVICE
    IntVect maxIndex (const Box& subbox,
                      int        comp = 0) const noexcept;
    /**
    * \return return maximum value and location to allow
    * efficient looping over multiple boxes.
    */
    AMREX_GPU_HOST_DEVICE
    void  maxIndex (const Box&     subbox,
                          Real&    max_value,
                          IntVect& max_index,
                          int      comp = 0) const noexcept;


    /**
    * \brief Compute mask array with value of 1 in cells where
    * BaseFab has value less than val, 0 otherwise.
    * mask is resized by this function.
    * The number of cells marked with 1 returned.
    */
    AMREX_GPU_HOST_DEVICE
    int maskLT (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const noexcept;
    //! Same as above except mark cells with value less than or equal to val.
    AMREX_GPU_HOST_DEVICE
    int maskLE (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const noexcept;

    //! Same as above except mark cells with value equal to val.
    AMREX_GPU_HOST_DEVICE
    int maskEQ (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const noexcept;
    //! Same as above except mark cells with value greater than val.
    AMREX_GPU_HOST_DEVICE
    int maskGT (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const noexcept;
    //! Same as above except mark cells with value greater than or equal to val.
    AMREX_GPU_HOST_DEVICE
    int maskGE (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const noexcept;
    //! Returns sum of given component of FAB state vector.
    AMREX_GPU_HOST_DEVICE
    T sum (int comp,
           int numcomp = 1) const noexcept;
    //! Compute sum of given component of FAB state vector in given subbox.
    AMREX_GPU_HOST_DEVICE
    T sum (const Box& subbox,
           int        comp,
           int        numcomp = 1) const noexcept;
    //! Most general version, specify subbox and which components.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& invert (T          v,
                                  const Box& subbox,
                                  int        comp=0,
                                  int        numcomp=1) noexcept;
    //! As above except on entire domain.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& invert (T   v,
                                  int comp,
                                  int numcomp=1) noexcept;

    //! Negate BaseFab, most general.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& negate (const Box& subbox,
                                  int        comp=0,
                                  int        numcomp=1) noexcept;
    //! As above, except on entire domain.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& negate (int comp,
                                  int numcomp=1) noexcept;

    //! Scalar addition (a[i] <- a[i] + r), most general.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& plus (T          r,
                                const Box& b,
                                int        comp=0,
                                int        numcomp=1) noexcept;

    //! As above, except on entire domain.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& plus (T   r,
                                int comp,
                                int numcomp=1) noexcept;

    /**
    * \brief Add src components (srccomp:srccomp+numcomp-1) to
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& plus (const BaseFab<T,Allocator>& src,
                                int               srccomp,
                                int               destcomp,
                                int               numcomp=1) noexcept;
    /**
    * \brief Same as above except addition is restricted to intersection
    * of subbox and src FAB. NOTE: subbox must be contained in this
    * FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& plus (const BaseFab<T,Allocator>& src,
                                const Box&        subbox,
                                int               srccomp,
                                int               destcomp,
                                int               numcomp=1) noexcept;
    /**
    * \brief Add srcbox region of src FAB to destbox region of this FAB.
    * The srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& plus (const BaseFab<T,Allocator>& src,
                                const Box&        srcbox,
                                const Box&        destbox,
                                int               srccomp,
                                int               destcomp,
                                int               numcomp=1) noexcept;

    //! Atomic FAB addition (a[i] <- a[i] + b[i]).
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& atomicAdd (const BaseFab<T,Allocator>& src) noexcept;

    /**
    * \brief Atomically add src components (srccomp:srccomp+numcomp-1) to
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& atomicAdd (const BaseFab<T,Allocator>& src,
                                     int               srccomp,
                                     int               destcomp,
                                     int               numcomp=1) noexcept;
    /**
    * \brief Same as above except addition is restricted to intersection
    * of subbox and src FAB. NOTE: subbox must be contained in this
    * FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& atomicAdd (const BaseFab<T,Allocator>& src,
                                     const Box&        subbox,
                                     int               srccomp,
                                     int               destcomp,
                                     int               numcomp=1) noexcept;
    /**
    * \brief Atomically add srcbox region of src FAB to destbox region of this FAB.
    * The srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& atomicAdd (const BaseFab<T,Allocator>& src,
                                     const Box&        srcbox,
                                     const Box&        destbox,
                                     int               srccomp,
                                     int               destcomp,
                                     int               numcomp=1) noexcept;

    //! FAB SAXPY (y[i] <- y[i] + a * x[i]), in place.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& saxpy (T a, const BaseFab<T,Allocator>& x,
                                 const Box&        srcbox,
                                 const Box&        destbox,
                                 int               srccomp,
                                 int               destcomp,
                                 int               numcomp=1) noexcept;
    //! FAB SAXPY (y[i] <- y[i] + a * x[i]), in place.  All components.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& saxpy (T a, const BaseFab<T,Allocator>& x) noexcept;

    //! FAB XPAY (y[i] <- x[i] + a * y[i])
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& xpay (T a, const BaseFab<T,Allocator>& x,
                                const Box&        srcbox,
                                const Box&        destbox,
                                int               srccomp,
                                int               destcomp,
                                int               numcomp=1) noexcept;

    //! y[i] <- y[i] + x1[i] * x2[i])
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& addproduct (const Box&        destbox,
                                      int               destcomp,
                                      int               numcomp,
                                      const BaseFab<T,Allocator>& src1,
                                      int               comp1,
                                      const BaseFab<T,Allocator>& src2,
                                      int               comp2) noexcept;

    /**
    * \brief Subtract src components (srccomp:srccomp+numcomp-1) to
    * this FABs components (destcomp:destcomp+numcomp-1) where
    * the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& minus (const BaseFab<T,Allocator>& src,
                                 int               srccomp,
                                 int               destcomp,
                                 int               numcomp=1) noexcept;
    /**
    * \brief Same as above except subtraction is restricted to intersection
    * of subbox and src FAB.  NOTE: subbox must be contained in
    * this FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& minus (const BaseFab<T,Allocator>& src,
                                 const Box&        subbox,
                                 int               srccomp,
                                 int               destcomp,
                                 int               numcomp=1) noexcept;
    /**
    * \brief Subtract srcbox region of src FAB from destbox region
    * of this FAB. srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& minus (const BaseFab<T,Allocator>& src,
                                 const Box&        srcbox,
                                 const Box&        destbox,
                                 int               srccomp,
                                 int               destcomp,
                                 int               numcomp=1) noexcept;

    //! Scalar multiplication, except control which components are multiplied.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& mult (T   r,
                                int comp,
                                int numcomp=1) noexcept;
    /**
    * \brief As above, except specify sub-box.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& mult (T          r,
                                const Box& b,
                                int        comp=0,
                                int        numcomp=1) noexcept;

    /**
    * \brief Multiply src components (srccomp:srccomp+numcomp-1) with
    * this FABs components (destcomp:destcomp+numcomp-1) where
    * the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& mult (const BaseFab<T,Allocator>& src,
                                int               srccomp,
                                int               destcomp,
                                int               numcomp=1) noexcept;

    /**
    * \brief Same as above except multiplication is restricted to
    * intersection of subbox and src FAB.  NOTE: subbox must be
    * contained in this FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& mult (const BaseFab<T,Allocator>& src,
                                const Box&        subbox,
                                int               srccomp,
                                int               destcomp,
                                int               numcomp=1) noexcept;

    /**
    * \brief Multiply srcbox region of src FAB with destbox region
    * of this FAB. The srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& mult (const BaseFab<T,Allocator>& src,
                                const Box&        srcbox,
                                const Box&        destbox,
                                int               srccomp,
                                int               destcomp,
                                int               numcomp=1) noexcept;

    //! As above except specify which components.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& divide (T   r,
                                  int comp,
                                  int numcomp=1) noexcept;

    //! As above except specify sub-box.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& divide (T          r,
                                  const Box& b,
                                  int        comp=0,
                                  int        numcomp=1) noexcept;

    /**
    * \brief This FAB is numerator, src FAB is denominator
    * divide src components (srccomp:srccomp+numcomp-1) into
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& divide (const BaseFab<T,Allocator>& src,
                                  int               srccomp,
                                  int               destcomp,
                                  int               numcomp=1) noexcept;
    /**
    * \brief Same as above except division is restricted to
    * intersection of subbox and src FAB.  NOTE: subbox must be
    * contained in this FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& divide (const BaseFab<T,Allocator>& src,
                                  const Box&        subbox,
                                  int               srccomp,
                                  int               destcomp,
                                  int               numcomp=1) noexcept;
    /**
    * \brief destbox region of this FAB is numerator. srcbox regions of
    * src FAB is denominator. srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& divide (const BaseFab<T,Allocator>& src,
                                  const Box&        srcbox,
                                  const Box&        destbox,
                                  int               srccomp,
                                  int               destcomp,
                                  int               numcomp=1) noexcept;
    /**
    * \brief Divide wherever "src" is "true" or "non-zero".
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& protected_divide (const BaseFab<T,Allocator>& src) noexcept;

    /**
    * \brief Divide wherever "src" is "true" or "non-zero".
    * This FAB is numerator, src FAB is denominator
    * divide src components (srccomp:srccomp+numcomp-1) into
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& protected_divide (const BaseFab<T,Allocator>& src,
                                            int               srccomp,
                                            int               destcomp,
                                            int               numcomp=1) noexcept;

    /**
    * \brief Divide wherever "src" is "true" or "non-zero".
    * Same as above except division is restricted to
    * intersection of subbox and src FAB.  NOTE: subbox must be
    * contained in this FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& protected_divide (const BaseFab<T,Allocator>& src,
                                            const Box&        subbox,
                                            int               srccomp,
                                            int               destcomp,
                                            int               numcomp=1) noexcept;

    /**
    * Divide wherever "src" is "true" or "non-zero".
    * destbox region of this FAB is numerator. srcbox regions of
    * src FAB is denominator. srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& protected_divide (const BaseFab<T,Allocator>& src,
                                            const Box&        srcbox,
                                            const Box&        destbox,
                                            int               srccomp,
                                            int               destcomp,
                                            int               numcomp=1) noexcept;

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
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& linInterp (const BaseFab<T,Allocator>& f1,
                                     const Box&        b1,
                                     int               comp1,
                                     const BaseFab<T,Allocator>& f2,
                                     const Box&        b2,
                                     int               comp2,
                                     Real              t1,
                                     Real              t2,
                                     Real              t,
                                     const Box&        b,
                                     int               comp,
                                     int               numcomp = 1) noexcept;

    AMREX_GPU_HOST_DEVICE
    //! Version of linInterp() in which b, b1, & b2 are the same.
    BaseFab<T,Allocator>& linInterp (const BaseFab<T,Allocator>& f1,
                                     int               comp1,
                                     const BaseFab<T,Allocator>& f2,
                                     int               comp2,
                                     Real              t1,
                                     Real              t2,
                                     Real              t,
                                     const Box&        b,
                                     int               comp,
                                     int               numcomp = 1) noexcept;

    /**
    * \brief Linear combination.  Result is alpha*f1 + beta*f2.
    * Data is taken from b1 region of f1, b2 region of f2
    * and stored in b region of this FAB.
    * Boxes b, b1 and b2 must be the same size.
    * Data is taken from component comp1 of f1, comp2 of f2,
    * and stored in component comp of this FAB.
    * This FAB is returned as a reference for chaining.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& linComb (const BaseFab<T,Allocator>& f1,
                                   const Box&        b1,
                                   int               comp1,
                                   const BaseFab<T,Allocator>& f2,
                                   const Box&        b2,
                                   int               comp2,
                                   Real              alpha,
                                   Real              beta,
                                   const Box&        b,
                                   int               comp,
                                   int               numcomp = 1) noexcept;

    //! Dot product of x (i.e.,this) and y
    AMREX_GPU_HOST_DEVICE
    T dot (const Box& xbx, int xcomp,
	   const BaseFab<T,Allocator>& y, const Box& ybx, int ycomp,
	   int numcomp = 1) const noexcept;

    AMREX_GPU_HOST_DEVICE
    T dotmask (const BaseFab<int>& mask, const Box& xbx, int xcomp,
               const BaseFab<T,Allocator>& y, const Box& ybx, int ycomp,
               int numcomp) const noexcept;

    //! Change the Box type without change the length
    AMREX_GPU_HOST_DEVICE
    void SetBoxType (const IndexType& typ) noexcept { this->domain.setType(typ); }

    template <class F> //!< AMREX_GPU_HOST_DEVICE
    void ForEach (const Box& b, int c, int nc, F f) noexcept {
        ForEachImpl(*this, b, c, nc, f);
    }
    template <class F> //!< AMREX_GPU_HOST_DEVICE
    void ForEach (const Box& b, int c, int nc, F f) const noexcept {
        ForEachImpl(*this, b, c, nc, f);
    }
    template <class F> //!< AMREX_GPU_HOST_DEVICE
    void ForEachIV (const Box& b, int c, int nc, F f) noexcept {
        ForEachIVImpl(*this, b, c, nc, f);
    }
    template <class F> //!< AMREX_GPU_HOST_DEVICE
    void ForEachIV (const Box& b, int c, int nc, F f) const noexcept {
        ForEachIVImpl(*this, b, c, nc, f);
    }
    template <class F> //!< AMREX_GPU_HOST_DEVICE
    void ForEach (const Box& dstbox, int dstcomp, int numcomp,
                  const BaseFab<T,Allocator>& src, int srccomp, F f) noexcept;

    template <class F> //!< AMREX_GPU_HOST_DEVICE
    void ForEach (const Box& dstbox, int dstcomp, int numcomp,
                  const BaseFab<T,Allocator>& src, const Box& srcbox, int srccomp, F f) noexcept;

    template <typename P, class F> //!< AMREX_GPU_HOST_DEVICE
    P Accumulate (const Box& b, int c, int nc, P init, F f) const noexcept;

    template <class F> //!< AMREX_GPU_HOST_DEVICE
    void Transform (T* dst, const Box& b, int c, int nc, F f) const noexcept;

    template <class F> //!< AMREX_GPU_HOST_DEVICE
    void Transform (const Box& b, int c, int nc, T const* src, F f) noexcept;

    //
    // New interfaces
    //

    //! Set value on the whole domain and all components
    AMREX_GPU_HOST_DEVICE
    void setVal (T val) noexcept;
    //
    //! Do nothing if bx is empty.
    AMREX_GPU_HOST_DEVICE
    void setVal (T val, Box bx, DestComp dcomp, NumComps ncomp) noexcept;

    AMREX_GPU_HOST_DEVICE
    void setValIf (T val, const BaseFab<int>& mask) noexcept;
    //
    //! Do nothing if bx is empty.
    AMREX_GPU_HOST_DEVICE
    void setValIf (T val, Box bx, const BaseFab<int>& mask, DestComp dcomp, NumComps ncomp) noexcept;

    AMREX_GPU_HOST_DEVICE
    void setValIfNot (T val, const BaseFab<int>& mask) noexcept;
    //
    //! Do nothing if bx is empty.
    AMREX_GPU_HOST_DEVICE
    void setValIfNot (T val, Box bx, const BaseFab<int>& mask, DestComp dcomp, NumComps ncomp) noexcept;

    //! setVal on the complement of bx in the fab's domain
    void setComplement (T val, Box const& bx, DestComp dcomp, NumComps ncomp) noexcept;

    /**
    * copy is performed on the intersection of dest and src fabs.
    * All components of dest fab are copied. src fab must have enough
    * components (more is OK).
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& copy (const BaseFab<T,Allocator>& src) noexcept;
    //
    //! Do nothing if bx does not intersect with src fab.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& copy (const BaseFab<T,Allocator>& src, Box bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept;

    //! Scalar addition on the whole domain and all components
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& plus (T val) noexcept;
    //
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& operator+= (T val) noexcept;
    //
    //! Do nothing if bx is empty.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& plus (T val, Box bx, DestComp dcomp, NumComps ncomp) noexcept;
    /**
    * Fab addition is performed on the intersection of dest and src fabs.
    * All components of dest fab are copied. src fab must have enough
    * components (more is OK).
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& plus (const BaseFab<T,Allocator>& src) noexcept;
    //
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& operator+= (const BaseFab<T,Allocator>& src) noexcept;
    //
    //! Do nothing if bx does not intersect with src fab.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& plus (const BaseFab<T,Allocator>& src, Box bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept;

    //! Scalar subtraction on the whole domain and all components
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& minus (T val) noexcept;
    //
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& operator-= (T val) noexcept;
    //
    //! Do nothing if bx is empty.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& minus (T val, Box bx, DestComp dcomp, NumComps ncomp) noexcept;
    /**
    * Fab subtraction is performed on the intersection of dest and src fabs.
    * All components of dest fab are copied. src fab must have enough
    * components (more is OK).
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& minus (const BaseFab<T,Allocator>& src) noexcept;
    //
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& operator-= (const BaseFab<T,Allocator>& src) noexcept;
    //
    //! Do nothing if bx does not intersect with src fab.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& minus (const BaseFab<T,Allocator>& src, Box bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept;

    //! Scalar multiplication on the whole domain and all components
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& mult (T val) noexcept;
    //
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& operator*= (T val) noexcept;
    //
    //! Do nothing if bx is empty.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& mult (T val, Box bx, DestComp dcomp, NumComps ncomp) noexcept;
    /**
    * Fab multiplication is performed on the intersection of dest and src fabs.
    * All components of dest fab are copied. src fab must have enough
    * components (more is OK).
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& mult (const BaseFab<T,Allocator>& src) noexcept;
    //
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& operator*= (const BaseFab<T,Allocator>& src) noexcept;
    //
    //! Do nothing if bx does not intersect with src fab.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& mult (const BaseFab<T,Allocator>& src, Box bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept;

    //! Scalar division on the whole domain and all components
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& divide (T val) noexcept;
    //
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& operator/= (T val) noexcept;
    //
    //! Do nothing if bx is empty.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& divide (T val, Box bx, DestComp dcomp, NumComps ncomp) noexcept;
    /**
    * Fab division is performed on the intersection of dest and src fabs.
    * All components of dest fab are copied. src fab must have enough
    * components (more is OK).
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& divide (const BaseFab<T,Allocator>& src) noexcept;
    //
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& operator/= (const BaseFab<T,Allocator>& src) noexcept;
    //
    //! Do nothing if bx does not intersect with src fab.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& divide (const BaseFab<T,Allocator>& src, Box bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept;

    //! on the whole domain and all components
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& negate () noexcept;
    //
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& negate (const Box& bx, DestComp dcomp, NumComps ncomp) noexcept;

    //! Fab <- Fab/r on the whole domain and all components
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& invert (T r) noexcept;
    //
    AMREX_GPU_HOST_DEVICE
    BaseFab<T,Allocator>& invert (T r, const Box& bx, DestComp dcomp, NumComps ncomp) noexcept;

    //! Sum
    AMREX_GPU_HOST_DEVICE
    T sum (const Box& bx, DestComp dcomp, NumComps ncomp) const noexcept;

    //! Dot product of two Fabs
    AMREX_GPU_HOST_DEVICE
    T dot (const BaseFab<T,Allocator>& src, const Box& bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) const noexcept;

    //! Dot product of two Fabs with mask
    AMREX_GPU_HOST_DEVICE
    T dotmask (const BaseFab<T,Allocator>& src, const Box& bx, const BaseFab<int>& mask,
               SrcComp scomp, DestComp dcomp, NumComps ncomp) const noexcept;

protected:
    //! Allocates memory for the BaseFab<T,Allocator>.
    void define ();

#ifdef USE_PERILLA
public:
    LocalConnection l_con;
    RemoteConnection r_con;
    bool fireable;
    int padding[1024];
#endif

};


template <class B, class F>
//AMREX_GPU_HOST_DEVICE
inline
void
ForEachImpl (B& fab, const Box& b, int c, int nc, F f) noexcept
{
    AMREX_ASSERT(fab.contains(b));
    AMREX_ASSERT(c >= 0 && c + nc <= fab.nComp());

    const auto len3 = b.length3d();
    const int* blo = b.loVect();
    for (int n = c; n < c+nc; ++n) {
        for     (int k = 0; k < len3[2]; ++k) {
            for (int j = 0; j < len3[1]; ++j) {
                const IntVect line_begin{AMREX_D_DECL(blo[0],
                                                      blo[1]+j,
                                                      blo[2]+k)};
                auto d = fab.dataPtr(line_begin, n);
                for (int i = 0; i < len3[0]; ++i) {
                    f(*(d+i));
                }
            }
        }
    }
}

template <class B, class F>
inline
//AMREX_GPU_HOST_DEVICE
void
ForEachIVImpl (B& fab, const Box& b, int c, int nc, F f) noexcept
{
    AMREX_ASSERT(fab.contains(b));
    AMREX_ASSERT(c >= 0 && c + nc <= fab.nComp());

    const auto len3 = b.length3d();
    const int* blo = b.loVect();
    for (int n = c; n < c+nc; ++n) {
        for     (int k = 0; k < len3[2]; ++k) {
            for (int j = 0; j < len3[1]; ++j) {
                const IntVect line_begin{AMREX_D_DECL(blo[0],
                                                      blo[1]+j,
                                                      blo[2]+k)};
                auto d = fab.dataPtr(line_begin, n);
                for (int i = 0; i < len3[0]; ++i) {
                    f(*(d+i), IntVect{AMREX_D_DECL(blo[0]+i,
                                                   blo[1]+j,
                                                   blo[2]+k)});
                }
            }
        }
    }
}

template <class T, class Allocator>
template <class F>
void
BaseFab<T,Allocator>::ForEach (const Box& dstbox, int dstcomp, int numcomp,
                               const BaseFab<T,Allocator>& src, const Box& srcbox, int srccomp, F f) noexcept
{
    AMREX_ASSERT(contains(dstbox));
    AMREX_ASSERT(src.contains(srcbox));
    AMREX_ASSERT(dstbox.size() == srcbox.size());
    AMREX_ASSERT(dstcomp >= 0 && dstcomp+numcomp <= nComp());
    AMREX_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());

    const auto len3 = dstbox.length3d();
    const int* dlo = dstbox.loVect();
    const int* slo = srcbox.loVect();
    for (int n = 0; n < numcomp; ++n) {
        for     (int k = 0; k < len3[2]; ++k) {
            for (int j = 0; j < len3[1]; ++j) {
                const IntVect dline_begin{AMREX_D_DECL(dlo[0],
                                                       dlo[1]+j,
                                                       dlo[2]+k)};
                auto d = dataPtr(dline_begin, n+dstcomp);
                const IntVect sline_begin{AMREX_D_DECL(slo[0],
                                                       slo[1]+j,
                                                       slo[2]+k)};
                auto s = src.dataPtr(sline_begin, n+srccomp);
                for (int i = 0; i < len3[0]; ++i) {
                    f(*(d+i), *(s+i));
                }
            }
        }
    }
}

template <class T, class Allocator>
template <class F>
void
BaseFab<T,Allocator>::ForEach (const Box& dstbox, int dstcomp, int numcomp,
                               const BaseFab<T,Allocator>& src, int srccomp, F f) noexcept
{
    AMREX_ASSERT(contains(dstbox));
    AMREX_ASSERT(src.contains(dstbox));
    AMREX_ASSERT(dstcomp >= 0 && dstcomp+numcomp <= nComp());
    AMREX_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());

    const auto len3 = dstbox.length3d();
    const int* dlo = dstbox.loVect();
    for (int n = 0; n < numcomp; ++n) {
        for     (int k = 0; k < len3[2]; ++k) {
            for (int j = 0; j < len3[1]; ++j) {
                const IntVect dline_begin{AMREX_D_DECL(dlo[0],
                                                       dlo[1]+j,
                                                       dlo[2]+k)};
                auto d = dataPtr(dline_begin, n+dstcomp);
                auto s = src.dataPtr(dline_begin, n+srccomp);
                for (int i = 0; i < len3[0]; ++i) {
                    f(*(d+i), *(s+i));
                }
            }
        }
    }
}

template <class T, class Allocator>
template <typename P, class F>
P
BaseFab<T,Allocator>::Accumulate (const Box& b, int c, int nc, P init, F f) const noexcept
{
    const auto len3 = b.length3d();
    const int* blo = b.loVect();
    for (int n = c; n < c+nc; ++n) {
        for     (int k = 0; k < len3[2]; ++k) {
            for (int j = 0; j < len3[1]; ++j) {
                const IntVect line_begin{AMREX_D_DECL(blo[0],
                                                      blo[1]+j,
                                                      blo[2]+k)};
                const T* d = dataPtr(line_begin, n);
                for (int i = 0; i < len3[0]; ++i) {
                    init = f(init, *(d+i));
                }
            }
        }
    }
    return init;
}

template <class T, class Allocator>
template <class F>
void
BaseFab<T,Allocator>::Transform (T* dst, const Box& b, int c, int nc, F f) const noexcept
{
    AMREX_ASSERT(contains(b));
    AMREX_ASSERT(c >= 0 && c + nc <= nComp());

    const auto len3 = b.length3d();
    const int* blo = b.loVect();
    for (int n = c; n < c+nc; ++n) {
        for     (int k = 0; k < len3[2]; ++k) {
            for (int j = 0; j < len3[1]; ++j) {
                const IntVect line_begin{AMREX_D_DECL(blo[0],
                                                      blo[1]+j,
                                                      blo[2]+k)};
                auto s = dataPtr(line_begin, n);
                for (int i = 0; i < len3[0]; ++i) {
                    f(*(dst++), *(s+i));
                }
            }
        }
    }
}

template <class T, class Allocator>
template <class F>
void
BaseFab<T,Allocator>::Transform (const Box& b, int c, int nc, T const* src, F f) noexcept
{
    AMREX_ASSERT(contains(b));
    AMREX_ASSERT(c >= 0 && c + nc <= nComp());

    const auto len3 = b.length3d();
    const int* blo = b.loVect();
    for (int n = c; n < c+nc; ++n) {
        for     (int k = 0; k < len3[2]; ++k) {
            for (int j = 0; j < len3[1]; ++j) {
                const IntVect line_begin{AMREX_D_DECL(blo[0],
                                                      blo[1]+j,
                                                      blo[2]+k)};
                auto d = dataPtr(line_begin, n);
                for (int i = 0; i < len3[0]; ++i) {
                    f(*(d+i), *(src++));
                }
            }
        }
    }
}

template <class T, class Allocator>
inline
T*
BaseFab<T,Allocator>::dataPtr (const IntVect& p, int n) noexcept
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr + (this->domain.index(p)+n*this->domain.numPts());
}

template <class T, class Allocator>
inline
const T*
BaseFab<T,Allocator>::dataPtr (const IntVect& p, int n) const noexcept
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr + (this->domain.index(p)+n*this->domain.numPts());
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::prefetchToHost () const noexcept
{
#ifdef AMREX_USE_CUDA
    std::size_t s = sizeof(T)*this->nvar*this->domain.numPts();
    AMREX_GPU_SAFE_CALL(cudaMemPrefetchAsync(this->dptr, s,
                                             cudaCpuDeviceId,
                                             Cuda::Device::cudaStream()));
#endif
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::prefetchToDevice () const noexcept
{
#ifdef AMREX_USE_CUDA
    std::size_t s = sizeof(T)*this->nvar*this->domain.numPts();
    AMREX_GPU_SAFE_CALL(cudaMemPrefetchAsync(this->dptr, s,
                                             Cuda::Device::deviceId(),
                                             Cuda::Device::cudaStream()));
#endif
}

template <class T, class Allocator>
inline
T&
BaseFab<T,Allocator>::operator() (const IntVect& p, int n) noexcept
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr[this->domain.index(p)+n*this->domain.numPts()];
}

template <class T, class Allocator>
inline
T&
BaseFab<T,Allocator>::operator() (const IntVect& p) noexcept
{
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr[this->domain.index(p)];
}

template <class T, class Allocator>
inline
const T&
BaseFab<T,Allocator>::operator() (const IntVect& p, int n) const noexcept
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr[this->domain.index(p)+n*this->domain.numPts()];
}

template <class T, class Allocator>
inline
const T&
BaseFab<T,Allocator>::operator() (const IntVect& p) const noexcept
{
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr[this->domain.index(p)];
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::getVal  (T*             data,
                               const IntVect& pos,
                               int            n,
                               int            numcomp) const noexcept
{
    const int loc      = this->domain.index(pos);
    const long sz    = this->domain.numPts();

    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(n >= 0 && n + numcomp <= this->nvar);

    for (int k = 0; k < numcomp; k++)
        data[k] = this->dptr[loc+(n+k)*sz];
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::getVal (T*             data,
                              const IntVect& pos) const noexcept
{
    getVal(data,pos,0,this->nvar);
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::shift (const IntVect& v) noexcept
{
    this->domain += v;
    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::shift (int idir, int n_cell) noexcept
{
    this->domain.shift(idir,n_cell);
    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator> &
BaseFab<T,Allocator>::shiftHalf (const IntVect& v) noexcept
{
    this->domain.shiftHalf(v);
    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator> &
BaseFab<T,Allocator>::shiftHalf (int idir, int n_cell) noexcept
{
    this->domain.shiftHalf(idir,n_cell);
    return *this;
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setVal (T x, const Box& bx, int n) noexcept
{
    this->setVal(x, bx, DestComp{n}, NumComps{1});
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setVal (T x, int n) noexcept
{
    this->setVal(x, this->domain, DestComp{n}, NumComps{1});
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setVal (T x, const Box& bx, int dcomp, int ncomp) noexcept
{
    this->setVal(x, bx, DestComp{dcomp}, NumComps{ncomp});
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setValIfNot (T val, const Box& bx, const BaseFab<int>& mask, int ns, int num) noexcept
{
    this->setValIfNot(val, bx, mask, DestComp{ns}, NumComps{num});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::copy (const BaseFab<T,Allocator>& src,
                            const Box&        srcbox,
                            int               srccomp,
                            const Box&        destbox,
                            int               destcomp,
                            int               numcomp) noexcept
{
    AMREX_ASSERT(destbox.ok());
    AMREX_ASSERT(srcbox.sameSize(destbox));
    AMREX_ASSERT(src.box().contains(srcbox));
    AMREX_ASSERT(this->domain.contains(destbox));
    AMREX_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    AMREX_ASSERT(destcomp >= 0 && destcomp+numcomp <= this->nvar);

    const auto len = amrex::length(destbox);
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const auto dp  =     view(dlo, destcomp);
    const auto sp  = src.view(slo, srccomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) = sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::copy (const BaseFab<T,Allocator>& src, const Box& destbox) noexcept
{
    return this->copy(src, destbox, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::copy (const BaseFab<T,Allocator>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    return copy(src, this->domain, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::define ()
{
    AMREX_ASSERT(this->nvar > 0);
    AMREX_ASSERT(this->dptr == 0);
    AMREX_ASSERT(this->domain.numPts() > 0);
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
    AMREX_ASSERT(std::numeric_limits<long>::max()/this->nvar > this->domain.numPts());
#else
    AMREX_ASSERT(LONG_MAX/this->nvar > this->domain.numPts());
#endif

    this->truesize  = this->nvar*this->domain.numPts();
    this->ptr_owner = true;
    this->dptr = static_cast<T*>(Allocator{}.alloc(this->truesize*sizeof(T)));
    //
    // Now call T::T() on the raw memory so we have valid Ts.
    //
    T* ptr = this->dptr;
    //
    // Note this must be long not int for very large (e.g.,1024^3) boxes.
    //
    for (long i = 0; i < this->truesize; i++, ptr++)
    {
        new (ptr) T;
    }

    amrex::update_fab_stats(this->domain.numPts(), this->truesize, sizeof(T));
}

template <class T, class Allocator>
BaseFab<T,Allocator>::BaseFab ()
    : BaseFabData<T>{nullptr, Box(), 0, 0L, false, false}
{}

template <class T, class Allocator>
BaseFab<T,Allocator>::BaseFab (const Box& bx,
                               int        n,
                               bool       alloc,
                               bool       shared)
    : BaseFabData<T>{nullptr, bx, n, 0L, false, shared}
{
    if (!this->shared_memory && alloc) define();
}

template <class T, class Allocator>
BaseFab<T,Allocator>::BaseFab (const BaseFab<T,Allocator>& rhs, MakeType make_type, int scomp, int ncomp)
    : BaseFabData<T>{const_cast<T*>(rhs.dataPtr(scomp)), rhs.domain, ncomp,
                     ncomp*rhs.domain.numPts(), false, false}
{
    AMREX_ASSERT(scomp+ncomp <= rhs.nComp());
    if (make_type == amrex::make_deep_copy)
    {
        this->dptr = nullptr;
        define();
        this->copy(rhs, this->domain, scomp, this->domain, 0, ncomp);
    } else if (make_type == amrex::make_alias) {
        ; // nothing to do
    } else {
        amrex::Abort("BaseFab: unknown MakeType");
    }
}

template<class T, class Allocator>
BaseFab<T,Allocator>::BaseFab (const Box& bx, int ncomp, T* p)
    : BaseFabData<T>{p, bx, ncomp, bx.numPts()*ncomp, false, false}
{
}

template <class T, class Allocator>
BaseFab<T,Allocator>::~BaseFab ()
{
    clear();
}

template <class T, class Allocator>
BaseFab<T,Allocator>::BaseFab (BaseFab<T,Allocator>&& rhs) noexcept
    : BaseFabData<T>{rhs.dptr, rhs.domain, rhs.nvar,rhs.truesize,
                     rhs.ptr_owner, rhs.shared_memory}
{
    rhs.dptr = nullptr;
    rhs.ptr_owner = false;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::operator= (T t) noexcept
{
    setVal(t);
    return *this;
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::resize (const Box& b, int n)
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

template <class T, class Allocator>
template <class,class>
Elixir
BaseFab<T,Allocator>::elixir () noexcept
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
    return Elixir((o ? this->dptr : nullptr), Allocator{}.arena());
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::clear () noexcept
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

	    T* ptr = this->dptr;

	    for (long i = 0; i < this->truesize; i++, ptr++)
	    {
		ptr->~T();
	    }

	    Allocator{}.free(this->dptr);

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

template <class T, class Allocator>
std::size_t
BaseFab<T,Allocator>::copyToMem (const Box& srcbox,
                                 int        srccomp,
                                 int        numcomp,
                                 void*      dst) const noexcept
{
    BL_ASSERT(box().contains(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= nComp());

    if (srcbox.ok())
    {
        T* AMREX_RESTRICT dp = static_cast<T*>(dst);

        const auto len = amrex::length(srcbox);
        const auto lo  = amrex::lbound(srcbox);
        const auto sp  = view(lo, srccomp);

        std::size_t offset = 0;
        for (int n = 0; n < numcomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = 0; i < len.x; ++i) {
                        dp[offset+i] = sp(i,j,k,n);
                    }
                    offset += len.x;
                }
            }
        }
        return sizeof(T)*offset;
    }
    else
    {
        return 0;
    }
}

template <class T, class Allocator>
std::size_t
BaseFab<T,Allocator>::copyFromMem (const Box&  dstbox,
                                   int         dstcomp,
                                   int         numcomp,
                                   const void* src) noexcept
{
    BL_ASSERT(box().contains(dstbox));
    BL_ASSERT(dstcomp >= 0 && dstcomp+numcomp <= nComp());

    if (dstbox.ok())
    {
        T const* AMREX_RESTRICT sp = static_cast<T const*>(src);

        const auto len = amrex::length(dstbox);
        const auto lo  = amrex::lbound(dstbox);
        const auto dp  = view(lo, dstcomp);

        std::size_t offset = 0;
        for (int n = 0; n < numcomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = 0; i < len.x; ++i) {
                        dp(i,j,k,n) = sp[offset+i];
                    }
                    offset += len.x;
                }
            }
        }
        return sizeof(T)*offset;
    }
    else
    {
        return 0;
    }
}

template <class T, class Allocator>
std::size_t
BaseFab<T,Allocator>::addFromMem (const Box&  dstbox,
                                  int         dstcomp,
                                  int         numcomp,
                                  const void* src) noexcept
{
    BL_ASSERT(box().contains(dstbox));
    BL_ASSERT(dstcomp >= 0 && dstcomp+numcomp <= nComp());

    if (dstbox.ok())
    {
        T const* AMREX_RESTRICT sp = static_cast<T const*>(src);

        const auto len = amrex::length(dstbox);
        const auto lo  = amrex::lbound(dstbox);
        const auto dp  = view(lo, dstcomp);

        std::size_t offset = 0;
        for (int n = 0; n < numcomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = 0; i < len.x; ++i) {
                        dp(i,j,k,n) += sp[offset+i];
                    }
                    offset += len.x;
                }
            }
        }
        return sizeof(T)*offset;
    }
    else
    {
        return 0;
    }
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setComplement (T x, const Box& b, int ns, int num) noexcept
{
    setComplement(x, b, DestComp{ns}, NumComps{num});
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::abs () noexcept
{
    this->abs(this->domain,0,this->nvar);
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::abs (int comp, int numcomp) noexcept
{
    this->abs(this->domain,comp,numcomp);
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::abs (const Box& subbox,
                           int        comp,
                           int        numcomp) noexcept
{
    const auto len = amrex::length(subbox);
    const auto lo  = amrex::lbound(subbox);
    const auto dp  = view(lo, comp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) = std::abs(dp(i,j,k,n));
                }
            }
        }
    }
}

template <class T, class Allocator>
Real
BaseFab<T,Allocator>::norminfmask (const Box& subbox, const BaseFab<int>& mask,
                                   int scomp, int ncomp) const noexcept
{
    BL_ASSERT(this->domain.contains(subbox));
    BL_ASSERT(scomp >= 0 && scomp + ncomp <= this->nvar);

    Real r = 0.0;

    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   =      view(lo, scomp);
    const auto mp   = mask.view(lo);

    for (int n = 0; n < ncomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    if (mp(i,j,k,0)) {
                        Real t = static_cast<Real>(std::abs(dp(i,j,k,n)));
                        r = (r > t) ? r : t;
                    }
                }
            }
        }
    }

    return r;
}

template <class T, class Allocator>
Real
BaseFab<T,Allocator>::norm (int p,
                            int comp,
                            int numcomp) const noexcept
{
    return norm(this->domain,p,comp,numcomp);
}

template <class T, class Allocator>
Real
BaseFab<T,Allocator>::norm (const Box& subbox,
                            int        p,
                            int        comp,
                            int        numcomp) const noexcept
{
    BL_ASSERT(this->domain.contains(subbox));
    BL_ASSERT(comp >= 0 && comp + numcomp <= this->nvar);

    Real nrm = 0.;

    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

    if (p == 0)
    {
        for (int n = 0; n < numcomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        Real t = static_cast<Real>(std::abs(dp(i,j,k,n)));
                        nrm = (nrm > t) ? nrm : t;
                    }
                }
            }
        }
    }
    else if (p == 1)
    {
        for (int n = 0; n < numcomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        nrm += std::abs(dp(i,j,k,n));
                    }
                }
            }
        }
    }
    else if (p == 2)
    {
        for (int n = 0; n < numcomp; ++n) {
            for         (int k = 0; k < len.z; ++k) {
                for     (int j = 0; j < len.y; ++j) {
                    for (int i = 0; i < len.x; ++i) {
                        nrm += dp(i,j,k,n)*dp(i,j,k,n);
                    }
                }
            }
        }
        nrm = std::sqrt(nrm);
    }
    else
    {
        amrex::Error("BaseFab<T,Allocator>::norm: wrong p");
    }

    return nrm;
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::min (int comp) const noexcept
{
    return this->min(this->domain,comp);
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::min (const Box& subbox, int comp) const noexcept
{
    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
    T r = std::numeric_limits<T>::max();
#else
    T r = *dp.p;
#endif
    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                T t = dp(i,j,k,0);
                r = (r < t) ? r : t;
            }
        }
    }

    return r;
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::max (int comp) const noexcept
{
    return this->max(this->domain,comp);
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::max (const Box& subbox, int comp) const noexcept
{
    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
    T r = std::numeric_limits<T>::lowest();
#else
    T r = *dp.p;
#endif

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                T t = dp(i,j,k,0);
                r = (r > t) ? r : t;
            }
        }
    }

    return r;
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::maxabs (int comp) const noexcept
{
    return this->maxabs(this->domain,comp);
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::maxabs (const Box& subbox, int comp) const noexcept
{
    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
    T r = std::numeric_limits<T>::lowest();
#else
    T r = *dp.p;
#endif

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                T t = std::abs(dp(i,j,k,0));
                r = (r > t) ? r : t;
            }
        }
    }

    return r;
}


template <class T, class Allocator>
IntVect
BaseFab<T,Allocator>::minIndex (int comp) const noexcept
{
    return this->minIndex(this->domain,comp);
}

template <class T, class Allocator>
IntVect
BaseFab<T,Allocator>::minIndex (const Box& subbox, int comp) const noexcept
{
    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
    AMREX_D_DECL(int imin, jmin, kmin);
    T min_val = std::numeric_limits<T>::max();
#else
    AMREX_D_DECL(int imin=0, jmin=0, kmin=0);
    T min_val = *dp.p;
#endif

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                if (dp(i,j,k,0) < min_val) {
                    min_val = dp(i,j,k,0);
                    AMREX_D_TERM(imin = i;,  jmin = j;,  kmin = k;)
                }
            }
        }
    }

    return IntVect(AMREX_D_DECL(imin+lo.x,jmin+lo.y,kmin+lo.z));
}

template <class T, class Allocator>
IntVect
BaseFab<T,Allocator>::indexFromValue (Real value, const Box& subbox, int comp) const noexcept
{
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
    AMREX_D_DECL(int imin = std::numeric_limits<int>::lowest(),
                     jmin = std::numeric_limits<int>::lowest(),
                     kmin = std::numeric_limits<int>::lowest());
#else
    AMREX_D_DECL(int imin = INT_MIN, jmin = INT_MIN, kmin = INT_MIN);
#endif

    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                if (dp(i,j,k,0) == value) {
                    AMREX_D_TERM(imin = i;,  jmin = j;,  kmin = k;)
                }
            }
        }
    }

    return IntVect(AMREX_D_DECL(imin+lo.x,jmin+lo.y,kmin+lo.z));
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::minIndex (const Box& subbox, Real& min_val, IntVect& min_idx, int comp) const noexcept
{
    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                if (dp(i,j,k,0) < min_val) {
                    min_val = dp(i,j,k,0);
                    min_idx = IntVect(AMREX_D_DECL(i,j,k));
                }
            }
        }
    }
    AMREX_D_TERM(min_idx[0] += lo.x;, min_idx[1] += lo.y;, min_idx[2] += lo.z);
}

template <class T, class Allocator>
IntVect
BaseFab<T,Allocator>::maxIndex (int comp) const noexcept
{
    return this->maxIndex(this->domain,comp);
}

template <class T, class Allocator>
IntVect
BaseFab<T,Allocator>::maxIndex (const Box& subbox,
                                int        comp) const noexcept
{
    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
    AMREX_D_DECL(int imax, jmax, kmax);
    T max_val = std::numeric_limits<T>::lowest();
#else
    AMREX_D_DECL(int imax=0, jmax=0, kmax=0);
    T max_val = *dp.p;
#endif

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                if (dp(i,j,k,0) > max_val) {
                    max_val = dp(i,j,k,0);
                    AMREX_D_TERM(imax = i;, jmax = j;, kmax = k;)
                }
            }
        }
    }

    return IntVect(AMREX_D_DECL(imax+lo.x,jmax+lo.y,kmax+lo.z));
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::maxIndex (const Box& subbox, Real& max_val, IntVect& max_idx, int comp) const noexcept
{
    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                if (dp(i,j,k,0) > max_val) {
                    max_val = dp(i,j,k,0);
                    max_idx = IntVect(AMREX_D_DECL(i,j,k));
                }
            }
        }
    }
    AMREX_D_TERM(max_idx[0] += lo.x;, max_idx[1] += lo.y;, max_idx[2] += lo.z);
}

template <class T, class Allocator>
int
BaseFab<T,Allocator>::maskLT (BaseFab<int>& mask,
                              T             val,
                              int           comp) const noexcept
{
    mask.resize(this->domain,1);
    mask.setVal(0);

    int cnt = 0;

    const auto len  = amrex::length(this->domain);
    const auto lo   = amrex::lbound(this->domain);
    const auto dp   =      view(lo, comp);
    const auto mp   = mask.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                if (dp(i,j,k,0) < val) {
                    mp(i,j,k,0) = 1;
                    ++cnt;
                }
            }
        }
    }

    return cnt;
}

template <class T, class Allocator>
int
BaseFab<T,Allocator>::maskLE (BaseFab<int>& mask,
                              T             val,
                              int           comp) const noexcept
{
    mask.resize(this->domain,1);
    mask.setVal(0);

    int cnt = 0;

    const auto len  = amrex::length(this->domain);
    const auto lo   = amrex::lbound(this->domain);
    const auto dp   =      view(lo, comp);
    const auto mp   = mask.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                if (dp(i,j,k,0) <= val) {
                    mp(i,j,k,0) = 1;
                    ++cnt;
                }
            }
        }
    }

    return cnt;
}

template <class T, class Allocator>
int
BaseFab<T,Allocator>::maskEQ (BaseFab<int>& mask,
                              T             val,
                              int           comp) const noexcept
{
    mask.resize(this->domain,1);
    mask.setVal(0);

    int cnt = 0;

    const auto len  = amrex::length(this->domain);
    const auto lo   = amrex::lbound(this->domain);
    const auto dp   =      view(lo, comp);
    const auto mp   = mask.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                if (dp(i,j,k,0) == val) {
                    mp(i,j,k,0) = 1;
                    ++cnt;
                }
            }
        }
    }

    return cnt;
}

template <class T, class Allocator>
int
BaseFab<T,Allocator>::maskGT (BaseFab<int>& mask,
                              T             val,
                              int           comp) const noexcept
{
    mask.resize(this->domain,1);
    mask.setVal(0);

    int cnt = 0;

    const auto len  = amrex::length(this->domain);
    const auto lo   = amrex::lbound(this->domain);
    const auto dp   =      view(lo, comp);
    const auto mp   = mask.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                if (dp(i,j,k,0) > val) {
                    mp(i,j,k,0) = 1;
                    ++cnt;
                }
            }
        }
    }

    return cnt;
}

template <class T, class Allocator>
int
BaseFab<T,Allocator>::maskGE (BaseFab<int>& mask,
                              T             val,
                              int           comp) const noexcept
{
    mask.resize(this->domain,1);
    mask.setVal(0);

    int cnt = 0;

    const auto len  = amrex::length(this->domain);
    const auto lo   = amrex::lbound(this->domain);
    const auto dp   =      view(lo, comp);
    const auto mp   = mask.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                if (dp(i,j,k,0) >= val) {
                    mp(i,j,k,0) = 1;
                    ++cnt;
                }
            }
        }
    }

    return cnt;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::atomicAdd (const BaseFab<T,Allocator>& x) noexcept
{
    Box ovlp(this->domain);
    ovlp &= x.domain;
    return ovlp.ok() ? this->atomicAdd(x,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::saxpy (T a, const BaseFab<T,Allocator>& x,
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

    const auto len = amrex::length(destbox);
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const auto dp  =   view(dlo, destcomp);
    const auto sp  = x.view(slo, srccomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) += a * sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::saxpy (T a, const BaseFab<T,Allocator>& x) noexcept
{
    Box ovlp(this->domain);
    ovlp &= x.domain;
    return ovlp.ok() ? saxpy(a,x,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::xpay (T a, const BaseFab<T,Allocator>& x,
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

    const auto len = amrex::length(destbox);
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const auto dp  =   view(dlo, destcomp);
    const auto sp  = x.view(slo, srccomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) = sp(i,j,k,n) + a * dp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::addproduct (const Box& destbox, int destcomp, int numcomp,
                                  const BaseFab<T,Allocator>& src1, int comp1,
                                  const BaseFab<T,Allocator>& src2, int comp2) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(   comp1 >= 0 &&    comp1+numcomp <= src1.nComp());
    BL_ASSERT(   comp2 >= 0 &&    comp2+numcomp <= src2.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <=      nComp());

    const auto len = amrex::length(destbox);
    const auto lo  = amrex::lbound(destbox);
    const auto dp  =      view(lo, destcomp);
    const auto sp1 = src1.view(lo, comp1);
    const auto sp2 = src2.view(lo, comp2);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) += sp1(i,j,k,n) * sp2(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::linComb (const BaseFab<T,Allocator>& f1, const Box& b1, int comp1,
                               const BaseFab<T,Allocator>& f2, const Box& b2, int comp2,
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

    const auto len = amrex::length(b);
    const auto dlo = amrex::lbound(b);
    const auto slo1 = amrex::lbound(b1);
    const auto slo2 = amrex::lbound(b2);
    const auto dp  =    view(dlo,   comp);
    const auto sp1 = f1.view(slo1, comp1);
    const auto sp2 = f2.view(slo2, comp2);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) = alpha*sp1(i,j,k,n) + beta*sp2(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::dot (const Box& xbx, int xcomp,
                           const BaseFab<T,Allocator>& y, const Box& ybx, int ycomp,
                           int numcomp) const noexcept
{
    BL_ASSERT(xbx.ok());
    BL_ASSERT(box().contains(xbx));
    BL_ASSERT(y.box().contains(ybx));
    BL_ASSERT(xbx.sameSize(ybx));
    BL_ASSERT(xcomp >= 0 && xcomp+numcomp <=   nComp());
    BL_ASSERT(ycomp >= 0 && ycomp+numcomp <= y.nComp());

    T r = 0;

    const auto len = amrex::length(xbx);
    const auto xlo = amrex::lbound(xbx);
    const auto ylo = amrex::lbound(ybx);
    const auto xp  =   view(xlo, xcomp);
    const auto yp  = y.view(ylo, ycomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                T const * const AMREX_RESTRICT xpr = &(xp(0,j,k,n));
                T const * const AMREX_RESTRICT ypr = &(yp(0,j,k,n));
                for (int i = 0; i < len.x; ++i) {
//                    r += xp(i,j,k,n) * yp(i,j,k,n);
                    r += xpr[i] * ypr[i];
                }
            }
        }
    }

    return r;
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::dotmask (const BaseFab<int>& mask, const Box& xbx, int xcomp,
                               const BaseFab<T,Allocator>& y, const Box& ybx, int ycomp,
                               int numcomp) const noexcept
{
    BL_ASSERT(xbx.ok());
    BL_ASSERT(box().contains(xbx));
    BL_ASSERT(y.box().contains(ybx));
    BL_ASSERT(xbx.sameSize(ybx));
    BL_ASSERT(xcomp >= 0 && xcomp+numcomp <=   nComp());
    BL_ASSERT(ycomp >= 0 && ycomp+numcomp <= y.nComp());

    T r = 0;

    const auto len = amrex::length(xbx);
    const auto xlo = amrex::lbound(xbx);
    const auto ylo = amrex::lbound(ybx);
    const auto xp  =      view(xlo, xcomp);
    const auto yp  =    y.view(ylo, ycomp);
    const auto mp  = mask.view(xlo);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    int m = static_cast<int>(static_cast<bool>(mp(i,j,k,0)));
                    r += xp(i,j,k,n) * yp(i,j,k,n) * m;
                }
            }
        }
    }

    return r;
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::sum (int comp, int numcomp) const noexcept
{
    return this->sum(this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::sum (const Box& subbox, int comp, int numcomp) const noexcept
{
    return this->sum(subbox, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::negate (int comp, int numcomp) noexcept
{
    return this->negate(this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::negate (const Box& b, int comp, int numcomp) noexcept
{
    return this->negate(b, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::invert (T r, int comp, int numcomp) noexcept
{
    return this->invert(r, this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::invert (T r, const Box& b, int comp, int numcomp) noexcept
{
    return this->invert(r, b, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::plus (T r, int comp, int numcomp) noexcept
{
    return this->plus(r, this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::plus (T r, const Box& b, int comp, int numcomp) noexcept
{
    return this->plus(r, b, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::plus (const BaseFab<T,Allocator>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->plus(src, this->domain, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::atomicAdd (const BaseFab<T,Allocator>& src,
                                 int               srccomp,
                                 int               destcomp,
                                 int               numcomp) noexcept
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->atomicAdd(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::plus (const BaseFab<T,Allocator>& src,
                            const Box&        subbox,
                            int               srccomp,
                            int               destcomp,
                            int               numcomp) noexcept
{
    return this->plus(src, subbox, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::atomicAdd (const BaseFab<T,Allocator>& src,
                                 const Box&        subbox,
                                 int               srccomp,
                                 int               destcomp,
                                 int               numcomp) noexcept
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    ovlp &= subbox;
    return ovlp.ok() ? this->atomicAdd(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::plus (const BaseFab<T,Allocator>& src,
                            const Box&        srcbox,
                            const Box&        destbox,
                            int               srccomp,
                            int               destcomp,
                            int               numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    const auto len = amrex::length(destbox);
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const auto dp  =     view(dlo, destcomp);
    const auto sp  = src.view(slo, srccomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) += sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::atomicAdd (const BaseFab<T,Allocator>& src,
                                 const Box&        srcbox,
                                 const Box&        destbox,
                                 int               srccomp,
                                 int               destcomp,
                                 int               numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    const auto len = amrex::length(destbox);
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const auto dp  =     view(dlo, destcomp);
    const auto sp  = src.view(slo, srccomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
#ifdef _OPENMP
                T       * dp1d = &(dp(0,j,k,n));
                T const * sp1d = &(sp(0,j,k,n));
#else
                T       * AMREX_RESTRICT dp1d = &(dp(0,j,k,n));
                T const * AMREX_RESTRICT sp1d = &(sp(0,j,k,n));
                AMREX_PRAGMA_SIMD
#endif
                for (int i = 0; i < len.x; ++i) {
#ifndef __CUDA_ARCH__
#ifdef _OPENMP
#pragma omp atomic update
#endif
                    dp1d[i] += sp1d[i];
#else
                    amrex::Gpu::Atomic::Add(&(dp1d[i]), sp1d[i]);
#endif
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::minus (const BaseFab<T,Allocator>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->minus(src, this->domain, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::minus (const BaseFab<T,Allocator>& src, const Box& subbox, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->minus(src, subbox, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::minus (const BaseFab<T,Allocator>& src,
                             const Box&        srcbox,
                             const Box&        destbox,
                             int               srccomp,
                             int               destcomp,
                             int               numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    const auto len = amrex::length(destbox);
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const auto dp  =     view(dlo, destcomp);
    const auto sp  = src.view(slo, srccomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) -= sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::mult (T r, int comp, int numcomp) noexcept
{
    return this->mult(r, this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::mult (T r, const Box& b, int comp, int numcomp) noexcept
{
    return this->mult(r, b, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::mult (const BaseFab<T,Allocator>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->mult(src, this->domain, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::mult (const BaseFab<T,Allocator>& src, const Box& subbox, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->mult(src, subbox, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::mult (const BaseFab<T,Allocator>& src,
                            const Box&        srcbox,
                            const Box&        destbox,
                            int               srccomp,
                            int               destcomp,
                            int               numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    const auto len = amrex::length(destbox);
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const auto dp  =     view(dlo, destcomp);
    const auto sp  = src.view(slo, srccomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) *= sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::divide (T r, int comp, int numcomp) noexcept
{
    return this->divide(r, this->domain, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::divide (T r, const Box& b, int comp, int numcomp) noexcept
{
    return this->divide(r, b, DestComp{comp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::divide (const BaseFab<T,Allocator>& src, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->divide(src, this->domain, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::divide (const BaseFab<T,Allocator>& src, const Box& subbox, int srccomp, int destcomp, int numcomp) noexcept
{
    return this->divide(src, subbox, SrcComp{srccomp}, DestComp{destcomp}, NumComps{numcomp});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::divide (const BaseFab<T,Allocator>& src,
                              const Box&        srcbox,
                              const Box&        destbox,
                              int               srccomp,
                              int               destcomp,
                              int               numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    const auto len = amrex::length(destbox);
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const auto dp  =     view(dlo, destcomp);
    const auto sp  = src.view(slo, srccomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) /= sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::protected_divide (const BaseFab<T,Allocator>& src) noexcept
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->protected_divide(src,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::protected_divide (const BaseFab<T,Allocator>& src,
                                        int               srccomp,
                                        int               destcomp,
                                        int               numcomp) noexcept
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->protected_divide(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::protected_divide (const BaseFab<T,Allocator>& src,
                                        const Box&        subbox,
                                        int               srccomp,
                                        int               destcomp,
                                        int               numcomp) noexcept
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    ovlp &= subbox;
    return ovlp.ok() ? this->protected_divide(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::protected_divide (const BaseFab<T,Allocator>& src,
                                        const Box&        srcbox,
                                        const Box&        destbox,
                                        int               srccomp,
                                        int               destcomp,
                                        int               numcomp) noexcept
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    const auto len = amrex::length(destbox);
    const auto dlo = amrex::lbound(destbox);
    const auto slo = amrex::lbound(srcbox);
    const auto dp  =     view(dlo, destcomp);
    const auto sp  = src.view(slo, srccomp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    if (sp(i,j,k,n) != 0) {
                        dp(i,j,k,n) /=  sp(i,j,k,n);
                    }
                }
            }
        }
    }

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

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::linInterp (const BaseFab<T,Allocator>& f1,
                                 const Box&        b1,
                                 int               comp1,
                                 const BaseFab<T,Allocator>& f2,
                                 const Box&        b2,
                                 int               comp2,
                                 Real              t1,
                                 Real              t2,
                                 Real              t,
                                 const Box&        b,
                                 int               comp,
                                 int               numcomp) noexcept
{
    Real alpha = (t2-t)/(t2-t1);
    Real beta = (t-t1)/(t2-t1);
    return linComb(f1,b1,comp1,f2,b2,comp2,alpha,beta,b,comp,numcomp);
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::linInterp (const BaseFab<T,Allocator>& f1,
                                 int               comp1,
                                 const BaseFab<T,Allocator>& f2,
                                 int               comp2,
                                 Real              t1,
                                 Real              t2,
                                 Real              t,
                                 const Box&        b,
                                 int               comp,
                                 int               numcomp) noexcept
{
  Real tol = 1.0e-16;
  if(std::abs(t2-t1) > tol)
  {
    Real alpha = (t2-t)/(t2-t1);
    Real beta = (t-t1)/(t2-t1);

    return linComb(f1,b,comp1,f2,b,comp2,alpha,beta,b,comp,numcomp);
  }
  else
  {
      copy(f1,b,comp1,b,comp,numcomp);
  }
  return *this;
}

//
// New interfaces
//

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setVal (T val) noexcept
{
    this->setVal(val, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setVal (T x, Box bx, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, dcomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) = x;
                }
            }
        }
    }
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setValIf (T val, const BaseFab<int>& mask) noexcept
{
    this->setValIf(val, this->domain, mask, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setValIf (T val, Box bx, const BaseFab<int>& mask, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, dcomp.i);
    const auto mp  = mask.view(lo);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    if (mp(i,j,k,0)) { dp(i,j,k,n) = val; }
                }
            }
        }
    }
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setValIfNot (T val, const BaseFab<int>& mask) noexcept
{
    this->setValIfNot(val, this->domain, mask, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setValIfNot (T val, Box bx, const BaseFab<int>& mask, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, dcomp.i);
    const auto mp  = mask.view(lo);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    if (!mp(i,j,k,0)) { dp(i,j,k,n) = val; }
                }
            }
        }
    }
}

template <class T, class Allocator>
void
BaseFab<T,Allocator>::setComplement (T x, const Box& bx, DestComp dcomp, NumComps ncomp) noexcept
{
    const BoxList b_lst = amrex::boxDiff(this->domain,bx);
    for (auto const& b : b_lst) {
        this->setVal(x, b, dcomp, ncomp);
    }
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::copy (const BaseFab<T,Allocator>& src) noexcept
{
    this->copy(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::copy (const BaseFab<T,Allocator>& src, Box bx,
                            SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    bx &= src.domain;

    const auto len = amrex::length(bx);
    const auto  lo = amrex::lbound(bx);
    const auto dp  =     view(lo, dcomp.i);
    const auto sp  = src.view(lo, scomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) = sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::plus (T val) noexcept
{
    return this->plus(val, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::operator+= (T val) noexcept
{
    return this->plus(val, this->domain,  DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::plus (T val, Box bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, dcomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) += val;
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::plus (const BaseFab<T,Allocator>& src) noexcept
{
    return this->plus(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::operator+= (const BaseFab<T,Allocator>& src) noexcept
{
    return this->plus(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::plus (const BaseFab<T,Allocator>& src, Box bx,
                            SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    bx &= src.domain;

    const auto len = amrex::length(bx);
    const auto  lo = amrex::lbound(bx);
    const auto dp  =     view(lo, dcomp.i);
    const auto sp  = src.view(lo, scomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) += sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::minus (T val) noexcept
{
    return this->minus(val, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::operator-= (T val) noexcept
{
    return this->minus(val, this->domain,  DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::minus (T val, Box bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, dcomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) -= val;
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::minus (const BaseFab<T,Allocator>& src) noexcept
{
    return this->minus(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::operator-= (const BaseFab<T,Allocator>& src) noexcept
{
    return this->minus(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::minus (const BaseFab<T,Allocator>& src, Box bx,
                             SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    bx &= src.domain;

    const auto len = amrex::length(bx);
    const auto  lo = amrex::lbound(bx);
    const auto dp  =     view(lo, dcomp.i);
    const auto sp  = src.view(lo, scomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) -= sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::mult (T val) noexcept
{
    return this->mult(val, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::operator*= (T val) noexcept
{
    return this->mult(val, this->domain,  DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::mult (T val, Box bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, dcomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) *= val;
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::mult (const BaseFab<T,Allocator>& src) noexcept
{
    return this->mult(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::operator*= (const BaseFab<T,Allocator>& src) noexcept
{
    return this->mult(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::mult (const BaseFab<T,Allocator>& src, Box bx,
                            SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    bx &= src.domain;

    const auto len = amrex::length(bx);
    const auto  lo = amrex::lbound(bx);
    const auto dp  =     view(lo, dcomp.i);
    const auto sp  = src.view(lo, scomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) *= sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::divide (T val) noexcept
{
    return this->divide(val, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::operator/= (T val) noexcept
{
    return this->divide(val, this->domain,  DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::divide (T val, Box bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, dcomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) /= val;
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::divide (const BaseFab<T,Allocator>& src) noexcept
{
    return this->divide(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::operator/= (const BaseFab<T,Allocator>& src) noexcept
{
    return this->divide(src, this->domain, SrcComp{0}, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::divide (const BaseFab<T,Allocator>& src, Box bx,
                              SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    bx &= src.domain;

    const auto len = amrex::length(bx);
    const auto  lo = amrex::lbound(bx);
    const auto dp  =     view(lo, dcomp.i);
    const auto sp  = src.view(lo, scomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) /= sp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::negate () noexcept
{
    return this->negate(this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::negate (const Box& bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, dcomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                   dp(i,j,k,n) = -dp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::invert (T r) noexcept
{
    return this->invert(r, this->domain, DestComp{0}, NumComps{this->nvar});
}

template <class T, class Allocator>
BaseFab<T,Allocator>&
BaseFab<T,Allocator>::invert (T r, const Box& bx, DestComp dcomp, NumComps ncomp) noexcept
{
    BL_ASSERT(dcomp.i >= 0 && dcomp.i + ncomp.n <= this->nvar);

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, dcomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) = r / dp(i,j,k,n);
                }
            }
        }
    }

    return *this;
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::sum (const Box& bx, DestComp dcomp, NumComps ncomp) const noexcept
{
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    T r = 0;

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, dcomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                T const * const AMREX_RESTRICT dpr = &(dp(0,j,k,n));
                for (int i = 0; i < len.x; ++i) {
                    r += dpr[i];
                }
            }
        }
    }

    return r;
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::dot (const BaseFab<T,Allocator>& src, const Box& bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) const noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    T r = 0;

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  =     view(lo, dcomp.i);
    const auto sp  = src.view(lo, scomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                T const * const AMREX_RESTRICT dpr = &(dp(0,j,k,n));
                T const * const AMREX_RESTRICT spr = &(sp(0,j,k,n));
                for (int i = 0; i < len.x; ++i) {
                    r += dpr[i] * spr[i];
                }
            }
        }
    }

    return r;
}

template <class T, class Allocator>
T
BaseFab<T,Allocator>::dotmask (const BaseFab<T,Allocator>& src, const Box& bx, const BaseFab<int>& mask,
                               SrcComp scomp, DestComp dcomp, NumComps ncomp) const noexcept
{
    AMREX_ASSERT(this->domain.sameType(src.domain));
    AMREX_ASSERT(this->domain.sameType(mask.domain));
    AMREX_ASSERT(scomp.i >= 0 && scomp.i+ncomp.n <= src.nvar);
    AMREX_ASSERT(dcomp.i >= 0 && dcomp.i+ncomp.n <= this->nvar);

    T r = 0;

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  =      view(lo, dcomp.i);
    const auto sp  =  src.view(lo, scomp.i);
    const auto mp  = mask.view(lo);

    for (int n = 0; n < ncomp.n; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                T const * const AMREX_RESTRICT dpr = &(dp(0,j,k,n));
                T const * const AMREX_RESTRICT spr = &(sp(0,j,k,n));
                int const * const AMREX_RESTRICT mpr = &(mp(0,j,k));
                for (int i = 0; i < len.x; ++i) {
                    int m = static_cast<int>(static_cast<bool>(mpr[i]));
                    r += dpr[i] * spr[i] * m;
                }
            }
        }
    }

    return r;
}

}

#endif /*BL_BASEFAB_H*/
