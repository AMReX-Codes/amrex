#ifndef AMREX_BASEFAB_H_
#define AMREX_BASEFAB_H_

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <array>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Extension.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxList.H>
#include <AMReX_BArena.H>
#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_MakeType.H>

#include <AMReX_Gpu.H>

#ifdef USE_PERILLA
#include <LocalConnection.H>
#include <RemoteConnection.H>
#endif

namespace amrex
{
    extern long private_total_bytes_allocated_in_fabs;     // total bytes at any given time
    extern long private_total_bytes_allocated_in_fabs_hwm; // high-water-mark over a given interval
    extern long private_total_cells_allocated_in_fabs;     // total cells at any given time
    extern long private_total_cells_allocated_in_fabs_hwm; // high-water-mark over a given interval
#ifdef _OPENMP
#pragma omp threadprivate(private_total_bytes_allocated_in_fabs)
#pragma omp threadprivate(private_total_bytes_allocated_in_fabs_hwm)
#pragma omp threadprivate(private_total_cells_allocated_in_fabs)
#pragma omp threadprivate(private_total_cells_allocated_in_fabs_hwm)
#endif

    long TotalBytesAllocatedInFabs();
    long TotalBytesAllocatedInFabsHWM();
    long TotalCellsAllocatedInFabs();
    long TotalCellsAllocatedInFabsHWM();
    void ResetTotalBytesAllocatedInFabsHWM();
    void update_fab_stats (long n, long s, std::size_t szt);

    void BaseFab_Initialize ();
    void BaseFab_Finalize ();

template <typename T>
struct FabView
{
    T* p;
    long jstride;
    long kstride;
    long nstride;
#ifdef AMREX_DEBUG
    int iend, jend, kend, nend;
#endif
    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE T& operator() (int i, int j, int k, int n=0) const {
#ifdef AMREX_DEBUG
        AMREX_ASSERT(i<iend && j<jend && k<kend && n<nend);
#endif
        return *(p + (i+j*jstride+k*kstride+n*nstride));
    }
};

template <typename T>
struct BaseFabData : public Gpu::Managed
{
    AMREX_GPU_HOST_DEVICE
    BaseFabData (Box const& a_domain, IntVect const& a_dlen,
                 int a_nvar, long a_numpts, long a_truesize,
                 T* a_dptr, bool a_ptr_owner, bool a_shared_memory)
        : domain(a_domain), dlen(a_dlen),
          nvar(a_nvar), numpts(a_numpts), truesize(a_truesize),
          dptr(a_dptr), ptr_owner(a_ptr_owner), shared_memory(a_shared_memory)
        {}
    BaseFabData () = default;
    Box     domain;   // My index space.
    IntVect dlen;     // Length of domain in each direction.
    int     nvar;     // Number components.
    long    numpts;   // Cached number of points in FAB.
    long    truesize; // nvar*numpts that was allocated on heap.
    T*      dptr;     // The data pointer.
    bool    ptr_owner;// Owner of T*?
    bool    shared_memory;  // Is the memory allocated in shared memory?
    void setOwner (bool b) { ptr_owner = b; }
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
class BaseFab : public BaseFabData<T>
{
public:

    template <class U> friend class BaseFab;

    typedef T value_type;
    //! Construct an empty BaseFab, which must be resized (see BaseFab::resize) before use.
    BaseFab ();

    //!  Make BaseFab with desired domain (box) and number of components.
    explicit BaseFab (const Box& bx,
                      int        n = 1,
		      bool       alloc = true,
		      bool       shared = false);

    BaseFab (const BaseFab<T>& rhs, MakeType make_type, int scomp, int ncomp);

#ifdef AMREX_USE_GPU
    BaseFab (const BaseFab<T>& rhs, MakeType make_type);
#endif

    /**
     * \brief Create an NON-OWNING BaseFab.  Thus BaseFab is not
     * responsible for memory management.  And it's caller's responsibility that
     * p points to a chunk of memory large enough.
     */
    BaseFab (const Box& bx, int ncomp, T* p);

    //! The destructor deletes the array memory.
    ~BaseFab ();

    BaseFab (const BaseFab<T>& rhs) = delete;
    BaseFab<T>& operator= (const BaseFab<T>& rhs) = delete;
    BaseFab<T>& operator= (BaseFab<T>&& rhs) = delete;

    BaseFab (BaseFab<T>&& rhs) noexcept;

    AMREX_GPU_HOST_DEVICE
    BaseFab& operator= (const T&);

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

    /**
     * \brief The function returns the BaseFab to the invalid state.
     * The memory is freed.
     */
    void clear ();

    //! Returns how many bytes used
    AMREX_GPU_HOST_DEVICE
    std::size_t nBytes () const { return this->truesize*sizeof(T); }

    static bool preAllocatable () { return true; }

    static bool isCopyOMPSafe () { return true; }

    //! Returns bytes used in the Box for those components
    AMREX_GPU_HOST_DEVICE
    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const
        { return bx.numPts() * sizeof(T) * ncomps; }

    //! Returns the number of components
    AMREX_GPU_HOST_DEVICE
    int nComp () const { return this->nvar; }

    //! for calls to fortran.
    AMREX_GPU_HOST_DEVICE
    const int* nCompPtr() const {
        return &(this->nvar);
    }

    //! Returns the number of points
    AMREX_GPU_HOST_DEVICE
    long nPts () const { return this->numpts; }

    //! Returns the total number of points of all components
    AMREX_GPU_HOST_DEVICE
    long size () const { return this->nvar*this->numpts; }

    //! Returns the domain (box) where the array is defined
    AMREX_GPU_HOST_DEVICE
    const Box& box () const { return this->domain; }

    /**
    * \brief Returns a pointer to an array of SPACEDIM integers
    * giving the length of the domain in each direction
    */
    AMREX_GPU_HOST_DEVICE
    const int* length () const { return this->dlen.getVect(); }

    /**
    * \brief Returns the lower corner of the domain
    * See class Box for analogue.
    */
    AMREX_GPU_HOST_DEVICE
    const IntVect& smallEnd () const { return this->domain.smallEnd(); }

    //!  Returns the upper corner of the domain.  See class Box for analogue.
    AMREX_GPU_HOST_DEVICE
    const IntVect& bigEnd () const { return this->domain.bigEnd(); }

    /**
    * \brief Returns the lower corner of the domain.

    *Instead of returning them in the form of INTVECTs, as in smallEnd and
    * bigEnd, it returns the values as a pointer to an array of
    * constant integers.  This is useful when interfacing to
    * Fortran subroutines.
    */
    AMREX_GPU_HOST_DEVICE
    const int* loVect () const { return this->domain.loVect(); }

    /**
    * \brief Returns the upper corner of the domain.

    *Instead of returning them in the form of INTVECTs, as in smallEnd and
    * bigEnd, it returns the values as a pointer to an array of
    * constant integers.  This is useful when interfacing to
    * Fortran subroutines.
    */
    AMREX_GPU_HOST_DEVICE
    const int* hiVect () const { return this->domain.hiVect(); }

    /**
    * \brief Returns true if the domain of fab is totally contained within
    * the domain of this BaseFab.
    */
    AMREX_GPU_HOST_DEVICE
    bool contains (const BaseFab<T>& fab) const
    {
        return box().contains(fab.box()) && this->nvar <= fab.nvar;
    }

    /**
    * \brief Returns true if bx is totally contained
    * within the domain of this BaseFab.
    */
    AMREX_GPU_HOST_DEVICE
    bool contains (const Box& bx) const { return box().contains(bx); }

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
    T* dataPtr (int n = 0) { AMREX_ASSERT(!(this->dptr == 0)); return &(this->dptr[n*this->numpts]); }

    //! Same as above except works on const FABs.
    AMREX_GPU_HOST_DEVICE
    const T* dataPtr (int n = 0) const { AMREX_ASSERT(!(this->dptr == 0)); return &(this->dptr[n*this->numpts]); }

    AMREX_GPU_HOST_DEVICE
    T* dataPtr (const IntVect& iv, int n = 0);
    AMREX_GPU_HOST_DEVICE
    const T* dataPtr (const IntVect& iv, int n = 0) const;

    void setPtr (T* p, long sz) { AMREX_ASSERT(this->dptr == 0 && this->truesize == 0); this->dptr = p; this->truesize = sz; }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T const> view (int n = 0) const
    {
        return view(this->domain.loVect3d(),n);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T> view (int n = 0)
    {
        return view(this->domain.loVect3d(),n);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T const> view (const IntVect& iv, int n = 0) const
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
    FabView<T> view (const IntVect& iv, int n = 0)
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
    FabView<T const> view (const Box& subbox, int n = 0) const
    {
        return view(subbox.loVect3d(),n);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T> view (const Box& subbox, int n = 0)
    {
        return view(subbox.loVect3d(),n);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_INLINE
    FabView<T const> view (const GpuArray<int,3>& slo, int n = 0) const
    {
        const auto& len = this->domain.length3d();
        long jstride = len[0];
        long kstride = jstride*len[1];
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
    FabView<T> view (const GpuArray<int,3>& slo, int n = 0)
    {
        const auto& len = this->domain.length3d();
        long jstride = len[0];
        long kstride = jstride*len[1];
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
    FabView<T const> view (const Dim3& slo, int n = 0) const
    {
        const auto len = amrex::length(this->domain);
        const auto dlo = amrex::lbound(this->domain);
        long jstride = len.x;
        long kstride = jstride*len.y;
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
    FabView<T> view (const Dim3& slo, int n = 0)
    {
        const auto len = amrex::length(this->domain);
        const auto dlo = amrex::lbound(this->domain);
        long jstride = len.x;
        long kstride = jstride*len.y;
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
    bool isAllocated () const { return this->dptr != 0; }

    /**
    * \brief Returns a reference to the Nth component value
    * defined at position p in the domain.  This operator may be
    * inefficient if the C++ compiler is unable to optimize the
    * C++ code.
    */
    AMREX_GPU_HOST_DEVICE
    T& operator() (const IntVect& p, int N);

    //! Same as above, except returns component 0.
    AMREX_GPU_HOST_DEVICE
    T& operator() (const IntVect& p);

    //! Same as above except works on const FABs.
    AMREX_GPU_HOST_DEVICE
    const T& operator() (const IntVect& p, int N) const;

    //! Same as above, except returns component 0.
    AMREX_GPU_HOST_DEVICE
    const T& operator() (const IntVect& p) const;

    /**
    * \brief This function puts numcomp component values, starting at
    * component N, from position pos in the domain into array data,
    * that must be allocated by the user.
    */
    AMREX_GPU_HOST_DEVICE
    void getVal (T*             data,
                 const IntVect& pos,
                 int            N,
                 int            numcomp) const;
    //! Same as above, except that starts at component 0 and copies all comps.
    AMREX_GPU_HOST_DEVICE
    void getVal (T*             data,
                 const IntVect& pos) const;
    /**
    * \brief The setVal functions set sub-regions in the BaseFab to a
    * constant value.  This most general form specifies the sub-box,
    * the starting component number, and the number of components
    * to be set.
    */
    AMREX_GPU_HOST_DEVICE
    void setVal (T          x,
                 const Box& bx,
                 int        nstart,
                 int        ncomp);
    //! Same as above, except the number of modified components is one. N is the component to be modified.
    AMREX_GPU_HOST_DEVICE
    void setVal (T          x,
                 const Box& bx,
                 int        N = 0);
    //! Same as above, except the sub-box defaults to the entire domain.
    AMREX_GPU_HOST_DEVICE
    void setVal (T   x,
                 int N);
    //! Same as above, except all components are set.
    AMREX_GPU_HOST_DEVICE
    void setVal (T x);

    AMREX_GPU_HOST_DEVICE
    void setValIfNot (T x, const Box& bx, const BaseFab<int>& mask, int nstart, int ncomp);

    /**
    * \brief This function is analogous to the fourth form of
    * setVal above, except that instead of setting values on the
    * Box b, values are set on the complement of b in the domain.
    */
    void setComplement (T          x,
                        const Box& b,
                        int        ns,
                        int        num);
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
    * destcomp. Note that the actual copy is made by the
    * function performCopy of this class.  The results are
    * UNDEFINED if the src and dest are the same and the srcbox and
    * destbox overlap.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& copy (const BaseFab<T>& src,
                      const Box&        srcbox,
                      int               srccomp,
                      const Box&        destbox,
                      int               destcomp,
                      int               numcomp);

    /**
    * \brief As above, except the destination Box and the source Box
    * are taken to be the entire domain of the destination.   A copy
    * of the intersecting region is performed.  Note that the actual
    * copy is made by the function performCopy() of this
    * class.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& copy (const BaseFab<T>& src,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp = 1);
    /**
    * \brief As above, except that the destination Box is specified,
    * but the source Box is taken to the equal to the source
    * Box, and all components of the destination BaseFab are
    * copied.  Note that the actual copy is made by the
    * function performCopy of this class.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& copy (const BaseFab<T>& src,
                      const Box&        destbox);

    /**
    * As above, except that the destbox defaults to the entire domain
    * of the destination BaseFab, and all components are copied.
    * Note that the actual copy is made by the function
    * performCopy() of this class.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& copy (const BaseFab<T>& src);

    //! Copy from the srcbox of this Fab to raw memory and return the number of bytes copied
    AMREX_GPU_HOST_DEVICE
    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const;

    AMREX_GPU_HOST_DEVICE
    //! Copy from raw memory to the dstbox of this Fab and return the number of bytes copied
    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);

    AMREX_GPU_HOST_DEVICE
    //! Add from raw memory to the dstbox of this Fab and return the number of bytes copied
    std::size_t addFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);

    /**
    * \brief Perform shifts upon the domain of the BaseFab. They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& shift (const IntVect& v);
    /**
    * \brief Perform shifts upon the domain of the BaseFab.  They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& shift (int idir,
                       int n_cell);
    /**
    * \brief Perform shifts upon the domain of the BaseFab.  They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& shiftHalf (int dir,
                           int num_halfs);
    /**
    * \brief Perform shifts upon the domain of the BaseFab. They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& shiftHalf (const IntVect& num_halfs);

    AMREX_GPU_HOST_DEVICE
    Real norminfmask (const Box& subbox, const BaseFab<int>& mask, int scomp=0, int ncomp=1) const;

    /**
    * \brief Compute the Lp-norm of this FAB using components (scomp : scomp+ncomp-1).
    *   p < 0  -> ERROR
    *   p = 0  -> infinity norm (max norm)
    *   p = 1  -> sum of ABS(FAB)
    */
    AMREX_GPU_HOST_DEVICE
    Real norm (int p,
               int scomp = 0,
               int ncomp = 1) const;

    //! Same as above except only on given subbox.
    AMREX_GPU_HOST_DEVICE
    Real norm (const Box& subbox,
               int        p,
               int        scomp = 0,
               int        ncomp = 1) const;
    //!Compute absolute value for all components of this FAB.
    AMREX_GPU_HOST_DEVICE
    void abs ();
    //! Same as above except only for components (comp: comp+numcomp-1)
    AMREX_GPU_HOST_DEVICE
    void abs (int comp,
              int numcomp=1);
    /**
    * \brief Calculate abs() on subbox for given component range.
    */
    AMREX_GPU_HOST_DEVICE
    void abs (const Box& subbox,
              int        comp = 0,
              int        numcomp=1);
    /**
    * \return Minimum value of given component.
    */
    AMREX_GPU_HOST_DEVICE
    T min (int comp = 0) const;
    /**
    * \return Minimum value of given component in given subbox.
    */
    AMREX_GPU_HOST_DEVICE
    T min (const Box& subbox,
           int        comp = 0) const;
    /**
    * \return Maximum value of given component.
    */
    AMREX_GPU_HOST_DEVICE
    T max (int comp = 0) const;
    /**
    * \return Maximum value of given component in given subbox.
    */
    AMREX_GPU_HOST_DEVICE
    T max (const Box& subbox,
           int        comp = 0) const;
    /**
    * \return Maximum of the absolute value of given component.
    */
    AMREX_GPU_HOST_DEVICE
    T maxabs (int comp = 0) const;
    /**
    * \return Maximum of the absolute value of given component in given subbox.
    */
    AMREX_GPU_HOST_DEVICE
    T maxabs (const Box& subbox,
              int        comp = 0) const;

    /**
    * \return location of a cell containing the specified value
    * given subbox. Returns IntVect outside box if value not present.
    */
    AMREX_GPU_HOST_DEVICE
    IntVect indexFromValue (Real       value,
                            const Box& subbox,
                            int        comp = 0) const;

    /**
    * \return location of minimum value in given component.
    */
    AMREX_GPU_HOST_DEVICE
    IntVect minIndex (int comp = 0) const;
    /**
    * \return location of minimum value in given component in
    * given subbox.
    */
    AMREX_GPU_HOST_DEVICE
    IntVect minIndex (const Box& subbox,
                      int        comp = 0) const;
    /**
    * \return return mininum value and location to allow
    * efficient looping over multiple boxes.
    */
    AMREX_GPU_HOST_DEVICE
    void  minIndex (const Box&     subbox,
                          Real&    min_val,
                          IntVect& min_idx,
                          int      comp = 0) const;

    /**
    * \return location of maximum value in given component.
    */
    AMREX_GPU_HOST_DEVICE
    IntVect maxIndex (int comp = 0) const;
    /**
    * \return location of maximum value in given component in given
    * subbox.
    */
    AMREX_GPU_HOST_DEVICE
    IntVect maxIndex (const Box& subbox,
                      int        comp = 0) const;
    /**
    * \return return maximum value and location to allow
    * efficient looping over multiple boxes.
    */
    AMREX_GPU_HOST_DEVICE
    void  maxIndex (const Box&     subbox,
                          Real&    max_value,
                          IntVect& max_index,
                          int      comp = 0) const;


    /**
    * \brief Compute mask array with value of 1 in cells where
    * BaseFab has value less than val, 0 otherwise.
    * mask is resized by this function.
    * The number of cells marked with 1 returned.
    */
    AMREX_GPU_HOST_DEVICE
    int maskLT (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
    //! Same as above except mark cells with value less than or equal to val.
    AMREX_GPU_HOST_DEVICE
    int maskLE (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;

    //! Same as above except mark cells with value equal to val.
    AMREX_GPU_HOST_DEVICE
    int maskEQ (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
    //! Same as above except mark cells with value greater than val.
    AMREX_GPU_HOST_DEVICE
    int maskGT (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
    //! Same as above except mark cells with value greater than or equal to val.
    AMREX_GPU_HOST_DEVICE
    int maskGE (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
    //! Returns sum of given component of FAB state vector.
    AMREX_GPU_HOST_DEVICE
    T sum (int comp,
           int numcomp = 1) const;
    //! Compute sum of given component of FAB state vector in given subbox.
    AMREX_GPU_HOST_DEVICE
    T sum (const Box& subbox,
           int        comp,
           int        numcomp = 1) const;
    //! Most general version, specify subbox and which components.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& invert (T          v,
                        const Box& subbox,
                        int        comp=0,
                        int        numcomp=1);
    //! As above except on entire domain.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& invert (T   v,
                        int comp,
                        int numcomp=1);
    //! As above except on entire domain, all components.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& invert (T v);

    //! Negate BaseFab, most general.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& negate (const Box& subbox,
                        int        comp=0,
                        int        numcomp=1);
    //! As above, except on entire domain.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& negate (int comp,
                        int numcomp=1);
    //! As above, except on entire domain and all components.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& negate ();

    //! Scalar addition (a[i] <- a[i] + r), most general.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& plus (T          r,
                      const Box& b,
                      int        comp=0,
                      int        numcomp=1);

    //! As above, except on entire domain.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& plus (T   r,
                      int comp,
                      int numcomp=1);
    //! As above, except on entire domain and all components.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& plus (T r);

    //! Addition in place.  This will often be more efficient than  making new BaseFab for result.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& operator+= (T r);

    //! FAB addition (a[i] <- a[i] + b[i]) in place.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& operator+= (const BaseFab<T>& f);

    //! FAB addition (a[i] <- a[i] + b[i]). The same as += operator.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& plus (const BaseFab<T>& src);

    /**
    * \brief Add src components (srccomp:srccomp+numcomp-1) to
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& plus (const BaseFab<T>& src,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);
    /**
    * \brief Same as above except addition is restricted to intersection
    * of subbox and src FAB. NOTE: subbox must be contained in this
    * FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& plus (const BaseFab<T>& src,
                      const Box&        subbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);
    /**
    * \brief Add srcbox region of src FAB to destbox region of this FAB.
    * The srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& plus (const BaseFab<T>& src,
                      const Box&        srcbox,
                      const Box&        destbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);

    //! Atomic FAB addition (a[i] <- a[i] + b[i]).
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& atomicAdd (const BaseFab<T>& src);

    /**
    * \brief Atomically add src components (srccomp:srccomp+numcomp-1) to
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& atomicAdd (const BaseFab<T>& src,
                           int               srccomp,
                           int               destcomp,
                           int               numcomp=1);
    /**
    * \brief Same as above except addition is restricted to intersection
    * of subbox and src FAB. NOTE: subbox must be contained in this
    * FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& atomicAdd (const BaseFab<T>& src,
                           const Box&        subbox,
                           int               srccomp,
                           int               destcomp,
                           int               numcomp=1);
    /**
    * \brief Atomically add srcbox region of src FAB to destbox region of this FAB.
    * The srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& atomicAdd (const BaseFab<T>& src,
                           const Box&        srcbox,
                           const Box&        destbox,
                           int               srccomp,
                           int               destcomp,
                           int               numcomp=1);

    //! FAB SAXPY (y[i] <- y[i] + a * x[i]), in place.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& saxpy (T a, const BaseFab<T>& x,
                       const Box&        srcbox,
                       const Box&        destbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);
    //! FAB SAXPY (y[i] <- y[i] + a * x[i]), in place.  All components.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& saxpy (T a, const BaseFab<T>& x);

    //! FAB XPAY (y[i] <- x[i] + a * y[i])
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& xpay (T a, const BaseFab<T>& x,
		      const Box&        srcbox,
		      const Box&        destbox,
		      int               srccomp,
		      int               destcomp,
		      int               numcomp=1);

    //! y[i] <- y[i] + x1[i] * x2[i])
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& addproduct (const Box&        destbox,
			    int               destcomp,
			    int               numcomp,
			    const BaseFab<T>& src1,
			    int               comp1,
			    const BaseFab<T>& src2,
			    int               comp2);
    /*
    * \brief Scalar subtraction (a[i] <- a[i] - r).
    * Note: use plus(-r) for more general operations.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& operator-= (T r);

    //! FAB subtraction (a[i] <- a[i] - b[i]), in place.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& operator-= (const BaseFab<T>& f);

    //! FAB subtraction (a[i] <- a[i] - b[i]). The same as -= operator.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& minus (const BaseFab<T>& src);

    /**
    * \brief Subtract src components (srccomp:srccomp+numcomp-1) to
    * this FABs components (destcomp:destcomp+numcomp-1) where
    * the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& minus (const BaseFab<T>& src,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);
    /**
    * \brief Same as above except subtraction is restricted to intersection
    * of subbox and src FAB.  NOTE: subbox must be contained in
    * this FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& minus (const BaseFab<T>& src,
                       const Box&        subbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);
    /**
    * \brief Subtract srcbox region of src FAB from destbox region
    * of this FAB. srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& minus (const BaseFab<T>& src,
                       const Box&        srcbox,
                       const Box&        destbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);
    //! Scalar multiplication (a[i] <- a[i] * r), in place.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& operator*= (T r);

    //! Scalar multiplication (a[i] <- a[i] * r).  The same as *=.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& mult (T r);

    //! Scalar multiplication, except control which components are multiplied.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& mult (T   r,
                      int comp,
                      int numcomp=1);
    /**
    * \brief As above, except specify sub-box.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& mult (T          r,
                      const Box& b,
                      int        comp=0,
                      int        numcomp=1);

    //! FAB multiplication (a[i] <- a[i] * b[i]), in place.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& operator*= (const BaseFab<T>& f);

    //! As above.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& mult (const BaseFab<T>& src);

    /**
    * \brief Multiply src components (srccomp:srccomp+numcomp-1) with
    * this FABs components (destcomp:destcomp+numcomp-1) where
    * the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& mult (const BaseFab<T>& src,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);

    /**
    * \brief Same as above except multiplication is restricted to
    * intersection of subbox and src FAB.  NOTE: subbox must be
    * contained in this FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& mult (const BaseFab<T>& src,
                      const Box&        subbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);

    /**
    * \brief Multiply srcbox region of src FAB with destbox region
    * of this FAB. The srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& mult (const BaseFab<T>& src,
                      const Box&        srcbox,
                      const Box&        destbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);
    //! Scalar division (a[i] <- a[i] / r), in place.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& operator/= (T r);

    //! Scalar division (a[i] <- a[i] / r), in place.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& divide (T r);

    //! As above except specify which components.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& divide (T   r,
                        int comp,
                        int numcomp=1);

    //! As above except specify sub-box.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& divide (T          r,
                        const Box& b,
                        int        comp=0,
                        int        numcomp=1);

    //! FAB division, in place.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& operator/= (const BaseFab<T>& src);

    //! Same as above.
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& divide (const BaseFab<T>& src);

    /**
    * \brief This FAB is numerator, src FAB is denominator
    * divide src components (srccomp:srccomp+numcomp-1) into
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& divide (const BaseFab<T>& src,
                        int               srccomp,
                        int               destcomp,
                        int               numcomp=1);
    /**
    * \brief Same as above except division is restricted to
    * intersection of subbox and src FAB.  NOTE: subbox must be
    * contained in this FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& divide (const BaseFab<T>& src,
                        const Box&        subbox,
                        int               srccomp,
                        int               destcomp,
                        int               numcomp=1);
    /**
    * \brief destbox region of this FAB is numerator. srcbox regions of
    * src FAB is denominator. srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& divide (const BaseFab<T>& src,
                        const Box&        srcbox,
                        const Box&        destbox,
                        int               srccomp,
                        int               destcomp,
                        int               numcomp=1);
    /**
    * \brief Divide wherever "src" is "true" or "non-zero".
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& protected_divide (const BaseFab<T>& src);

    /**
    * \brief Divide wherever "src" is "true" or "non-zero".
    * This FAB is numerator, src FAB is denominator
    * divide src components (srccomp:srccomp+numcomp-1) into
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& protected_divide (const BaseFab<T>& src,
                                  int               srccomp,
                                  int               destcomp,
                                  int               numcomp=1);

    /**
    * \brief Divide wherever "src" is "true" or "non-zero".
    * Same as above except division is restricted to
    * intersection of subbox and src FAB.  NOTE: subbox must be
    * contained in this FAB.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& protected_divide (const BaseFab<T>& src,
                                  const Box&        subbox,
                                  int               srccomp,
                                  int               destcomp,
                                  int               numcomp=1);

    /**
    * Divide wherever "src" is "true" or "non-zero".
    * destbox region of this FAB is numerator. srcbox regions of
    * src FAB is denominator. srcbox and destbox must be same size.
    */
    AMREX_GPU_HOST_DEVICE
    BaseFab<T>& protected_divide (const BaseFab<T>& src,
                                  const Box&        srcbox,
                                  const Box&        destbox,
                                  int               srccomp,
                                  int               destcomp,
                                  int               numcomp=1);

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
    BaseFab<T>& linInterp (const BaseFab<T>& f1,
                           const Box&        b1,
                           int               comp1,
                           const BaseFab<T>& f2,
                           const Box&        b2,
                           int               comp2,
                           Real              t1,
                           Real              t2,
                           Real              t,
                           const Box&        b,
                           int               comp,
                           int               numcomp = 1);

    AMREX_GPU_HOST_DEVICE
    //! Version of linInterp() in which b, b1, & b2 are the same.
    BaseFab<T>& linInterp (const BaseFab<T>& f1,
                           int               comp1,
                           const BaseFab<T>& f2,
                           int               comp2,
                           Real              t1,
                           Real              t2,
                           Real              t,
                           const Box&        b,
                           int               comp,
                           int               numcomp = 1);

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
    BaseFab<T>& linComb (const BaseFab<T>& f1,
                         const Box&        b1,
                         int               comp1,
                         const BaseFab<T>& f2,
                         const Box&        b2,
                         int               comp2,
                         Real              alpha,
                         Real              beta,
                         const Box&        b,
                         int               comp,
                         int               numcomp = 1);

    //! Dot product of x (i.e.,this) and y
    AMREX_GPU_HOST_DEVICE
    T dot (const Box& xbx, int xcomp,
	   const BaseFab<T>& y, const Box& ybx, int ycomp,
	   int numcomp = 1) const;

    AMREX_GPU_HOST_DEVICE
    T dotmask (const BaseFab<int>& mask, const Box& xbx, int xcomp,
               const BaseFab<T>& y, const Box& ybx, int ycomp,
               int numcomp) const;

    //! Change the Box type without change the length
    AMREX_GPU_HOST_DEVICE
    void SetBoxType (const IndexType& typ) { this->domain.setType(typ); }

    template <class F> // AMREX_GPU_HOST_DEVICE
    void ForEach (const Box& b, int c, int nc, F f) {
        ForEachImpl(*this, b, c, nc, f);
    }
    template <class F> // AMREX_GPU_HOST_DEVICE
    void ForEach (const Box& b, int c, int nc, F f) const {
        ForEachImpl(*this, b, c, nc, f);
    }
    template <class F> // AMREX_GPU_HOST_DEVICE
    void ForEachIV (const Box& b, int c, int nc, F f) {
        ForEachIVImpl(*this, b, c, nc, f);
    }
    template <class F> // AMREX_GPU_HOST_DEVICE
    void ForEachIV (const Box& b, int c, int nc, F f) const {
        ForEachIVImpl(*this, b, c, nc, f);
    }
    template <class F> // AMREX_GPU_HOST_DEVICE
    void ForEach (const Box& dstbox, int dstcomp, int numcomp,
                  const BaseFab<T>& src, int srccomp, F f);

    template <class F> // AMREX_GPU_HOST_DEVICE
    void ForEach (const Box& dstbox, int dstcomp, int numcomp,
                  const BaseFab<T>& src, const Box& srcbox, int srccomp, F f);

    template <typename P, class F> // AMREX_GPU_HOST_DEVICE
    P Accumulate (const Box& b, int c, int nc, P init, F f) const;

    template <class F> // AMREX_GPU_HOST_DEVICE
    void Transform (T* dst, const Box& b, int c, int nc, F f) const;

    template <class F> // AMREX_GPU_HOST_DEVICE
    void Transform (const Box& b, int c, int nc, T const* src, F f);

protected:
    //! Allocates memory for the BaseFab<T>.
    void define ();

    //! The function called by BaseFab copy operations.
    AMREX_GPU_HOST_DEVICE
    void performCopy (const BaseFab<T>& src,
                      const Box&        srcbox,
                      int               srccomp,
                      const Box&        destbox,
                      int               destcomp,
                      int               numcomp);

    //! This function is called by the BaseFab setVal operations.
    AMREX_GPU_HOST_DEVICE
    void performSetVal (T          x,
                        const Box& bx,
                        int        nstart,
                        int        numcomp);

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
ForEachImpl (B& fab, const Box& b, int c, int nc, F f)
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
ForEachIVImpl (B& fab, const Box& b, int c, int nc, F f)
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

template <class T>
template <class F>
void
BaseFab<T>::ForEach (const Box& dstbox, int dstcomp, int numcomp,
                     const BaseFab<T>& src, const Box& srcbox, int srccomp, F f)
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

template <class T>
template <class F>
void
BaseFab<T>::ForEach (const Box& dstbox, int dstcomp, int numcomp,
                     const BaseFab<T>& src, int srccomp, F f)
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

template <class T>
template <typename P, class F>
P
BaseFab<T>::Accumulate (const Box& b, int c, int nc, P init, F f) const
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

template <class T>
template <class F>
void
BaseFab<T>::Transform (T* dst, const Box& b, int c, int nc, F f) const
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

template <class T>
template <class F>
void
BaseFab<T>::Transform (const Box& b, int c, int nc, T const* src, F f)
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

template <class T>
inline
T*
BaseFab<T>::dataPtr (const IntVect& p, int n)
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr + (this->domain.index(p)+n*this->numpts);
}

template <class T>
inline
const T*
BaseFab<T>::dataPtr (const IntVect& p, int n) const
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr + (this->domain.index(p)+n*this->numpts);
}

template <class T>
inline
T&
BaseFab<T>::operator() (const IntVect& p,
                        int            n)
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr[this->domain.index(p)+n*this->numpts];
}

template <class T>
inline
T&
BaseFab<T>::operator() (const IntVect& p)
{
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr[this->domain.index(p)];
}

template <class T>
inline
const T&
BaseFab<T>::operator() (const IntVect& p,
                        int            n) const
{
    AMREX_ASSERT(n >= 0);
    AMREX_ASSERT(n < this->nvar);
    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(this->domain.contains(p));

    return this->dptr[this->domain.index(p)+n*this->numpts];
}

template <class T>
inline
const T&
BaseFab<T>::operator() (const IntVect& p) const
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
                     int            numcomp) const
{
    const int loc      = this->domain.index(pos);
    const long sz    = this->domain.numPts();

    AMREX_ASSERT(!(this->dptr == 0));
    AMREX_ASSERT(n >= 0 && n + numcomp <= this->nvar);

    for (int k = 0; k < numcomp; k++)
        data[k] = this->dptr[loc+(n+k)*sz];
}

template <class T>
void
BaseFab<T>::getVal (T*             data,
                    const IntVect& pos) const
{
    getVal(data,pos,0,this->nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::shift (const IntVect& v)
{
    this->domain += v;
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::shift (int idir,
                   int n_cell)
{
    this->domain.shift(idir,n_cell);
    return *this;
}

template <class T>
BaseFab<T> &
BaseFab<T>::shiftHalf (const IntVect& v)
{
    this->domain.shiftHalf(v);
    return *this;
}

template <class T>
BaseFab<T> &
BaseFab<T>::shiftHalf (int idir,
                       int n_cell)
{
    this->domain.shiftHalf(idir,n_cell);
    return *this;
}

template <class T>
void
BaseFab<T>::setVal (T val)
{
    performSetVal(val,box(), 0, this->nvar);
}

template <class T>
void
BaseFab<T>::setVal (T          x,
                    const Box& bx,
                    int        n)
{
    performSetVal(x,bx,n,1);
}

template <class T>
void
BaseFab<T>::setVal (T   x,
                    int n)
{
    performSetVal(x,this->domain,n,1);
}

template <class T>
void
BaseFab<T>::setVal (T          x,
                    const Box& b,
                    int        ns,
                    int        num)
{
    performSetVal(x,b,ns,num);
}

template <class T>
void
BaseFab<T>::setValIfNot (T val, const Box& bx, const BaseFab<int>& mask, int ns, int num)
{
    AMREX_ASSERT(this->domain.contains(bx));
    AMREX_ASSERT(ns >= 0 && ns + num <= this->nvar);

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, ns);
    const auto mp  = mask.view(lo);

    for (int n = 0; n < num; ++n) {
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

template <class T>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src,
                  const Box&        srcbox,
                  int               srccomp,
                  const Box&        destbox,
                  int               destcomp,
                  int               numcomp)
{
    AMREX_ASSERT(destbox.ok());
    AMREX_ASSERT(srcbox.sameSize(destbox));
    AMREX_ASSERT(src.box().contains(srcbox));
    AMREX_ASSERT(this->domain.contains(destbox));
    AMREX_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    AMREX_ASSERT(destcomp >= 0 && destcomp+numcomp <= this->nvar);
    performCopy(src,srcbox,srccomp,destbox,destcomp,numcomp);
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src)
{
    AMREX_ASSERT(this->nvar <= src.nvar);
    AMREX_ASSERT(this->domain.sameType(src.domain));
    Box overlap(this->domain);
    overlap &= src.domain;
    if (overlap.ok())
        performCopy(src,overlap,0,overlap,0,this->nvar);
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src,
                  const Box&        destbox)
{
    AMREX_ASSERT(this->nvar <= src.nvar);
    AMREX_ASSERT(this->domain.contains(destbox));
    Box overlap(destbox);
    overlap &= src.domain;
    if (overlap.ok())
        performCopy(src,overlap,0,overlap,0,this->nvar);
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
{
    AMREX_ASSERT(srccomp >= 0 && srccomp + numcomp <= src.nvar);
    AMREX_ASSERT(destcomp >= 0 && destcomp + numcomp <= this->nvar);
    Box overlap(this->domain);
    overlap &= src.domain;
    if (overlap.ok())
        performCopy(src,overlap,srccomp,overlap,destcomp,numcomp);
    return *this;
}

template <class T>
void
BaseFab<T>::define ()
{
    AMREX_ASSERT(this->nvar > 0);
    AMREX_ASSERT(this->dptr == 0);
    AMREX_ASSERT(this->numpts > 0);
    AMREX_ASSERT(std::numeric_limits<long>::max()/this->nvar > this->numpts);

    this->truesize  = this->nvar*this->numpts;
    this->ptr_owner = true;
    this->dptr = static_cast<T*>(amrex::The_Arena()->alloc(this->truesize*sizeof(T)));
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

    amrex::update_fab_stats(this->numpts, this->truesize, sizeof(T));
}

template <class T>
BaseFab<T>::BaseFab ()
    : BaseFabData<T>{Box(),
                     IntVect(0),
                     0,
                     0L,
                     0L,
                     nullptr,
                     false,
                     false}
{}

template <class T>
BaseFab<T>::BaseFab (const Box& bx,
                     int        n,
		     bool       alloc,
		     bool       shared)
    : BaseFabData<T>{bx,
                     bx.size(),
                     n,
                     bx.numPts(),
                     0L,
                     nullptr,
                     false,
                     shared}
{
    if (!this->shared_memory && alloc) define();
}

template <class T>
BaseFab<T>::BaseFab (const BaseFab<T>& rhs, MakeType make_type, int scomp, int ncomp)
    : BaseFabData<T>{rhs.domain,
                     rhs.dlen,
                     ncomp,
                     rhs.numpts,
                     ncomp*rhs.numpts,
                     const_cast<T*>(rhs.dataPtr(scomp)),
                     false,
                     false}
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

#ifdef AMREX_USE_GPU
template <class T>
BaseFab<T>::BaseFab (const BaseFab<T>& rhs, MakeType make_type)
{
    auto dst = static_cast<BaseFabData<T>*>(this);
    auto src = static_cast<BaseFabData<T>const*>(&rhs);
    if (Gpu::isDevicePtr(&rhs)) {
        AMREX_GPU_SAFE_CALL(cudaMemcpy(dst, src, sizeof(BaseFabData<T>), cudaMemcpyDeviceToHost));
    } else {
        std::memcpy(dst, src, sizeof(BaseFabData<T>));
    }

    if (make_type == amrex::make_deep_copy)
    {
        this->dptr = nullptr;
        define();
        this->copy(rhs, this->domain, 0, this->domain, 0, this->nvar);
    } else if (make_type == amrex::make_alias) {
        this->setOwner(false);
    } else {
        amrex::Abort("BaseFab: unknown MakeType");
    }
}
#endif

template<class T>
BaseFab<T>::BaseFab (const Box& bx, int ncomp, T* p)
    : BaseFabData<T>{bx,
                     bx.size(),
                     ncomp,
                     bx.numPts(),
                     bx.numPts()*ncomp,
                     p,
                     false,
                     false}
{
}

template <class T>
BaseFab<T>::~BaseFab ()
{
    clear();
}

template <class T>
BaseFab<T>::BaseFab (BaseFab<T>&& rhs) noexcept
    : BaseFabData<T>{rhs.domain,
                     rhs.dlen,
                     rhs.nvar,
                     rhs.numpts,
                     rhs.truesize,
                     rhs.dptr,
                     rhs.ptr_owner,
                     rhs.shared_memory}
{
    rhs.dptr = nullptr;
    rhs.ptr_owner = false;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator= (const T& t)
{
    setVal(t);
    return *this;
}

template <class T>
void
BaseFab<T>::resize (const Box& b,
                    int        n)
{
    this->nvar   = n;
    this->domain = b;
    this->dlen   = b.size();
    this->numpts = this->domain.numPts();

    if (this->dptr == 0)
    {
	if (this->shared_memory)
	    amrex::Abort("BaseFab::resize: BaseFab in shared memory cannot increase size");

        define();
    }
    else if (this->nvar*this->numpts > this->truesize)
    {
	if (this->shared_memory)
	    amrex::Abort("BaseFab::resize: BaseFab in shared memory cannot increase size");

        clear();

        define();
    }
}

template <class T>
void
BaseFab<T>::clear ()
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

	    amrex::The_Arena()->free(this->dptr);

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
void
BaseFab<T>::performCopy (const BaseFab<T>& src,
                         const Box&        srcbox,
                         int               srccomp,
                         const Box&        destbox,
                         int               destcomp,
                         int               numcomp)
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
                    dp(i,j,k,n) = sp(i,j,k,n);
                }
            }
        }
    }
}

template <class T>
void
BaseFab<T>::performSetVal (T          val,
                           const Box& bx,
                           int        ns,
                           int        num)
{
    if (!bx.ok()) return;

    AMREX_ASSERT(this->domain.contains(bx));
    AMREX_ASSERT(ns >= 0 && ns + num <= this->nvar);

    const auto len = amrex::length(bx);
    const auto lo  = amrex::lbound(bx);
    const auto dp  = view(lo, ns);

    for (int n = 0; n < num; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) = val;
                }
            }
        }
    }
}

template <class T>
std::size_t
BaseFab<T>::copyToMem (const Box& srcbox,
                       int        srccomp,
                       int        numcomp,
                       void*      dst) const
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

template <class T>
std::size_t
BaseFab<T>::copyFromMem (const Box&  dstbox,
                         int         dstcomp,
                         int         numcomp,
                         const void* src)
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

template <class T>
std::size_t
BaseFab<T>::addFromMem (const Box&  dstbox,
                        int         dstcomp,
                        int         numcomp,
                        const void* src)
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

template <class T>
void
BaseFab<T>::setComplement (T          x,
                           const Box& b,
                           int        ns,
                           int        num)
{
    BoxList b_lst = amrex::boxDiff(this->domain,b);
    for (BoxList::iterator bli = b_lst.begin(), End = b_lst.end(); bli != End; ++bli) {
        performSetVal(x, *bli, ns, num);
    }
}

template <class T>
void
BaseFab<T>::abs ()
{
    this->abs(this->domain,0,this->nvar);
}

template <class T>
void
BaseFab<T>::abs (int comp,
                 int numcomp)
{
    this->abs(this->domain,comp,numcomp);
}

template <class T>
void
BaseFab<T>::abs (const Box& subbox,
                 int        comp,
                 int        numcomp)
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

template <class T>
Real
BaseFab<T>::norminfmask (const Box& subbox, const BaseFab<int>& mask,
                         int scomp, int ncomp) const
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
                    if (mp(i,j,k,n)) {
                        Real t = static_cast<Real>(std::abs(dp(i,j,k,n)));
                        r = (r > t) ? r : t;
                    }
                }
            }
        }
    }

    return r;
}

template <class T>
Real
BaseFab<T>::norm (int p,
                  int comp,
                  int numcomp) const
{
    return norm(this->domain,p,comp,numcomp);
}

template <class T>
Real
BaseFab<T>::norm (const Box& subbox,
                  int        p,
                  int        comp,
                  int        numcomp) const
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
        amrex::Error("BaseFab<T>::norm: wrong p");
    }

    return nrm;
}

template <class T>
T
BaseFab<T>::min (int comp) const
{
    return this->min(this->domain,comp);
}

template <class T>
T
BaseFab<T>::min (const Box& subbox, int comp) const
{
    T r = std::numeric_limits<T>::max();

    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

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

template <class T>
T
BaseFab<T>::max (int comp) const
{
    return this->max(this->domain,comp);
}

template <class T>
T
BaseFab<T>::max (const Box& subbox, int comp) const
{
    T r = std::numeric_limits<T>::lowest();

    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

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

template <class T>
T
BaseFab<T>::maxabs (int comp) const
{
    return this->maxabs(this->domain,comp);
}

template <class T>
T
BaseFab<T>::maxabs (const Box& subbox, int comp) const
{
    T r = std::numeric_limits<T>::lowest();

    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

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


template <class T>
IntVect
BaseFab<T>::minIndex (int comp) const
{
    return this->minIndex(this->domain,comp);
}

template <class T>
IntVect
BaseFab<T>::minIndex (const Box& subbox, int comp) const
{
    AMREX_D_DECL(int imin, jmin, kmin);
    T min_val = std::numeric_limits<T>::max();

    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);
 

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

template <class T>
IntVect
BaseFab<T>::indexFromValue (Real value, const Box& subbox, int comp) const
{
    AMREX_D_DECL(int imin = std::numeric_limits<int>::lowest(),
                     jmin = std::numeric_limits<int>::lowest(),
                     kmin = std::numeric_limits<int>::lowest());

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

template <class T>
void
BaseFab<T>::minIndex (const Box& subbox, Real& min_val, IntVect& min_idx, int comp) const
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

template <class T>
IntVect
BaseFab<T>::maxIndex (int comp) const
{
    return this->maxIndex(this->domain,comp);
}

template <class T>
IntVect
BaseFab<T>::maxIndex (const Box& subbox,
                      int        comp) const
{
    AMREX_D_DECL(int imax, jmax, kmax);
    T max_val = std::numeric_limits<T>::lowest();

    const auto len  = amrex::length(subbox);
    const auto lo   = amrex::lbound(subbox);
    const auto dp   = view(lo, comp);

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

template <class T>
void
BaseFab<T>::maxIndex (const Box& subbox, Real& max_val, IntVect& max_idx, int comp) const
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

template <class T>
int
BaseFab<T>::maskLT (BaseFab<int>& mask,
                    T             val,
                    int           comp) const
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

template <class T>
int
BaseFab<T>::maskLE (BaseFab<int>& mask,
                    T             val,
                    int           comp) const
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

template <class T>
int
BaseFab<T>::maskEQ (BaseFab<int>& mask,
                    T             val,
                    int           comp) const
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

template <class T>
int
BaseFab<T>::maskGT (BaseFab<int>& mask,
                    T             val,
                    int           comp) const
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

template <class T>
int
BaseFab<T>::maskGE (BaseFab<int>& mask,
                    T             val,
                    int           comp) const
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

template <class T>
BaseFab<T>&
BaseFab<T>::plus (T r)
{
    return operator+=(r);
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& x)
{
    return operator+=(x);
}

template <class T>
BaseFab<T>&
BaseFab<T>::atomicAdd (const BaseFab<T>& x)
{
    Box ovlp(this->domain);
    ovlp &= x.domain;
    return ovlp.ok() ? this->atomicAdd(x,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::saxpy (T a, const BaseFab<T>& x,
                   const Box& srcbox, const Box& destbox,
                   int srccomp, int destcomp, int numcomp)
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

template <class T>
BaseFab<T>&
BaseFab<T>::saxpy (T a, const BaseFab<T>& x)
{
    Box ovlp(this->domain);
    ovlp &= x.domain;
    return ovlp.ok() ? saxpy(a,x,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::xpay (T a, const BaseFab<T>& x,
                  const Box& srcbox, const Box& destbox,
                  int srccomp, int destcomp, int numcomp)
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

template <class T>
BaseFab<T>&
BaseFab<T>::addproduct (const Box& destbox, int destcomp, int numcomp,
                        const BaseFab<T>& src1, int comp1,
                        const BaseFab<T>& src2, int comp2)
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

template <class T>
BaseFab<T>&
BaseFab<T>::linComb (const BaseFab<T>& f1, const Box& b1, int comp1,
                     const BaseFab<T>& f2, const Box& b2, int comp2,
                     Real alpha, Real beta, const Box& b,
                     int comp, int numcomp)
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

template <class T>
T
BaseFab<T>::dot (const Box& xbx, int xcomp,
                 const BaseFab<T>& y, const Box& ybx, int ycomp,
                 int numcomp) const
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

template <class T>
T
BaseFab<T>::dotmask (const BaseFab<int>& mask, const Box& xbx, int xcomp,
                     const BaseFab<T>& y, const Box& ybx, int ycomp,
                     int numcomp) const
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
                    int m = static_cast<int>(static_cast<bool>(mp(i,j,k,n)));
                    r += xp(i,j,k,n) * yp(i,j,k,n) * m;
                }
            }
        }
    }

    return r;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator-= (T r)
{
    return operator+=(-r);
}

template <class T>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& x)
{
    return operator-=(x);
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (T r)
{
    return operator*=(r);
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& x)
{
    return operator*=(x);
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (T r)
{
    return operator/=(r);
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& x)
{
    return operator/=(x);
}

template <class T>
T
BaseFab<T>::sum (int comp,
                 int numcomp) const
{
    return this->sum(this->domain,comp,numcomp);
}

template <class T>
T
BaseFab<T>::sum (const Box& subbox,
                 int        comp,
                 int        numcomp) const
{
    BL_ASSERT(this->domain.contains(subbox));
    BL_ASSERT(comp >= 0 && comp + numcomp <= this->nvar);

    T r = 0;

    const auto len = amrex::length(subbox);
    const auto lo  = amrex::lbound(subbox);
    const auto dp  = view(lo, comp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                T const * const AMREX_RESTRICT dpr = &(dp(0,j,k,n));
                for (int i = 0; i < len.x; ++i) {
//                    r += dp(i,j,k,n);
                    r += dpr[i];
                }
            }
        }
    }

    return r;
}

template <class T>
BaseFab<T>&
BaseFab<T>::negate ()
{
    return this->negate(this->domain,0,this->nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::negate (int comp,
                    int numcomp)
{
    return this->negate(this->domain,comp,numcomp);
}

template <class T>
BaseFab<T>&
BaseFab<T>::negate (const Box& b,
                    int        comp,
                    int        numcomp)
{
    BL_ASSERT(this->domain.contains(b));
    BL_ASSERT(comp >= 0 && comp + numcomp <= this->nvar);

    const auto len = amrex::length(b);
    const auto lo  = amrex::lbound(b);
    const auto dp  = view(lo, comp);

    for (int n = 0; n < numcomp; ++n) {
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

template <class T>
BaseFab<T>&
BaseFab<T>::invert (T r)
{
    return this->invert(r,this->domain,0,this->nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::invert (T   r,
                    int comp,
                    int numcomp)
{
    return this->invert(r,this->domain,comp,numcomp);
}

template <class T>
BaseFab<T>&
BaseFab<T>::invert (T          r,
                    const Box& b,
                    int        comp,
                    int        numcomp)
{
    BL_ASSERT(this->domain.contains(b));
    BL_ASSERT(comp >= 0 && comp + numcomp <= this->nvar);

    const auto len = amrex::length(b);
    const auto lo  = amrex::lbound(b);
    const auto dp  = view(lo, comp);

    for (int n = 0; n < numcomp; ++n) {
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

template <class T>
BaseFab<T>&
BaseFab<T>::operator+= (T r)
{
    return this->plus(r,this->domain,0,this->nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (T   r,
                  int comp,
                  int numcomp)
{
    return this->plus(r,this->domain,comp,numcomp);
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (T          r,
                  const Box& b,
                  int        comp,
                  int        numcomp)
{
    AMREX_ASSERT(this->domain.contains(b));
    BL_ASSERT(comp >= 0 && comp + numcomp <= this->nvar);

    const auto len = amrex::length(b);
    const auto lo  = amrex::lbound(b);
    const auto dp  = view(lo, comp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) += r;
                }
            }
        }
    }

    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator+= (const BaseFab<T>& src)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->plus(src,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->plus(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::atomicAdd (const BaseFab<T>& src,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->atomicAdd(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src,
                  const Box&        subbox,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    ovlp &= subbox;
    return ovlp.ok() ? this->plus(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::atomicAdd (const BaseFab<T>& src,
                       const Box&        subbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    ovlp &= subbox;
    return ovlp.ok() ? this->atomicAdd(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src,
                  const Box&        srcbox,
                  const Box&        destbox,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
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

template <class T>
BaseFab<T>&
BaseFab<T>::atomicAdd (const BaseFab<T>& src,
                       const Box&        srcbox,
                       const Box&        destbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp)
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

template <class T>
BaseFab<T>&
BaseFab<T>::operator-= (const BaseFab<T>& src)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->minus(src,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src,
                   int               srccomp,
                   int               destcomp,
                   int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->minus(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src,
                   const Box&        subbox,
                   int               srccomp,
                   int               destcomp,
                   int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    ovlp &= subbox;
    return ovlp.ok() ? this->minus(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src,
                   const Box&        srcbox,
                   const Box&        destbox,
                   int               srccomp,
                   int               destcomp,
                   int               numcomp)
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

template <class T>
BaseFab<T>&
BaseFab<T>::operator*= (T r)
{
    return this->mult(r,this->domain,0,this->nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (T   r,
                  int comp,
                  int numcomp)
{
    return this->mult(r,this->domain,comp,numcomp);
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (T          r,
                  const Box& b,
                  int        comp,
                  int        numcomp)
{
    const auto len = amrex::length(b);
    const auto lo  = amrex::lbound(b);
    const auto dp  = view(lo, comp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD 
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) *= r;
                }
            }
        }
    }

    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator*= (const BaseFab<T>& src)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->mult(src,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->mult(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src,
                  const Box&        subbox,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    ovlp &= subbox;
    return ovlp.ok() ? this->mult(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src,
                  const Box&        srcbox,
                  const Box&        destbox,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
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

template <class T>
BaseFab<T>&
BaseFab<T>::operator/= (T r)
{
    return this->divide(r,this->domain,0,this->nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (T   r,
                    int comp,
                    int numcomp)
{
    return this->divide(r,this->domain,comp,numcomp);
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (T          r,
                    const Box& b,
                    int        comp,
                    int        numcomp)
{
    const auto len = amrex::length(b);
    const auto lo  = amrex::lbound(b);
    const auto dp  = view(lo, comp);

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dp(i,j,k,n) /= r;
                }
            }
        }
    }

    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator/= (const BaseFab<T>& src)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->divide(src,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src,
                    int               srccomp,
                    int               destcomp,
                    int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->divide(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src,
                    const Box&        subbox,
                    int               srccomp,
                    int               destcomp,
                    int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    ovlp &= subbox;
    return ovlp.ok() ? this->divide(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src,
                    const Box&        srcbox,
                    const Box&        destbox,
                    int               srccomp,
                    int               destcomp,
                    int               numcomp)
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

template <class T>
BaseFab<T>&
BaseFab<T>::protected_divide (const BaseFab<T>& src)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->protected_divide(src,ovlp,ovlp,0,0,this->nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::protected_divide (const BaseFab<T>& src,
                              int               srccomp,
                              int               destcomp,
                              int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->protected_divide(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::protected_divide (const BaseFab<T>& src,
                              const Box&        subbox,
                              int               srccomp,
                              int               destcomp,
                              int               numcomp)
{
    Box ovlp(this->domain);
    ovlp &= src.domain;
    ovlp &= subbox;
    return ovlp.ok() ? this->protected_divide(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::protected_divide (const BaseFab<T>& src,
                              const Box&        srcbox,
                              const Box&        destbox,
                              int               srccomp,
                              int               destcomp,
                              int               numcomp)
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

//
// Linear Interpolation / Extrapolation
// Result is (t2-t)/(t2-t1)*f1 + (t-t1)/(t2-t1)*f2
// Data is taken from b1 region of f1, b2 region of f2
// and stored in b region of this FAB.
// Boxes b, b1 and b2 must be the same size.
// Data is taken from component comp1 of f1, comp2 of f2,
// and stored in component comp of this FAB.
// This fab is returned as a reference for chaining.
//

template <class T>
BaseFab<T>&
BaseFab<T>::linInterp (const BaseFab<T>& f1,
                       const Box&        b1,
                       int               comp1,
                       const BaseFab<T>& f2,
                       const Box&        b2,
                       int               comp2,
                       Real              t1,
                       Real              t2,
                       Real              t,
                       const Box&        b,
                       int               comp,
                       int               numcomp)
{
    Real alpha = (t2-t)/(t2-t1);
    Real beta = (t-t1)/(t2-t1);
    return linComb(f1,b1,comp1,f2,b2,comp2,alpha,beta,b,comp,numcomp);
}

template <class T>
BaseFab<T>&
BaseFab<T>::linInterp (const BaseFab<T>& f1,
                       int               comp1,
                       const BaseFab<T>& f2,
                       int               comp2,
                       Real              t1,
                       Real              t2,
                       Real              t,
                       const Box&        b,
                       int               comp,
                       int               numcomp)
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
    copy(f1);
  }
  return *this;
}

}

#endif /*BL_BASEFAB_H*/
