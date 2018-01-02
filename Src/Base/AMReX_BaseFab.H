
#ifndef BL_BASEFAB_H
#define BL_BASEFAB_H

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxList.H>
#include <AMReX_CArena.H>
#include <AMReX_Looping.H>
#include <AMReX_REAL.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_MakeType.H>

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
class BaseFab
{
public:

    typedef T value_type;
    //! Construct an empty BaseFab, which must be resized (see BaseFab::resize) before use.
    BaseFab ();

    //!  Make BaseFab with desired domain (box) and number of components.
    explicit BaseFab (const Box& bx,
                      int        n = 1,
		      bool       alloc = true,
		      bool       shared = false);

    BaseFab (const BaseFab<T>& rhs, MakeType make_type, int scomp, int ncomp);

    /**
     * \brief Create an NON-OWNING BaseFab.  Thus BaseFab is not
     * responsible for memory management.  And it's caller's responsibility that
     * p points to a chunk of memory large enough.
     */
    BaseFab (const Box& bx, int ncomp, T* p);

    //! The destructor deletes the array memory.
    virtual ~BaseFab ();

    BaseFab (const BaseFab<T>& rhs) = delete;
    BaseFab<T>& operator= (const BaseFab<T>& rhs) = delete;
    BaseFab<T>& operator= (BaseFab<T>&& rhs) = delete;

    BaseFab (BaseFab<T>&& rhs) noexcept;

    BaseFab& operator= (const T&);

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

  ///for AggStencil and its minions
  virtual long offset(const IntVect& a_iv, const int& a_ivar) const
    {
      BL_ASSERT(a_ivar >= 0);
      BL_ASSERT(a_ivar < nvar);
      BL_ASSERT(domain.contains(a_iv));

      IntVect ivDiff = a_iv - domain.smallEnd();
#if BL_SPACEDIM > 1
      IntVect ivSize = domain.size();
#endif
      long retval = ivDiff[0];

#if BL_SPACEDIM > 1
      retval += ivDiff[1]*ivSize[0] ;
#endif
#if BL_SPACEDIM > 2
      retval += ivDiff[2]*ivSize[0]*ivSize[1];
#endif
      retval += numpts * a_ivar;

      return retval;
    }

    /**
     * \brief The function returns the BaseFab to the invalid state.
     * The memory is freed.
     */
    void clear ();

    //! Returns how many bytes used
    std::size_t nBytes () const { return truesize*sizeof(T); }

    static bool preAllocatable () { return true; }

    static bool isCopyOMPSafe () { return true; }

    //! Returns bytes used in the Box for those components
    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const
        { return bx.numPts() * sizeof(T) * ncomps; }

    //! Returns the number of components
    int nComp () const { return nvar; }

    //! for calls to fortran.  
    const int* nCompPtr() const {
#ifdef AMREX_USE_CUDA
        return nvar_d;
#else
        return &nvar;
#endif
    }

    //! Returns the number of points 
    long nPts () const { return numpts; }

    //! Returns the total number of points of all components
    long size () const { return nvar*numpts; }

    //! Returns the domain (box) where the array is defined
    const Box& box () const { return domain; }

    /**
    * \brief Returns a pointer to an array of SPACEDIM integers
    * giving the length of the domain in each direction
    */
    const int* length () const { return dlen.getVect(); }

    /**
    * \brief Returns the lower corner of the domain
    * See class Box for analogue.
    */
    const IntVect& smallEnd () const { return domain.smallEnd(); }

    //!  Returns the upper corner of the domain.  See class Box for analogue.
    const IntVect& bigEnd () const { return domain.bigEnd(); }

    /**
    * \brief Returns the lower corner of the domain.  
  
    *Instead of returning them in the form of INTVECTs, as in smallEnd and
    * bigEnd, it returns the values as a pointer to an array of
    * constant integers.  This is useful when interfacing to
    * Fortran subroutines.
    */
    const int* loVect () const { return domain.loVect(); }

    /**
    * \brief Returns the upper corner of the domain.  

    *Instead of returning them in the form of INTVECTs, as in smallEnd and
    * bigEnd, it returns the values as a pointer to an array of
    * constant integers.  This is useful when interfacing to
    * Fortran subroutines.
    */
    const int* hiVect () const { return domain.hiVect(); }

    /**
    * \brief Returns true if the domain of fab is totally contained within
    * the domain of this BaseFab.
    */
    bool contains (const BaseFab<T>& fab) const
    {
        return box().contains(fab.box()) && nvar <= fab.nvar;
    }

    /**
    * \brief Returns true if bx is totally contained
    * within the domain of this BaseFab.
    */
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
    T* dataPtr (int n = 0) { BL_ASSERT(!(dptr == 0)); return &dptr[n*numpts]; }

    //! Same as above except works on const FABs.
    const T* dataPtr (int n = 0) const { BL_ASSERT(!(dptr == 0)); return &dptr[n*numpts]; }

    void setPtr (T* p, long sz) { BL_ASSERT(dptr == 0 && truesize == 0); dptr = p; truesize = sz; }

    //! Returns true if the data for the FAB has been allocated.
    bool isAllocated () const { return dptr != 0; }

    /**
    * \brief Returns a reference to the Nth component value
    * defined at position p in the domain.  This operator may be
    * inefficient if the C++ compiler is unable to optimize the
    * C++ code.
    */
    T& operator() (const IntVect& p, int N);

    //! Same as above, except returns component 0.
    T& operator() (const IntVect& p);

    //! Same as above except works on const FABs.
    const T& operator() (const IntVect& p, int N) const;

    //! Same as above, except returns component 0.
    const T& operator() (const IntVect& p) const;

    /**
    * \brief This function puts numcomp component values, starting at
    * component N, from position pos in the domain into array data,
    * that must be allocated by the user.
    */
    void getVal (T*             data,
                 const IntVect& pos,
                 int            N,
                 int            numcomp) const;
    //! Same as above, except that starts at component 0 and copies all comps.
    void getVal (T*             data,
                 const IntVect& pos) const;
    /**
    * \brief The setVal functions set sub-regions in the BaseFab to a
    * constant value.  This most general form specifies the sub-box,
    * the starting component number, and the number of components
    * to be set.
    */
    void setVal (T          x,
                 const Box& bx,
                 int        nstart,
                 int        ncomp);
    //! Same as above, except the number of modified components is one. N is the component to be modified.
    void setVal (T          x,
                 const Box& bx,
                 int        N = 0);
    //! Same as above, except the sub-box defaults to the entire domain.
    void setVal (T   x,
                 int N);
    //! Same as above, except all components are set.
    void setVal (T x);

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
    BaseFab<T>& copy (const BaseFab<T>& src,
                      const Box&        srcbox,
                      int               srccomp,
                      const Box&        destbox,
                      int               destcomp,
                      int               numcomp);

  //for debugging
  BaseFab<T>&
  slowCopy(const BaseFab<T>& src,
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
    BaseFab<T>& copy (const BaseFab<T>& src,
                      const Box&        destbox);

    /**
    * As above, except that the destbox defaults to the entire domain
    * of the destination BaseFab, and all components are copied.
    * Note that the actual copy is made by the function
    * performCopy() of this class.
    */
    BaseFab<T>& copy (const BaseFab<T>& src);

    //! Copy from the srcbox of this Fab to raw memory and return the number of bytes copied
    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const;
    //! Copy from raw memory to the dstbox of this Fab and return the number of bytes copied
    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);

    /**
    * \brief Perform shifts upon the domain of the BaseFab. They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    BaseFab<T>& shift (const IntVect& v);
    /**
    * \brief Perform shifts upon the domain of the BaseFab.  They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    BaseFab<T>& shift (int idir,
                       int n_cell);
    /**
    * \brief Perform shifts upon the domain of the BaseFab.  They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    BaseFab<T>& shiftHalf (int dir,
                           int num_halfs);
    /**
    * \brief Perform shifts upon the domain of the BaseFab. They are
    * completely analogous to the corresponding Box functions.
    * There is no effect upon the array memory.
    */
    BaseFab<T>& shiftHalf (const IntVect& num_halfs);

    Real norminfmask (const Box& subbox, const BaseFab<int>& mask, int scomp=0, int ncomp=1) const;

    /**
    * \brief Compute the Lp-norm of this FAB using components (scomp : scomp+ncomp-1).
    *   p < 0  -> ERROR
    *   p = 0  -> infinity norm (max norm)
    *   p = 1  -> sum of ABS(FAB)
    */
    Real norm (int p,
               int scomp = 0,
               int ncomp = 1) const;
    
    //! Same as above except only on given subbox.
    Real norm (const Box& subbox,
               int        p,
               int        scomp = 0,
               int        ncomp = 1) const;
    //!Compute absolute value for all components of this FAB.
    void abs ();
    //! Same as above except only for components (comp: comp+numcomp-1)
    void abs (int comp,
              int numcomp=1);
    /**
    * \brief Calculate abs() on subbox for given component range.
    */
    void abs (const Box& subbox,
              int        comp = 0,
              int        numcomp=1);
    /**
    * \return Minimum value of given component.
    */
    T min (int comp = 0) const;
    /**
    * \return Minimum value of given component in given subbox.
    */
    T min (const Box& subbox,
           int        comp = 0) const;
    /**
    * \return Maximum value of given component.
    */
    T max (int comp = 0) const;
    /**
    * \return Maximum value of given component in given subbox.
    */
    T max (const Box& subbox,
           int        comp = 0) const;
    /**
    * \return location of minimum value in given component.
    */
    IntVect minIndex (int comp = 0) const;
    /**
    * \return location of minimum value in given component in
    * given subbox.
    */
    IntVect minIndex (const Box& subbox,
                      int        comp = 0) const;
    /**
    * \return location of maximum value in given component.
    */
    IntVect maxIndex (int comp = 0) const;
    /**
    * \return location of maximum value in given component in given
    * subbox.
    */
    IntVect maxIndex (const Box& subbox,
                      int        comp = 0) const;
    /**
    * \brief Compute mask array with value of 1 in cells where
    * BaseFab has value less than val, 0 otherwise.
    * mask is resized by this function.
    * The number of cells marked with 1 returned.
    */
    int maskLT (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
    //! Same as above except mark cells with value less than or equal to val.
    int maskLE (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
    
    //! Same as above except mark cells with value equal to val.
    int maskEQ (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
    //! Same as above except mark cells with value greater than val.
    int maskGT (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
    //! Same as above except mark cells with value greater than or equal to val.
    int maskGE (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
    //! Fill with a pattern of numbers.
    void patternFill (int mark = 0);
    /**
    * \brief Copies with index reversal from srcbox region of src
    * into destbox region of this FAB.  All components are
    * copied and multiplied by corresponding multiplier.
    */
    void copyRev (const Box&        destbox,
                  const BaseFab<T>& src,
                  const Box&        srcbox,
                  int               reversal_index,
                  T*                multiplier);
    //! Returns sum of given component of FAB state vector.
    T sum (int comp,
           int numcomp = 1) const;
    //! Compute sum of given component of FAB state vector in given subbox.
    T sum (const Box& subbox,
           int        comp,
           int        numcomp = 1) const;
    //! Most general version, specify subbox and which components.
    BaseFab<T>& invert (T          v,
                        const Box& subbox,
                        int        comp=0,
                        int        numcomp=1);
    //! As above except on entire domain.
    BaseFab<T>& invert (T   v,
                        int comp,
                        int numcomp=1);
    //! As above except on entire domain, all components.
    BaseFab<T>& invert (T v);

    //! Negate BaseFab, most general.
    BaseFab<T>& negate (const Box& subbox,
                        int        comp=0,
                        int        numcomp=1);
    //! As above, except on entire domain.
    BaseFab<T>& negate (int comp,
                        int numcomp=1);
    //! As above, except on entire domain and all components.
    BaseFab<T>& negate ();

    //! Scalar addition (a[i] <- a[i] + r), most general.
    BaseFab<T>& plus (T          r,
                      const Box& b,
                      int        comp=0,
                      int        numcomp=1);

    //! As above, except on entire domain.
    BaseFab<T>& plus (T   r,
                      int comp,
                      int numcomp=1);
    //! As above, except on entire domain and all components.
    BaseFab<T>& plus (T r);

    //! Addition in place.  This will often be more efficient than  making new BaseFab for result.
    BaseFab<T>& operator+= (T r);

    //! FAB addition (a[i] <- a[i] + b[i]) in place.
    BaseFab<T>& operator+= (const BaseFab<T>& f);

    //! FAB addition (a[i] <- a[i] + b[i]). The same as += operator.
    BaseFab<T>& plus (const BaseFab<T>& src);

    /**
    * \brief Add src components (srccomp:srccomp+numcomp-1) to
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
    BaseFab<T>& plus (const BaseFab<T>& src,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);
    /**
    * \brief Same as above except addition is restricted to intersection
    * of subbox and src FAB. NOTE: subbox must be contained in this
    * FAB.
    */
    BaseFab<T>& plus (const BaseFab<T>& src,
                      const Box&        subbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);
    /**
    * \brief Add srcbox region of src FAB to destbox region of this FAB.
    * The srcbox and destbox must be same size.
    */
    BaseFab<T>& plus (const BaseFab<T>& src,
                      const Box&        srcbox,
                      const Box&        destbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);
    //! FAB SAXPY (y[i] <- y[i] + a * x[i]), in place.
    BaseFab<T>& saxpy (T a, const BaseFab<T>& x,
                       const Box&        srcbox,
                       const Box&        destbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);
    //! FAB SAXPY (y[i] <- y[i] + a * x[i]), in place.  All components.
    BaseFab<T>& saxpy (T a, const BaseFab<T>& x);

    //! FAB XPAY (y[i] <- x[i] + a * y[i])
    BaseFab<T>& xpay (T a, const BaseFab<T>& x,
		      const Box&        srcbox,
		      const Box&        destbox,
		      int               srccomp,
		      int               destcomp,
		      int               numcomp=1);

    //! y[i] <- y[i] + x1[i] * x2[i])
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
    BaseFab<T>& operator-= (T r);

    //! FAB subtraction (a[i] <- a[i] - b[i]), in place.
    BaseFab<T>& operator-= (const BaseFab<T>& f);

    //! FAB subtraction (a[i] <- a[i] - b[i]). The same as -= operator.
    BaseFab<T>& minus (const BaseFab<T>& src);

    /**
    * \brief Subtract src components (srccomp:srccomp+numcomp-1) to
    * this FABs components (destcomp:destcomp+numcomp-1) where
    * the two FABs intersect.
    */
    BaseFab<T>& minus (const BaseFab<T>& src,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);
    /**
    * \brief Same as above except subtraction is restricted to intersection
    * of subbox and src FAB.  NOTE: subbox must be contained in
    * this FAB.
    */
    BaseFab<T>& minus (const BaseFab<T>& src,
                       const Box&        subbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);
    /**
    * \brief Subtract srcbox region of src FAB from destbox region
    * of this FAB. srcbox and destbox must be same size.
    */
    BaseFab<T>& minus (const BaseFab<T>& src,
                       const Box&        srcbox,
                       const Box&        destbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);
    //! Scalar multiplication (a[i] <- a[i] * r), in place.
    BaseFab<T>& operator*= (T r);

    //! Scalar multiplication (a[i] <- a[i] * r).  The same as *=.
    BaseFab<T>& mult (T r);

    //! Scalar multiplication, except control which components are multiplied.
    BaseFab<T>& mult (T   r,
                      int comp,
                      int numcomp=1);
    /**
    * \brief As above, except specify sub-box.
    */
    BaseFab<T>& mult (T          r,
                      const Box& b,
                      int        comp=0,
                      int        numcomp=1);

    //! FAB multiplication (a[i] <- a[i] * b[i]), in place.
    BaseFab<T>& operator*= (const BaseFab<T>& f);

    //! As above.
    BaseFab<T>& mult (const BaseFab<T>& src);

    /**
    * \brief Multiply src components (srccomp:srccomp+numcomp-1) with
    * this FABs components (destcomp:destcomp+numcomp-1) where
    * the two FABs intersect.
    */
    BaseFab<T>& mult (const BaseFab<T>& src,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);

    /**
    * \brief Same as above except multiplication is restricted to
    * intersection of subbox and src FAB.  NOTE: subbox must be
    * contained in this FAB.
    */
    BaseFab<T>& mult (const BaseFab<T>& src,
                      const Box&        subbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);

    /**
    * \brief Multiply srcbox region of src FAB with destbox region
    * of this FAB. The srcbox and destbox must be same size.
    */
    BaseFab<T>& mult (const BaseFab<T>& src,
                      const Box&        srcbox,
                      const Box&        destbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);
    //! Scalar division (a[i] <- a[i] / r), in place.
    BaseFab<T>& operator/= (T r);

    //! Scalar division (a[i] <- a[i] / r), in place.
    BaseFab<T>& divide (T r);

    //! As above except specify which components.
    BaseFab<T>& divide (T   r,
                        int comp,
                        int numcomp=1);

    //! As above except specify sub-box.
    BaseFab<T>& divide (T          r,
                        const Box& b,
                        int        comp=0,
                        int        numcomp=1);

    //! FAB division, in place.
    BaseFab<T>& operator/= (const BaseFab<T>& src);

    //! Same as above.
    BaseFab<T>& divide (const BaseFab<T>& src);

    /**
    * \brief This FAB is numerator, src FAB is denominator
    * divide src components (srccomp:srccomp+numcomp-1) into
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
    BaseFab<T>& divide (const BaseFab<T>& src,
                        int               srccomp,
                        int               destcomp,
                        int               numcomp=1);
    /**
    * \brief Same as above except division is restricted to
    * intersection of subbox and src FAB.  NOTE: subbox must be
    * contained in this FAB.
    */
    BaseFab<T>& divide (const BaseFab<T>& src,
                        const Box&        subbox,
                        int               srccomp,
                        int               destcomp,
                        int               numcomp=1);
    /**
    * \brief destbox region of this FAB is numerator. srcbox regions of
    * src FAB is denominator. srcbox and destbox must be same size.
    */
    BaseFab<T>& divide (const BaseFab<T>& src,
                        const Box&        srcbox,
                        const Box&        destbox,
                        int               srccomp,
                        int               destcomp,
                        int               numcomp=1);
    /**
    * \brief Divide wherever "src" is "true" or "non-zero".
    */
    BaseFab<T>& protected_divide (const BaseFab<T>& src);

    /**
    * \brief Divide wherever "src" is "true" or "non-zero".
    * This FAB is numerator, src FAB is denominator
    * divide src components (srccomp:srccomp+numcomp-1) into
    * this FABs components (destcomp:destcomp+numcomp-1)
    * where the two FABs intersect.
    */
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
    T dot (const Box& xbx, int xcomp,
	   const BaseFab<T>& y, const Box& ybx, int ycomp,
	   int numcomp = 1) const;

    //! Change the Box type without change the length
    void SetBoxType(const IndexType& typ) { domain.setType(typ); }
    
protected:
    //! Allocates memory for the BaseFab<T>.
    void define ();

    //! The function called by BaseFab copy operations.
    void performCopy (const BaseFab<T>& src,
                      const Box&        srcbox,
                      int               srccomp,
                      const Box&        destbox,
                      int               destcomp,
                      int               numcomp);

    //! This function is called by the BaseFab setVal operations.
    void performSetVal (T          x,
                        const Box& bx,
                        int        nstart,
                        int        numcomp);
protected:

    Box     domain;   // My index space.
    IntVect dlen;     // Length of domain in each direction.
    int     nvar;     // Number components.
    long    numpts;   // Cached number of points in FAB.
    long    truesize; // nvar*numpts that was allocated on heap.
    T*      dptr;     // The data pointer.
    bool    ptr_owner;// Owner of T*?
    bool    shared_memory;  // Is the memory allocated in shared memory?

#ifdef AMREX_USE_CUDA
    int* nvar_d;
#endif
};

template <class T>
inline
T&
BaseFab<T>::operator() (const IntVect& p,
                        int            n)
{
    BL_ASSERT(n >= 0);
    BL_ASSERT(n < nvar);
    BL_ASSERT(!(dptr == 0));
    BL_ASSERT(domain.contains(p));

    return dptr[domain.index(p)+n*numpts];
}

template <class T>
inline
T&
BaseFab<T>::operator() (const IntVect& p)
{
    BL_ASSERT(!(dptr == 0));
    BL_ASSERT(domain.contains(p));

    return dptr[domain.index(p)];
}

template <class T>
inline
const T&
BaseFab<T>::operator() (const IntVect& p,
                        int            n) const
{
    BL_ASSERT(n >= 0);
    BL_ASSERT(n < nvar);
    BL_ASSERT(!(dptr == 0));
    BL_ASSERT(domain.contains(p));

    return dptr[domain.index(p)+n*numpts];
}

template <class T>
inline
const T&
BaseFab<T>::operator() (const IntVect& p) const
{
    BL_ASSERT(!(dptr == 0));
    BL_ASSERT(domain.contains(p));

    return dptr[domain.index(p)];
}

template <class T>
void
BaseFab<T>::getVal  (T*             data,
                     const IntVect& pos,
                     int            n,
                     int            numcomp) const
{
    const int loc      = domain.index(pos);
    const long sz    = domain.numPts();

    BL_ASSERT(!(dptr == 0));
    BL_ASSERT(n >= 0 && n + numcomp <= nvar);

    for (int k = 0; k < numcomp; k++)
        data[k] = dptr[loc+(n+k)*sz];
}

template <class T>
void
BaseFab<T>::getVal (T*             data,
                    const IntVect& pos) const
{
    getVal(data,pos,0,nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::shift (const IntVect& v)
{
    domain += v;
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::shift (int idir,
                   int n_cell)
{
    domain.shift(idir,n_cell);
    return *this;
}

template <class T>
BaseFab<T> &
BaseFab<T>::shiftHalf (const IntVect& v)
{
    domain.shiftHalf(v);
    return *this;
}

template <class T>
BaseFab<T> &
BaseFab<T>::shiftHalf (int idir,
                       int n_cell)
{
    domain.shiftHalf(idir,n_cell);
    return *this;
}

template <class T>
void
BaseFab<T>::setVal (T val)
{
    performSetVal(val,box(), 0, nvar);
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
    performSetVal(x,domain,n,1);
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
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src,
                  const Box&        srcbox,
                  int               srccomp,
                  const Box&        destbox,
                  int               destcomp,
                  int               numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(srcbox.sameSize(destbox));
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(domain.contains(destbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nvar);
    performCopy(src,srcbox,srccomp,destbox,destcomp,numcomp);
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::slowCopy(const BaseFab<T>& src,
                     const Box&        srcbox,
                     int               srccomp,
                     const Box&        destbox,
                     int               destcomp,
                     int               numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(srcbox.sameSize(destbox));
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(domain.contains(destbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nvar);
    for(BoxIterator bit(destbox); bit.ok(); ++bit)
    {
      for(int icomp = 0; icomp < numcomp; icomp++)
      {
        int isrc =  srccomp  + icomp;
        int idst =  destcomp + icomp;
        (*this)(bit(), idst) = src(bit(), isrc);
      }
    }
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src)
{
    BL_ASSERT(nvar <= src.nvar);
    BL_ASSERT(domain.sameType(src.domain));
    Box overlap(domain);
    overlap &= src.domain;
    if (overlap.ok())
        performCopy(src,overlap,0,overlap,0,nvar);
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src,
                  const Box&        destbox)
{
    BL_ASSERT(nvar <= src.nvar);
    BL_ASSERT(domain.contains(destbox));
    Box overlap(destbox);
    overlap &= src.domain;
    if (overlap.ok())
        performCopy(src,overlap,0,overlap,0,nvar);
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
{
    BL_ASSERT(srccomp >= 0 && srccomp + numcomp <= src.nvar);
    BL_ASSERT(destcomp >= 0 && destcomp + numcomp <= nvar);
    Box overlap(domain);
    overlap &= src.domain;
    if (overlap.ok())
        performCopy(src,overlap,srccomp,overlap,destcomp,numcomp);
    return *this;
}

template <class T>
void
BaseFab<T>::define ()
{
    BL_ASSERT(nvar > 0);
    BL_ASSERT(dptr == 0);
    BL_ASSERT(numpts > 0);
    BL_ASSERT(std::numeric_limits<long>::max()/nvar > numpts);

    truesize  = nvar*numpts;
    dptr      = static_cast<T*>(amrex::The_Arena()->alloc(truesize*sizeof(T)));
    ptr_owner = true;
    //
    // Now call T::T() on the raw memory so we have valid Ts.
    //
    T* ptr = dptr;
    //
    // Note this must be long not int for very large (e.g.,1024^3) boxes.
    //
    for (long i = 0; i < truesize; i++, ptr++)
    {
        new (ptr) T;
    }

#ifdef AMREX_USE_CUDA
    nvar_d = static_cast<int*>(amrex::The_Nvar_Arena()->alloc(sizeof(int)));
    *nvar_d = nvar;
#endif

    amrex::update_fab_stats(numpts, truesize, sizeof(T));
}

template <class T>
BaseFab<T>::BaseFab ()
    :
    domain(Box()),
    nvar(0),
    numpts(0),
    truesize(0),
    dptr(0),
    ptr_owner(false),
    shared_memory(false)
{}

template <class T>
BaseFab<T>::BaseFab (const Box& bx,
                     int        n,
		     bool       alloc,
		     bool       shared)
    :
    domain(bx),
    dlen(bx.size()),
    nvar(n),
    numpts(bx.numPts()),
    truesize(0),
    dptr(0),
    ptr_owner(false),
    shared_memory(shared)
{
    if (!shared_memory && alloc) define();
}

template <class T>
BaseFab<T>::BaseFab (const BaseFab<T>& rhs, MakeType make_type, int scomp, int ncomp)
    :
    domain(rhs.domain),
    dlen(rhs.dlen),
    nvar(ncomp),
    numpts(rhs.numpts),
    truesize(ncomp*rhs.numpts),
    dptr(const_cast<T*>(rhs.dataPtr(scomp))),
    ptr_owner(false),
    shared_memory(false)
{
    BL_ASSERT(scomp+ncomp <= rhs.nComp());
    if (make_type == amrex::make_deep_copy)
    {
        dptr = nullptr;
        define();
        this->copy(rhs, domain, scomp, domain, 0, ncomp);
    } else if (make_type == amrex::make_alias) {
        ; // nothing to do
    } else {
        amrex::Abort("BaseFab: unknown MakeType");
    }
}

template<class T>
BaseFab<T>::BaseFab (const Box& bx, int ncomp, T* p)
    :
    domain(bx),
    dlen(bx.size()),
    nvar(ncomp),
    numpts(bx.numPts()),
    truesize(bx.numPts()*ncomp),
    dptr(p),
    ptr_owner(false),
    shared_memory(false)
{
}

template <class T>
BaseFab<T>::~BaseFab ()
{
    clear();
}

template <class T>
BaseFab<T>::BaseFab (BaseFab<T>&& rhs) noexcept
    :
    domain(rhs.domain),
    dlen(rhs.dlen),
    nvar(rhs.nvar),
    numpts(rhs.numpts),
    truesize(rhs.truesize),
    dptr(rhs.dptr),
    ptr_owner(rhs.ptr_owner),
    shared_memory(rhs.shared_memory)
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
    nvar   = n;
    domain = b;
    dlen   = b.size();
    numpts = domain.numPts();

    if (dptr == 0)
    {
	if (shared_memory)
	    amrex::Abort("BaseFab::resize: BaseFab in shared memory cannot increase size");

        define();
    }
    else if (nvar*numpts > truesize)
    {
	if (shared_memory)
	    amrex::Abort("BaseFab::resize: BaseFab in shared memory cannot increase size");

        clear();

        define();
    }
}

template <class T>
void
BaseFab<T>::clear ()
{
    if (dptr)
    {
        //
        // Call T::~T() on the to-be-destroyed memory.
        //
	if (ptr_owner)
	{
	    if (shared_memory)
	    {
		amrex::Abort("BaseFab::clear: BaseFab cannot be owner of shared memory");
	    }

	    T* ptr = dptr;

	    for (long i = 0; i < truesize; i++, ptr++)
	    {
		ptr->~T();
	    }
	    amrex::The_Arena()->free(dptr);

#ifdef AMREX_USE_CUDA
            amrex::The_Nvar_Arena()->free(nvar_d);
#endif

	    if (nvar > 1) {
		amrex::update_fab_stats(-truesize/nvar, -truesize, sizeof(T));
	    } else {
		amrex::update_fab_stats(0, -truesize, sizeof(T));
	    }
	}

	dptr = 0;
	truesize = 0;
    }
}

//
// performCopy() has been rewritten here so we can insert pragmas which
// will enhance vectorization on the Crays.  The downside is that the
// code is greatly expanded and rather incomprehensible.
//

template <class T>
void
BaseFab<T>::performCopy (const BaseFab<T>& src,
                         const Box&        srcbox,
                         int               srccomp,
                         const Box&        destbox,
                         int               destcomp,
                         int               numcomp)
{
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

#if (BL_SPACEDIM == 1)
{                                                                       
    BL_ASSERT((destcomp) >= 0 && (destcomp) + (numcomp) <= nComp());
    BL_ASSERT((srccomp) >= 0 && (srccomp) + (numcomp) <= (src).nComp());
    Box _subbox_ = box(); 
    _subbox_ &= destbox; 
    BL_ASSERT(srcbox.sameSize(_subbox_)); 
    if (_subbox_.ok()) 
    { 
        const int *_th_plo = loVect(); 
        const int *_th_plen = length(); 
        const int *_x_plo = (src).loVect(); 
        const int *_x_plen = (src).length(); 
        const int *_subbox_lo = _subbox_.loVect();
        IntVect _subbox_length = _subbox_.size();
        const int *_subbox_len = _subbox_length.getVect(); 
        const int *_bx_lo = (srcbox).loVect(); 
        T* _th_p = dataPtr(destcomp); 
        const T* _x_p  = (src).dataPtr(srccomp); 
        for (int _n = 0; _n < (numcomp); ++_n)
        { 
            T *_th_pp = _th_p + ((_subbox_lo[0]-_th_plo[0])+_n*_th_plen[0]); 
            const T *_x_pp = _x_p + ((_bx_lo[0]-_x_plo[0])+_n*_x_plen[0]);
            for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)
            { 
                int iR = _i + _subbox_lo[0]; iR += 0; 
                int ixR = _i + _bx_lo[0]; ixR += 0; 
                T &thisR = * _th_pp; const T & srcR = _x_pp[_i];
#elif (BL_SPACEDIM == 2)
{                                                                       
    BL_ASSERT((destcomp) >= 0 && (destcomp) + (numcomp) <= nComp()); 
    BL_ASSERT((srccomp) >= 0 && (srccomp) + (numcomp) <= (src).nComp()); 
    Box _subbox_ = box(); 
    _subbox_ &= destbox; 
    BL_ASSERT(srcbox.sameSize(_subbox_)); 
    if (_subbox_.ok())
    { 
        const int *_th_plo = loVect(); 
        const int *_th_plen = length(); 
        const int *_x_plo = (src).loVect(); 
        const int *_x_plen = (src).length(); 
        const int *_subbox_lo = _subbox_.loVect();
        IntVect _subbox_length = _subbox_.size();
        const int *_subbox_len = _subbox_length.getVect(); 
        const int *_bx_lo = (srcbox).loVect(); 
        T* _th_p = dataPtr(destcomp); 
        const T* _x_p  = (src).dataPtr(srccomp); 
        for (int _n = 0; _n < (numcomp); ++_n)
        { 
            for(int _j = 0; _j < _subbox_len[1]; ++_j)
            { 
                const int jR = _j + _subbox_lo[1]; 
                const int jxR = _j + _bx_lo[1]; 
                T *_th_pp = _th_p + ((_subbox_lo[0] - _th_plo[0]) 
                                     + _th_plen[0]*((jR - _th_plo[1]) 
                                                    + _n * _th_plen[1])); 
                const T *_x_pp = _x_p + ((_bx_lo[0] - _x_plo[0]) 
                                         + _x_plen[0]*((jxR - _x_plo[1]) 
                                                       + _n * _x_plen[1])); 
                for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)
                {  
                    T &thisR = * _th_pp; const T & srcR = _x_pp[_i]; 
#elif (BL_SPACEDIM == 3)
{ 
    BL_ASSERT((destcomp) >= 0 && (destcomp) + (numcomp) <= nComp()); 
    BL_ASSERT((srccomp) >= 0 && (srccomp) + (numcomp) <= (src).nComp()); 
    Box _subbox_(box()); 
    _subbox_ &= destbox; 
    BL_ASSERT((srcbox).sameSize(_subbox_)); 
    if (_subbox_.ok())
    { 
        const int *_th_plo = loVect(); 
        const int *_th_plen = length(); 
        const int *_x_plo = (src).loVect(); 
        const int *_x_plen = (src).length(); 
        const int *_subbox_lo = _subbox_.loVect();
        IntVect _subbox_length = _subbox_.size();
        const int *_subbox_len = _subbox_length.getVect(); 
        const int *_bx_lo = (srcbox).loVect(); 
        T* _th_p = dataPtr(destcomp); 
        const T* _x_p  = (src).dataPtr(srccomp); 
        for (int _n = 0; _n < (numcomp); ++_n)
        { 
            for (int _k = 0; _k < _subbox_len[2]; ++_k)
            { 
                const int kR = _k + _subbox_lo[2]; 
                const int kxR = _k + _bx_lo[2]; 
                for(int _j = 0; _j < _subbox_len[1]; ++_j)
                { 
                    const int jR = _j + _subbox_lo[1]; 
                    const int jxR = _j + _bx_lo[1]; 
                    T *_th_pp = _th_p + ((_subbox_lo[0] - _th_plo[0]) 
                                         + _th_plen[0]*((jR - _th_plo[1]) 
                                                        + _th_plen[1]*(
                                                            (kR - _th_plo[2]) 
                                                        + _n * _th_plen[2]))); 
                    const T *_x_pp = _x_p + ((_bx_lo[0] - _x_plo[0]) 
                                             + _x_plen[0]*((jxR - _x_plo[1]) 
                                                           + _x_plen[1]*( 
                                                            (kxR - _x_plo[2])
                                                         + _n * _x_plen[2])));
                    for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)
                    {
                        T &thisR = * _th_pp; const T & srcR = _x_pp[_i]; 
#endif
    {
        thisR = srcR;
    }
#if (BL_SPACEDIM == 1)
     }}}}
#elif (BL_SPACEDIM == 2)
     }}}}}
#elif (BL_SPACEDIM == 3)
     }}}}}}
#endif
}

template <class T>
void
BaseFab<T>::performSetVal (T          val,
                           const Box& bx,
                           int        ns,
                           int        num)
{
    BL_ASSERT(domain.contains(bx));
    BL_ASSERT(ns >= 0 && ns + num <= nvar);

    if (bx == domain)
    {
        T* data = &dptr[ns*numpts];

        for (long i = 0, N = num*numpts; i < N; i++)
        {
            *data++ = val;
        }
    }
    else
    {
        ForAllThisBNN(T,bx,ns,num)
        {
            thisR = val;
        } EndFor
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
        const int* _x_plo         = loVect(); 
        const int* _x_plen        = length(); 
        IntVect    _subbox_length = srcbox.size();
        const int* _subbox_len    = _subbox_length.getVect(); 
        const int* _bx_lo         = (srcbox).loVect(); 
        T*         _th_p          = static_cast<T*>(dst);
        const T*   _x_p           = dataPtr(srccomp); 

#if (BL_SPACEDIM == 1)

        for (int _n = 0; _n < (numcomp); ++_n)
        { 
            const T* _x_pp = _x_p + ((_bx_lo[0]-_x_plo[0])+_n*_x_plen[0]);

            for (int _i = 0; _i < _subbox_len[0]; ++_i)
            { 
                *_th_p++ = *_x_pp++;
            }
        }

#elif (BL_SPACEDIM == 2)

        for (int _n = 0; _n < (numcomp); ++_n)
        { 
            for (int _j = 0; _j < _subbox_len[1]; ++_j)
            { 
                const int jxR   = _j + _bx_lo[1]; 
                const T*  _x_pp = _x_p + ((_bx_lo[0] - _x_plo[0]) 
                                          + _x_plen[0]*((jxR - _x_plo[1]) 
                                                        + _n * _x_plen[1])); 
                for (int _i = 0; _i < _subbox_len[0]; ++_i)
                {  
                    *_th_p++ = *_x_pp++;
                }
            }
        }

#elif (BL_SPACEDIM == 3)

        for (int _n = 0; _n < (numcomp); ++_n)
        { 
            for (int _k = 0; _k < _subbox_len[2]; ++_k)
            { 
                const int kxR = _k + _bx_lo[2];
                for(int _j = 0; _j < _subbox_len[1]; ++_j)
                { 
                    const int jxR   = _j + _bx_lo[1]; 
                    const T*  _x_pp = _x_p + ((_bx_lo[0] - _x_plo[0]) 
                                              + _x_plen[0]*((jxR - _x_plo[1]) 
                                                           + _x_plen[1]*( 
                                                               (kxR - _x_plo[2])
                                                            + _n * _x_plen[2])));
                    for (int _i = 0; _i < _subbox_len[0]; ++_i)
                    {
                        *_th_p++ = *_x_pp++;
                    }
                }
            }
        }
#endif
        return sizeof(T)*numcomp*srcbox.numPts();
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
        const int* _x_plo         = loVect(); 
        const int* _x_plen        = length(); 
        IntVect    _subbox_length = dstbox.size();
        const int* _subbox_len    = _subbox_length.getVect(); 
        const int* _bx_lo         = (dstbox).loVect(); 
        const T*   _th_p          = static_cast<const T*>(src);
        T*         _x_p           = dataPtr(dstcomp); 

#if (BL_SPACEDIM == 1)

        for (int _n = 0; _n < (numcomp); ++_n)
        { 
            T* _x_pp = _x_p + ((_bx_lo[0]-_x_plo[0])+_n*_x_plen[0]);

            for (int _i = 0; _i < _subbox_len[0]; ++_i)
            { 
                *_x_pp++ = *_th_p++;
            }
        }

#elif (BL_SPACEDIM == 2)

        for (int _n = 0; _n < (numcomp); ++_n)
        { 
            for (int _j = 0; _j < _subbox_len[1]; ++_j)
            { 
                const int jxR   = _j + _bx_lo[1]; 
                T*        _x_pp = _x_p + ((_bx_lo[0] - _x_plo[0]) + _x_plen[0]*((jxR - _x_plo[1]) + _n * _x_plen[1])); 
                for (int _i = 0; _i < _subbox_len[0]; ++_i)
                {  
                    *_x_pp++ = *_th_p++;
                }
            }
        }

#elif (BL_SPACEDIM == 3)

        for (int _n = 0; _n < (numcomp); ++_n)
        { 
            for (int _k = 0; _k < _subbox_len[2]; ++_k)
            { 
                const int kxR = _k + _bx_lo[2];
                for(int _j = 0; _j < _subbox_len[1]; ++_j)
                { 
                    const int jxR   = _j + _bx_lo[1]; 
                    T*        _x_pp = _x_p + ((_bx_lo[0] - _x_plo[0]) 
                                              + _x_plen[0]*((jxR - _x_plo[1]) 
                                                            + _x_plen[1]*( (kxR - _x_plo[2]) + _n * _x_plen[2])));
                    for (int _i = 0; _i < _subbox_len[0]; ++_i)
                    {
                        *_x_pp++ = *_th_p++;
                    }
                }
            }
        }
#endif
        return sizeof(T)*numcomp*dstbox.numPts();
    }
    else
    {
        return 0;
    }
}

#if !defined(BL_NO_FORT)
//
// Forward declaration of template specializatons for Real.
// Definitions are found in BaseFab.cpp.
//
template <>
BaseFab<Real>&
BaseFab<Real>::invert (Real       val,
                       const Box& subbox,
                       int        comp,
                       int        numcomp);
template <>
void
BaseFab<Real>::performCopy (const BaseFab<Real>& src,
                            const Box&           srcbox,
                            int                  srccomp,
                            const Box&           destbox,
                            int                  destcomp,
                            int                  numcomp);
template <>
std::size_t
BaseFab<Real>::copyToMem (const Box& srcbox,
                          int        srccomp,
                          int        numcomp,
                          void*      dst) const;
template <>
std::size_t
BaseFab<Real>::copyFromMem (const Box&  dstbox,
                            int         dstcomp,
                            int         numcomp,
                            const void* src);
template <>
void
BaseFab<Real>::performSetVal (Real       val,
                              const Box& bx,
                              int        ns,
                              int        num);
template <>
Real
BaseFab<Real>::norminfmask (const Box& subbox, const BaseFab<int>& mask, int comp, int ncomp) const;

template <>
Real
BaseFab<Real>::norm (const Box& subbox,
                     int        p,
                     int        comp,
                     int        ncomp) const;

template <>
Real
BaseFab<Real>::sum (const Box& subbox,
                    int        comp,
                    int        ncomp) const;
template <>
BaseFab<Real>&
BaseFab<Real>::plus (const BaseFab<Real>& src,
                     const Box&           srcbox,
                     const Box&           destbox,
                     int                  srccomp,
                     int                  destcomp,
                     int                  numcomp);
template <>
BaseFab<Real>&
BaseFab<Real>::mult (const BaseFab<Real>& src,
                     const Box&           srcbox,
                     const Box&           destbox,
                     int                  srccomp,
                     int                  destcomp,
                     int                  numcomp);
template <>
BaseFab<Real>&
BaseFab<Real>::saxpy (Real a, const BaseFab<Real>& x,
                      const Box&        srcbox,
                      const Box&        destbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp);

template <>
BaseFab<Real>&
BaseFab<Real>::minus (const BaseFab<Real>& src,
                      const Box&           srcbox,
                      const Box&           destbox,
                      int                  srccomp,
                      int                  destcomp,
                      int                  numcomp);
template <>
BaseFab<Real>&
BaseFab<Real>::divide (const BaseFab<Real>& src,
                       const Box&           srcbox,
                       const Box&           destbox,
                       int                  srccomp,
                       int                  destcomp,
                       int                  numcomp);
template <>
BaseFab<Real>&
BaseFab<Real>::protected_divide (const BaseFab<Real>& src,
                                 const Box&           srcbox,
                                 const Box&           destbox,
                                 int                  srccomp,
                                 int                  destcomp,
                                 int                  numcomp);

template <>
BaseFab<Real>&
BaseFab<Real>::linComb (const BaseFab<Real>&  f1,
			const Box&         b1,
			int                comp1,
			const BaseFab<Real>&  f2,
			const Box&         b2,
			int                comp2,
			Real               alpha,
			Real               beta,
			const Box&         b,
			int                comp,
			int                numcomp);

#endif

template <class T>
void
BaseFab<T>::setComplement (T          x,
                           const Box& b,
                           int        ns,
                           int        num)
{
    BoxList b_lst = amrex::boxDiff(domain,b);
    for (BoxList::iterator bli = b_lst.begin(), End = b_lst.end(); bli != End; ++bli)
        performSetVal(x, *bli, ns, num);
}

template <class T>
void
BaseFab<T>::abs ()
{
    return this->abs(domain,0,nvar);
}

template <class T>
void
BaseFab<T>::abs (int comp,
                 int numcomp)
{
    return this->abs(domain,comp,numcomp);
}

template <class T>
void
BaseFab<T>::abs (const Box& subbox,
                 int        comp,
                 int        numcomp)
{
    ForAllThisBNN(Real,subbox,comp,numcomp)
    {
        thisR = std::abs(thisR);
    } EndFor
}

template <class T>
Real
BaseFab<T>::norm (int p,
                  int comp,
                  int numcomp) const
{
    return norm(domain,p,comp,numcomp);
}

template <class T>
Real
BaseFab<T>::norm (const Box& subbox,
                  int        p,
                  int        comp,
                  int        numcomp) const
{
    BL_ASSERT(comp >= 0 && comp+numcomp <= nComp());
    BL_ASSERT(p >= 0);

    Real  nrm    = 0;
    Real* tmp    = 0;
    int   tmplen = 0;

    if (p == 0)
    {
        ForAllThisCPencil(T,subbox,comp,numcomp)
        {
            const T* row = &thisR;
            if (tmp == 0)
            {
                tmp    = new Real[thisLen];
                tmplen = thisLen;
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = std::abs(row[i]);
            }
            else
            {
                for (int i = 0; i < thisLen; i++) {
		    // MLW: explicit cast to Real for SX6
		    Real rw = std::abs(row[i]);
                    tmp[i] = std::max(tmp[i],rw);
		}
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
            nrm = std::max(nrm, tmp[i]);
    }
    else if (p == 1)
    {
        ForAllThisCPencil(T,subbox,comp,numcomp)
        {
            const T* row = &thisR;
            if (tmp == 0)
            {
                tmp    = new Real[thisLen];
                tmplen = thisLen;
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = std::abs(row[i]);
            }
            else
            {
                for (int i = 0; i < thisLen; i++)
                    tmp[i] += std::abs(row[i]);
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
            nrm += tmp[i];
    }
    else
    {
      amrex::Error("BaseFab::norm(): only p == 0 or p == 1 are supported");
    }

    delete [] tmp;

    return nrm;
}

template <class T>
T
BaseFab<T>::min (int comp) const
{
    return this->min(domain,comp);
}

template <class T>
T
BaseFab<T>::min (const Box& subbox,
                 int        comp) const
{
    T *_min_row = 0;
    int _x_len = 0;
    ForAllThisCPencil(T,subbox,comp,1)
    {
        const T* _row = &thisR;
        if (_min_row == 0)
        {
            _min_row = new T[thisLen];
            _x_len = thisLen;
            for (int i = 0; i < thisLen; i++)
                _min_row[i] = _row[i];
        }
        else
        {
            for (int i = 0; i < thisLen; i++)
                _min_row[i] = std::min(_row[i],_min_row[i]);
        }
    } EndForPencil;

    T _min = _min_row[0];
    for (int i = 1; i < _x_len; i++)
        _min = std::min(_min,_min_row[i]);

    delete [] _min_row;

    return _min;
}

template <class T>
T
BaseFab<T>::max (int comp) const
{
    return this->max(domain,comp);
}

template <class T>
T
BaseFab<T>::max (const Box& subbox,
                 int        comp) const
{
    T*  _max_row = 0;
    int _x_len   = 0;
    ForAllThisCPencil(T,subbox,comp,1)
    {
        const T* _row = &thisR;
        if (_max_row == 0)
        {
            _max_row = new T[thisLen];
            _x_len = thisLen;
            for (int i = 0; i < thisLen; i++)
                _max_row[i] = _row[i];
        }
        else
        {
            for (int i = 0; i < thisLen; i++)
                _max_row[i] = std::max(_row[i],_max_row[i]);
        }
    } EndForPencil;

    T _max = _max_row[0];
    for (int i = 1; i < _x_len; i++)
        _max = std::max(_max,_max_row[i]);

    delete [] _max_row;

    return _max;
}

template <class T>
IntVect
BaseFab<T>::minIndex (int comp) const
{
    return this->minIndex(domain,comp);
}

template <class T>
IntVect
BaseFab<T>::minIndex (const Box& subbox,
                      int        comp) const
{
    IntVect _min_loc(subbox.smallEnd());
    T _min_val = (*this).operator()(_min_loc,comp);
    ForAllThisCBNN(T,subbox,comp,1)
    {
        if (thisR < _min_val)
        {
            _min_val = thisR;
            AMREX_D_EXPR(_min_loc[0] = iR,
                   _min_loc[1] = jR,
                   _min_loc[2] = kR);
        }
    } EndFor;

    return _min_loc;
}

template <class T>
IntVect
BaseFab<T>::maxIndex (int comp) const
{
    return this->maxIndex(domain,comp);
}

template <class T>
IntVect
BaseFab<T>::maxIndex (const Box& subbox,
                      int        comp) const
{
    IntVect _max_loc(subbox.smallEnd());
    T _max_val = (*this).operator()(_max_loc,comp);
    ForAllThisCBNN(T,subbox,comp,1)
    {
        if (thisR > _max_val)
        {
            _max_val = thisR;
            AMREX_D_EXPR(_max_loc[0] = iR,
                   _max_loc[1] = jR,
                   _max_loc[2] = kR);
        }
    } EndFor;

    return _max_loc;
}

template <class T>
int
BaseFab<T>::maskLT (BaseFab<int>& mask,
                    T             val,
                    int           comp) const
{
    mask.resize(domain,1);
    mask.setVal(0);

    int* mptr = mask.dataPtr();
    int  cnt  = 0;

    ForAllThisCBNN(T,domain,comp,1)
    {
        int ix = AMREX_D_TERM(_i, +_j*_b_len[0], +_k*_b_len[0]*_b_len[1]);
        if (thisR < val)
        {
            mptr[ix] = 1;
            cnt++;
        }
    } EndFor;

    return cnt;
}

template <class T>
int
BaseFab<T>::maskLE (BaseFab<int>& mask,
                    T             val,
                    int           comp) const
{
    mask.resize(domain,1);
    mask.setVal(0);

    int* mptr = mask.dataPtr();
    int  cnt  = 0;

    ForAllThisCBNN(T,domain,comp,1)
    {
        int ix = AMREX_D_TERM(_i, +_j*_b_len[0], +_k*_b_len[0]*_b_len[1]);
        if (thisR <= val)
        {
            mptr[ix] = 1;
            cnt++;
        }
    } EndFor;

    return cnt;
}

template <class T>
int
BaseFab<T>::maskEQ (BaseFab<int>& mask,
                    T             val,
                    int           comp) const
{
    mask.resize(domain,1);
    mask.setVal(0);

    int* mptr = mask.dataPtr();
    int  cnt  = 0;

    ForAllThisCBNN(T,domain,comp,1)
    {
        int ix = AMREX_D_TERM(_i, +_j*_b_len[0], +_k*_b_len[0]*_b_len[1]);
        if (thisR == val)
        {
            mptr[ix] = 1;
            cnt++;
        }
    } EndFor;

    return cnt;
}

template <class T>
int
BaseFab<T>::maskGT (BaseFab<int>& mask,
                    T             val,
                    int           comp) const
{
    mask.resize(domain,1);
    mask.setVal(0);

    int* mptr = mask.dataPtr();
    int  cnt  = 0;

    ForAllThisCBNN(T,domain,comp,1)
    {
        int ix = AMREX_D_TERM(_i, +_j*_b_len[0], +_k*_b_len[0]*_b_len[1]);
        if (thisR > val)
        {
            mptr[ix] = 1;
            cnt++;
        }
    } EndFor;

    return cnt;
}

template <class T>
int
BaseFab<T>::maskGE(BaseFab<int>& mask,
                   T             val,
                   int           comp) const
{
    mask.resize(domain,1);
    mask.setVal(0);

    int* mptr = mask.dataPtr();
    int  cnt  = 0;

    ForAllThisCBNN(T,domain,comp,1)
    {
        int ix = AMREX_D_TERM(_i, +_j*_b_len[0], +_k*_b_len[0]*_b_len[1]);
        if (thisR >= val)
        {
            mptr[ix] = 1;
            cnt++;
        }
    } EndFor;

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
BaseFab<T>::saxpy (T a, const BaseFab<T>& x)
{
    Box ovlp(domain);
    ovlp &= x.domain;
    return ovlp.ok() ? saxpy(a,x,ovlp,ovlp,0,0,nvar) : *this;
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
void
BaseFab<T>::patternFill (int mark)
{
    ForAllThis(T)
    {
        thisR = AMREX_D_TERM(iR*100, +jR*10, + kR) + 1000*nR + 10000*mark;
    } EndFor
}

template <class T>
void
BaseFab<T>::copyRev (const Box&        destbox,
                     const BaseFab<T>& src,
                     const Box&        srcbox,
                     int               reversal_index,
                     T*                multiplier)
{
    BaseFab<T>& dest = *this;

    ForAllRevXBNYCBNNN(T,dest,destbox,0,src,srcbox,0,nComp(),reversal_index)
    {
        destR = multiplier[_n]*srcR;
    } EndFor
}

template <class T>
T
BaseFab<T>::sum (int comp,
                 int numcomp) const
{
    return this->sum(domain,comp,numcomp);
}

template <class T>
T
BaseFab<T>::sum (const Box& subbox,
                 int        comp,
                 int        numcomp) const
{
    T*  _sum_row = 0;
    int _sum_len = 0;
    ForAllThisCPencil(T,subbox,comp,numcomp)
    {
        const T* _row = &thisR;
        if (_sum_row == 0)
        {
            _sum_row = new T[thisLen];
            _sum_len = thisLen;
            for (int i = 0; i < thisLen; i++)
                _sum_row[i] = _row[i];
        }
        else
        {
            for (int i = 0; i < thisLen; i++)
            {
                _sum_row[i] += _row[i];
            }
        }
    } EndForPencil;

    T _sum = _sum_row[0];
    for (int i = 1; i < _sum_len; i++)
        _sum += _sum_row[i];

    delete [] _sum_row;

    return _sum;
}

template <class T>
BaseFab<T>&
BaseFab<T>::negate ()
{
    return this->negate(domain,0,nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::negate (int comp,
                    int numcomp)
{
    return this->negate(domain,comp,numcomp);
}

template <class T>
BaseFab<T>&
BaseFab<T>::negate (const Box& b,
                    int        comp,
                    int        numcomp)
{
    ForAllThisBNN(T,b,comp,numcomp)
    {
        thisR = - thisR;
    } EndFor
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::invert (T r)
{
    return this->invert(r,domain,0,nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::invert (T   r,
                    int comp,
                    int numcomp)
{
    return this->invert(r,domain,comp,numcomp);
}

template <class T>
BaseFab<T>&
BaseFab<T>::invert (T          r,
                    const Box& b,
                    int        comp,
                    int        numcomp)
{
    ForAllThisBNN(T,b,comp,numcomp)
    {
        thisR = r/thisR;
    } EndFor
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator+= (T r)
{
    return this->plus(r,domain,0,nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (T   r,
                  int comp,
                  int numcomp)
{
    return this->plus(r,domain,comp,numcomp);
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (T          r,
                  const Box& b,
                  int        comp,
                  int        numcomp)
{
    ForAllThisBNN(T,b,comp,numcomp)
    {
        thisR += r;
    } EndFor
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator+= (const BaseFab<T>& src)
{
    Box ovlp(domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->plus(src,ovlp,ovlp,0,0,nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
{
    Box ovlp(domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->plus(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src,
                  const Box&        subbox,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
{
    Box ovlp(domain);
    ovlp &= src.domain;
    ovlp &= subbox;
    return ovlp.ok() ? this->plus(src,ovlp,ovlp,srccomp,destcomp,numcomp) : *this;
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
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        thisR += srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator-= (const BaseFab<T>& src)
{
    Box ovlp(domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->minus(src,ovlp,ovlp,0,0,nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src,
                   int               srccomp,
                   int               destcomp,
                   int               numcomp)
{
    Box ovlp(domain);
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
    Box ovlp(domain);
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
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        thisR -= srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator*= (T r)
{
    return this->mult(r,domain,0,nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (T   r,
                  int comp,
                  int numcomp)
{
    return this->mult(r,domain,comp,numcomp);
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (T          r,
                  const Box& b,
                  int        comp,
                  int        numcomp)
{
    ForAllThisBNN(T,b,comp,numcomp)
    {
        thisR *= r;
    } EndFor
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator*= (const BaseFab<T>& src)
{
    Box ovlp(domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->mult(src,ovlp,ovlp,0,0,nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
{
    Box ovlp(domain);
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
    Box ovlp(domain);
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
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        thisR *= srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator/= (T r)
{
    return this->divide(r,domain,0,nvar);
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (T   r,
                    int comp,
                    int numcomp)
{
    return this->divide(r,domain,comp,numcomp);
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (T          r,
                    const Box& b,
                    int        comp,
                    int        numcomp)
{
    ForAllThisBNN(T,b,comp,numcomp)
    {
        thisR /= r;
    } EndFor
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator/= (const BaseFab<T>& src)
{
    Box ovlp(domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->divide(src,ovlp,ovlp,0,0,nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src,
                    int               srccomp,
                    int               destcomp,
                    int               numcomp)
{
    Box ovlp(domain);
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
    Box ovlp(domain);
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
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        thisR /= srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::protected_divide (const BaseFab<T>& src)
{
    Box ovlp(domain);
    ovlp &= src.domain;
    return ovlp.ok() ? this->protected_divide(src,ovlp,ovlp,0,0,nvar) : *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::protected_divide (const BaseFab<T>& src,
                              int               srccomp,
                              int               destcomp,
                              int               numcomp)
{
    Box ovlp(domain);
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
    Box ovlp(domain);
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
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        if (srcR)
            thisR /= srcR;
    } EndForTX
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

//
// Linear combination, Result is alpha*f1 + beta*f2
// Data is taken from b1 region of f1, b2 region of f2
// and stored in b region of this FAB.
// Boxes b, b1 and b2 must be the same size.
// Data is taken from component comp1 of f1, comp2 of f2,
// and stored in component comp of this FAB.
// This fab is returned as a reference for chaining.
//

template <class T>
BaseFab<T>&
BaseFab<T>::linComb (const BaseFab<T>&  f1,
                     const Box&         b1,
                     int                comp1,
                     const BaseFab<T>&  f2,
                     const Box&         b2,
                     int                comp2,
                     Real               alpha,
                     Real               beta,
                     const Box&         b,
                     int                comp,
                     int                numcomp)
{
    ForAllThisBNNXCBNYCBN(T,b,comp,numcomp,f1,b1,comp1,f2,b2,comp2)
    {
        thisR = (T) (alpha*Real(f1R) + beta*Real(f2R));
    } EndForTX
    return *this;
}

    // \cond CODEGEN
    class BF_init
    {
    public:
        BF_init ();
        ~BF_init ();
    private:
        static int m_cnt;
    };
    // \endcond
}

static amrex::BF_init file_scope_BF_init_object;

#endif /*BL_BASEFAB_H*/
