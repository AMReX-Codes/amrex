//BL_COPYRIGHT_NOTICE

#ifndef BL_BASEFAB_H
#define BL_BASEFAB_H

//
// $Id: BaseFab.H,v 1.21 1999-04-29 20:41:45 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include <Assert.H>
#include <Box.H>
#include <BoxList.H>
#include <CArena.H>
#include <Looping.H>
#include <PArray.H>

#ifdef BL_USE_POINTLIB
#ifndef BL_CRAY_BUG_DEFARG
#include <PointBaseFab.H>
#include <PointDomain.H>
#endif /*BL_CRAY_BUG_DEFARG*/
#endif /*BL_USE_POINTLIB*/

#ifdef BL_USE_POINTLIB
#include <PointLooping.H>
#endif

//
//@Man:
//@Memo: A Fortran Array-like Object
/*@Doc:

  BaseFab emulates the Fortran array concept.  
  Useful operations can be performed upon
  BaseFab's in C++, and they provide a convenient interface to
  Fortran when it is necessary to retreat into that language.
        
  `BaseFab' is a template class.  Through use of the
  template, a `BaseFab' may be based upon any class.  So far at least,
  most applications have been based upon simple types like `integer''s,
  `real*4''s, or `real*8''s.  Most applications do not use BaseFab's 
  directly, but utilize specialized classes derived from BaseFab.
        
  It will be easier to use a `BaseFab' if you understand the following
  concepts.  `BaseFab' objects depend on the dimensionality of space
  (indirectly through the DOMAIN `Box' member).  It is
  typical to define the macro `SPACEDIM' to be 1, 2, or 3 to indicate
  the dimension of space.  See the discussion of class `Box' for more
  information.  A `BaseFab' contains a `Box' DOMAIN, which indicates the
  integer indexing space over which the array is defined.  A `BaseFab'
  also has NVAR components.  By components, we mean that for each
  point in the rectangular indexing space, there are NVAR values
  associated with that point.  A Fortran array corresponding to a
  `BaseFab' would have (SPACEDIM+1) dimensions.
        
  By design, the array layout in a `BaseFab' mirrors that of a
  Fortran array.  The first index (x direction for example) varies
  most rapidly, the next index (y direction), if any, varies next
  fastest. The component index varies last, after all the spatial
  indices.
        
  It is sometimes convenient to be able to treat a sub-array within an
  existing `BaseFab' as a `BaseFab' in its own right.  This is often
  referred to as 'aliasing' the `BaseFab'.  Note that when aliasing is 
  used, the BaseFab's domain will not, in general, be the same as the 
  parent BaseFab's domain, nor will the number of components.
  BaseFab is a dimension dependent class, so SPACEDIM must be 
  defined as either 1, 2, or 3 when compiling.

  This is NOT a polymorphic class.

  It does NOT provide a copy constructor or assignment operator.

  T MUST have a default constructor and an assignment operator.
*/

template <class T>
class BaseFab
{
public:

    /*@ManDoc: Construct an invalid `BaseFab'.  The domain is invalid, the
               number of components is zero, and no actually array memory is
               allocated. An invalid `BaseFab' must be resized
               (see `BaseFab::resize') before use.
    */
    BaseFab ();
    //
    //@ManDoc: Make BaseFab with desired domain and number of components.
    //
    explicit BaseFab (const Box& bx,
                      int        n = 1);

#ifdef BL_USE_POINTLIB
#ifndef BL_CRAY_BUG_DEFARG

    /*@ManDoc: This constructor constructs a BaseFab based on the points in a
               PointBaseFab f. This constructor is only defined
               if BL\_CRAY\_BUG\_DEFARG not defined and if BL\_USE\_POINTLIB
               is defined.
    */
    BaseFab (const PointBaseFab<PointDomain,T>& f,
             T                                  v);
#endif
#endif
    //
    //@ManDoc: The destructor deletes the array memory.
    //
    ~BaseFab ();

    /*@ManDoc: This function resizes a `BaseFab' so it covers the `Box' B
               with N components.  The default action
               is that under resize()ing, the memory allocated for the
               `BaseFab' only grows and never shrinks.  This function is
               particularly useful when a `BaseFab' is used as a temporary
               space which must be a different size whenever it is used.
               Resize()ing a temp will often be faster than re-allocating a
               `BaseFab' because memory allocation can often be avoided.
    */
    void resize (const Box& b,
                 int        N = 1);

    /*@ManDoc: The function returns the `BaseFab' to the invalid state.  (See
               comments for constructors above.)  The memory is freed.
    */
    void clear ();
    //
    //@ManDoc: Returns the number of components.
    //
    int nComp () const;
    //
    //@ManDoc: Returns the domain (box) where the array is defined.
    //
    const Box& box () const;

    /*@ManDoc: Returns a pointer to an array of SPACEDIM integers
               giving the length of the domain in each direction.
    */
    const int* length () const;

    /*@ManDoc: Returns the lower corner of the domain.  See
               class `Box' for analogue.
    */
    const IntVect& smallEnd () const;
 
    /*@ManDoc:  Returns the upper corner of the domain.  See
                class `Box' for analogue.
    */
    const IntVect& bigEnd () const;

    /*@ManDoc: Returns the lower corner of the domain.  Instead of
               returning them in the form of INTVECT's, as in smallEnd and
               bigEnd, it returns the values as a pointer to an array of
               constant integers.  This is useful when interfacing to
               Fortran subroutines.
    */
    const int* loVect () const;
 
    /*@ManDoc: Returns the upper corner of the domain.  Instead of
               returning them in the form of INTVECT's, as in smallEnd and
               bigEnd, it returns the values as a pointer to an array of
               constant integers.  This is useful when interfacing to
               Fortran subroutines.
    */
    const int* hiVect () const;
 
    /*@ManDoc: Returns true if the domain of fab is totally contained within
               the domain of this `BaseFab'.
    */
    bool contains (const BaseFab<T>& fab) const;
 
    /*@ManDoc: Returns true if bx is totally contained
               within the domain of this `BaseFab'.
    */
    bool contains (const Box& bx) const;
 
    /*@ManDoc: Returns a pointer to an object of type T that is the
               value of the Nth component associated with the cell at the
               low end of the domain.  This is commonly used to get a pointer
               to data in the array which is then handed off to a Fortran
               subroutine.  Remember that data is stored in Fortran array
               order, with the component index coming last.   In other words,
               `dataPtr' returns a pointer to all the Nth components.
    */
    T* dataPtr (int N = 0);
    //
    //@ManDoc: Same as above except works on const FAB's.
    //
    const T* dataPtr (int N = 0) const;

    /*@ManDoc: Returns a reference to the Nth component value
               defined at position p in the domain.  This operator may be
               inefficient if the C++ compiler is unable to optimize the
               C++ code.
    */
    T& operator() (const IntVect& p,
                   int            N);
    //
    //@ManDoc: Same as above, except returns component 0.
    //
    T& operator() (const IntVect& p);
    //
    //@ManDoc: Same as above except works on const FAB's.
    //
     const T& operator() (const IntVect& p,
                          int            N) const;
    //
    //@ManDoc: Same as above, except returns component 0.
    //
    const T& operator() (const IntVect& p) const;

    /*@ManDoc: This function puts numcomp component values, starting at
               component N, from position pos in the domain into array data,
               that must be allocated by the user.
    */
    void getVal (T*             data,
                 const IntVect& pos,
                 int            N,
                 int            numcomp) const;

    /*@ManDoc: Same as above, except that starts at component 0 and
               copies ALL comps.
    */
    void getVal (T*             data,
                 const IntVect& pos) const;

    /*@ManDoc: The setVal functions set subregions in the `BaseFab' to a
               constant value.  This most general form specifies the sub-box,
               the starting component number, and the number of components
               to be set.
    */
    void setVal (T          x,
                 const Box& bx,
                 int        nstart,
                 int        ncomp);

    /*@ManDoc: Same as above, except the number of modified components is one.
               N is the component to be modified.
    */
    void setVal (T          x,
                 const Box& bx,
                 int        N);
    //
    //@ManDoc: Same as above, except the subbox defaults to the entire domain.
    //
    void setVal (T   x,
                 int N);
    //
    //@ManDoc: Same as above, except ALL components are set.
    //
    void setVal (T x);

    /*@ManDoc: This function is analogous to the fourth form of
               setVal above, except that instead of setting values on the
               `Box' b, values are set on the complement of b in the domain.
    */
    void setComplement (T          x,
                        const Box& b,
                        int        ns,
                        int        num);

    /*@ManDoc: The `copy' functions copy the contents of one `BaseFab' into
               another.  The destination `BaseFab' is always the object which
               invokes the function.  This, the most general form of copy,
               specifies the contents of any sub-box srcbox in `BaseFab' src
               may be copied into a (possibly different) destbox in the
               destination `BaseFab'.  Note that although the srcbox and the
               destbox may be disjoint, they must be the same size and shape.
               If the sizes differ, the copy is undefined and a runtime error
               results.  This copy function is the only one of the copy
               functions to allow a copy between differing boxes. The user
               also specifies how many components are copied, starting at
               component srccomp in src and stored starting at component
               destcomp. Note that the actual copy is made by the
               function `performCopy' of this class.  The results are
               UNDEFINED if the src and dest are the same and the srcbox and
               destbox overlap.
    */
    BaseFab<T>& copy (const BaseFab<T>& src,
                      const Box&        srcbox,
                      int               srccomp,
                      const Box&        destbox,
                      int               destcomp,
                      int               numcomp);

    /*@ManDoc: As above, except the destination `Box' and the source `Box'
               are taken to be the entire domain of the destination.   A copy
               of the intersecting region is performed.  Note that the actual
               copy is made by the function `performCopy' of this
               class.
    */
    BaseFab<T>& copy (const BaseFab<T>& src,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp = 1);

    /*@ManDoc: As above, except that the destination `Box' is specified,
               but the source `Box' is taken to the equal to the source
               `Box', and all components of the destination `BaseFab' are
               copied.  Note that the actual copy is made by the
               function `performCopy' of this class.
    */
    BaseFab<T>& copy (const BaseFab<T>& src,
                      const Box&        destbox);

    /*@ManDoc: As above, except that the destbox defaults to the entire domain
               of the destination BaseFab, and all components are copied.
               Note that the actual copy is made by the function
               `performCopy' of this class.
    */
    BaseFab<T>& copy (const BaseFab<T>& src);

    /*@ManDoc: Perform shifts upon the domain of the `BaseFab'. They are
               completely analogous to the corresponding `Box' functions.
               There is no effect upon the array memory.
    */
    BaseFab<T>& shift (const IntVect& v);

    /*@ManDoc: Perform shifts upon the domain of the `BaseFab'.  They are
               completely analogous to the corresponding `Box' functions.
               There is no effect upon the array memory.
    */
    BaseFab<T>& shift (int idir,
                       int n_cell);

    /*@ManDoc: Perform shifts upon the domain of the `BaseFab'.  They are
               completely analogous to the corresponding `Box' functions.
               There is no effect upon the array memory.
    */
    BaseFab<T>& shiftHalf (int dir,
                           int num_halfs);

    /*@ManDoc: Perform shifts upon the domain of the `BaseFab'. They are
               completely analogous to the corresponding `Box' functions.
               There is no effect upon the array memory.
    */
    BaseFab<T>& shiftHalf (const IntVect& num_halfs);

protected:
    //
    // Allocates memory for the `BaseFab<T>'.
    //
    void define ();
    //
    // Deallocates memory for the `BaseFab<T>'.
    //
    void undefine ();
    //
    // The function called by BaseFab copy operations.
    //
    void performCopy (const BaseFab<T>& src,
                      const Box&        srcbox,
                      int               srccomp,
                      const Box&        destbox,
                      int               destcomp,
                      int               numcomp);

#ifdef BL_USE_POINTLIB
#ifndef BL_CRAY_BUG_DEFARG
    //
    // This function is called by the `PointBaseFab' copy
    // operations.  This function is only
    // defined if BL\_CRAY\_BUG\_DEFARG is NOT defined and if
    // BL_USE_POINTLIB is defined.
    //
    void performCopy (const PointBaseFab<PointDomain,T>& src,
                      const Box&                         srcbox,
                      int                                srccomp,
                      const Box&                         destbox,
                      int                                destcomp,
                      int                                numcomp);
#endif
#endif
    //
    // This function is called by the `BaseFab' setVal operations.
    //
    void performSetVal (T          x,
                        const Box& bx,
                        int        nstart,
                        int        numcomp);
private:
    //
    // These functions are made private to prevent use of the default 
    // functions provided by the C++ compiler.
    //
    BaseFab<T>& operator= (const BaseFab<T>&);
    BaseFab (const BaseFab<T>&);

protected:

    Box  domain;   // My index space.
    int  nvar;     // Number components.
    long numpts;   // Cached number of points in FAB.
    long truesize; // nvar*numpts that was allocated on heap.
    T*   dptr;     // The data pointer.
};

//
// Implementation.
//

template <class T>
inline
int
BaseFab<T>::nComp () const
{
    return nvar;
}

template <class T>
inline
const Box&
BaseFab<T>::box () const
{
    return domain;
}

template <class T>
inline
const int*
BaseFab<T>::length () const
{
    return domain.length().getVect();
}

template <class T>
inline
const IntVect&
BaseFab<T>::smallEnd () const
{
    return domain.smallEnd();
}

template <class T>
inline
const IntVect&
BaseFab<T>::bigEnd () const
{
    return domain.bigEnd();
}

template <class T>
inline
const int*
BaseFab<T>::loVect () const
{
    return domain.loVect();
}

template <class T>
inline
const int*
BaseFab<T>::hiVect () const
{
    return domain.hiVect();
}

template <class T>
inline
bool
BaseFab<T>::contains (const BaseFab<T>& fab) const
{
    return box().contains(fab.box()) && nvar <= fab.nvar;
}

template <class T>
inline
bool
BaseFab<T>::contains (const Box& bx) const
{
    return box().contains(bx);
}

template <class T>
inline
T*
BaseFab<T>::dataPtr (int n)
{
    assert(!(dptr == 0));
    return &dptr[n*numpts];
}

template <class T>
inline
const T*
BaseFab<T>::dataPtr (int n) const
{
    assert(!(dptr == 0));
    return &dptr[n*numpts];
}

template <class T>
inline
T&
BaseFab<T>::operator() (const IntVect& p,
                        int            n)
{
    assert(n >= 0);
    assert(n < nvar);
    assert(!(dptr == 0));
    assert(domain.contains(p));

    return dptr[domain.index(p)+n*numpts];
}

template <class T>
inline
T&
BaseFab<T>::operator() (const IntVect& p)
{
    assert(!(dptr == 0));
    assert(domain.contains(p));

    return dptr[domain.index(p)];
}

template <class T>
inline
const T&
BaseFab<T>::operator() (const IntVect& p,
                        int            n) const
{
    assert(n >= 0);
    assert(n < nvar);
    assert(!(dptr == 0));
    assert(domain.contains(p));

    return dptr[domain.index(p)+n*numpts];
}

template <class T>
inline
const T&
BaseFab<T>::operator() (const IntVect& p) const
{
    assert(!(dptr == 0));
    assert(domain.contains(p));

    return dptr[domain.index(p)];
}

template <class T>
inline
void
BaseFab<T>::getVal  (T*             data,
                     const IntVect& pos,
                     int            n,
                     int            numcomp) const
{
    const int loc      = domain.index(pos);
    const long size    = domain.numPts();

    assert(!(dptr == 0));
    assert(n >= 0 && n + numcomp <= nvar);

    for (int k = 0; k < numcomp; k++)
        data[k] = dptr[loc+(n+k)*size];
}

template <class T>
inline
void
BaseFab<T>::getVal (T*             data,
                    const IntVect& pos) const
{
    getVal(data,pos,0,nvar);
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::shift (const IntVect& v)
{
    domain += v;
    return *this;
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::shift (int idir,
                   int n_cell)
{
    domain.shift(idir,n_cell);
    return *this;
}

template <class T>
inline
BaseFab<T> &
BaseFab<T>::shiftHalf (const IntVect& v)
{
    domain.shiftHalf(v);
    return *this;
}

template <class T>
inline
BaseFab<T> &
BaseFab<T>::shiftHalf (int idir,
                       int n_cell)
{
    domain.shiftHalf(idir,n_cell);
    return *this;
}

template <class T>
inline
void
BaseFab<T>::setVal (T val)
{
    performSetVal(val,box(), 0, nvar);
}

template <class T>
inline
void
BaseFab<T>::setVal (T          x,
                    const Box& bx,
                    int        n)
{
    performSetVal(x,bx,n,1);
}

template <class T>
inline
void
BaseFab<T>::setVal (T   x,
                    int n)
{
    performSetVal(x,domain,n,1);
}

template <class T>
inline
void
BaseFab<T>::setVal (T          x,
                    const Box& b,
                    int        ns,
                    int        num)
{
    performSetVal(x,b,ns,num);
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src,
                  const Box&        srcbox,
                  int               srccomp,
                  const Box&        destbox,
                  int               destcomp,
                  int               numcomp)
{
    assert(srcbox.sameSize(destbox));
    assert(src.box().contains(srcbox));
    assert(domain.contains(destbox));
    assert(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= nvar);
    performCopy(src,srcbox,srccomp,destbox,destcomp,numcomp);
    return *this;
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src)
{
    assert(nvar <= src.nvar);
    assert(domain.sameType(src.domain));
    Box overlap(domain);
    overlap &= src.domain;
    if (overlap.ok())
        performCopy(src,overlap,0,overlap,0,nvar);
    return *this;
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src,
                  const Box&        destbox)
{
    assert(nvar <= src.nvar);
    assert(domain.contains(destbox));
    Box overlap(destbox);
    overlap &= src.domain;
    if (overlap.ok())
        performCopy(src,overlap,0,overlap,0,nvar);
    return *this;
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::copy (const BaseFab<T>& src,
                  int               srccomp,
                  int               destcomp,
                  int               numcomp)
{
    assert(srccomp >= 0 && srccomp + numcomp <= src.nvar);
    assert(destcomp >= 0 && destcomp + numcomp <= nvar);
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
    assert(nvar > 0);
    assert(dptr == 0);
    assert(numpts > 0);
    assert(!(The_FAB_Arena == 0));

    truesize = nvar*numpts;
    dptr     = static_cast<T*>(The_FAB_Arena->alloc(truesize*sizeof(T)));
    //
    // Now call T::T() on the raw memory so we have valid Ts.
    //
    T* ptr = dptr;

    for (int i = 0; i < truesize; i++, ptr++)
        new (ptr) T;
}

template <class T>
void
BaseFab<T>::undefine ()
{
    assert(!(The_FAB_Arena == 0));
    //
    // Call T::~T() on the to-be-destroyed memory.
    //
    T* ptr = dptr;

    for (int i = 0; i < truesize; i++, ptr++)
        ptr->~T();

    The_FAB_Arena->free(dptr);

    dptr = 0;
}

template <class T>
inline
BaseFab<T>::BaseFab ()
    :
    domain(Box()),
    nvar(0),
    numpts(0),
    truesize(0),
    dptr(0)
{}

template <class T>
inline
BaseFab<T>::BaseFab (const Box& bx,
                     int        n)
    :
    domain(bx),
    nvar(n),
    numpts(bx.numPts()),
    dptr(0)
{
    define();
}

#ifdef BL_USE_POINTLIB
#ifndef BL_CRAY_BUG_DEFARG
template <class T>
BaseFab<T>::BaseFab (const PointBaseFab<PointDomain,T>& pbf,
                     T                                  val)
    :
    domain(pbf.minimalBox()),
    nvar(pbf.nComp()),
    numpts(domain.numPts()),
    dptr(0)
{
    define();
    setVal(val);
    PointForAllCX(PointDomain,T,pbf)
    {
        this->operator()(ivR,nR) = pbfR;
    } EndPointForAll
}
#endif /*BL_CRAY_BUG_DEFARG*/
#endif /*BL_USE_POINTLIB*/

template <class T>
void
BaseFab<T>::resize (const Box& b,
                    int        n)
{
    nvar   = n;
    domain = b;
    numpts = domain.numPts();

    if (dptr == 0)
    {
        define();
    }
    else if (nvar*numpts > truesize)
    {
        undefine();

        define();
    }
}

template <class T>
inline
BaseFab<T>::~BaseFab ()
{
    undefine();
}

template <class T>
inline
void
BaseFab<T>::clear ()
{
    undefine();

    dptr   = 0;
    domain = Box();
    nvar   = 0;
    numpts = 0;
}

//
// performCopy() has been rewritten here so we can insert pragma's which
// will enhance vectorization on the Cray's.  The downside is that the
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
    assert(src.box().contains(srcbox));
    assert(box().contains(destbox));
    assert(destbox.sameSize(srcbox));
    assert(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= nComp());

#if (BL_SPACEDIM == 1)
{                                                                       
    assert((destcomp) >= 0 && (destcomp) + (numcomp) <= nComp());
    assert((srccomp) >= 0 && (srccomp) + (numcomp) <= (src).nComp());
    Box _subbox_ = box(); 
    _subbox_ &= destbox; 
    assert(srcbox.sameSize(_subbox_)); 
    if (_subbox_.ok()) 
    { 
        const int *_th_plo = loVect(); 
        const int *_th_plen = length(); 
        const int *_x_plo = (src).loVect(); 
        const int *_x_plen = (src).length(); 
        const int *_subbox_lo = _subbox_.loVect(); 
        const int *_subbox_len = _subbox_.length().getVect(); 
        const int *_bx_lo = (srcbox).loVect(); 
        const int *_bx_len = (srcbox).length().getVect(); 
        T* _th_p = dataPtr(destcomp); 
        const T* _x_p  = (src).dataPtr(srccomp); 
        for (int _n = 0; _n < (numcomp); ++_n)
        { 
            T *_th_pp = _th_p + ((_subbox_lo[0]-_th_plo[0])+_n*_th_plen[0]); 
            const T *_x_pp = _x_p + ((_bx_lo[0]-_x_plo[0])+_n*_x_plen[0]);
#ifdef BL_ARCH_CRAY
#pragma _CRI ivdep
#endif
            for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)
            { 
                int iR = _i + _subbox_lo[0]; iR += 0; 
                int ixR = _i + _bx_lo[0]; ixR += 0; 
                T &thisR = * _th_pp; const T & srcR = _x_pp[_i];
#elif (BL_SPACEDIM == 2)
{                                                                       
    assert((destcomp) >= 0 && (destcomp) + (numcomp) <= nComp()); 
    assert((srccomp) >= 0 && (srccomp) + (numcomp) <= (src).nComp()); 
    Box _subbox_ = box(); 
    _subbox_ &= destbox; 
    assert(srcbox.sameSize(_subbox_)); 
    if (_subbox_.ok())
    { 
        const int *_th_plo = loVect(); 
        const int *_th_plen = length(); 
        const int *_x_plo = (src).loVect(); 
        const int *_x_plen = (src).length(); 
        const int *_subbox_lo = _subbox_.loVect(); 
        const int *_subbox_len = _subbox_.length().getVect(); 
        const int *_bx_lo = (srcbox).loVect(); 
        const int *_bx_len = (srcbox).length().getVect(); 
        T* _th_p = dataPtr(destcomp); 
        const T* _x_p  = (src).dataPtr(srccomp); 
        for (int _n = 0; _n < (numcomp); ++_n)
        { 
            int nR = _n + destcomp; nR += 0; 
            int nxR = _n + srccomp; nxR += 0; 
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
#ifdef BL_ARCH_CRAY
#     pragma _CRI ivdep
#endif
                for (int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp)
                {  
                    T &thisR = * _th_pp; const T & srcR = _x_pp[_i]; 
#elif (BL_SPACEDIM == 3)
{ 
    assert((destcomp) >= 0 && (destcomp) + (numcomp) <= nComp()); 
    assert((srccomp) >= 0 && (srccomp) + (numcomp) <= (src).nComp()); 
    Box _subbox_(box()); 
    _subbox_ &= destbox; 
    assert((srcbox).sameSize(_subbox_)); 
    if (_subbox_.ok())
    { 
        const int *_th_plo = loVect(); 
        const int *_th_plen = length(); 
        const int *_x_plo = (src).loVect(); 
        const int *_x_plen = (src).length(); 
        const int *_subbox_lo = _subbox_.loVect(); 
        const int *_subbox_len = _subbox_.length().getVect(); 
        const int *_bx_lo = (srcbox).loVect(); 
        const int *_bx_len = (srcbox).length().getVect(); 
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
#ifdef BL_ARCH_CRAY
#     pragma _CRI ivdep
#endif
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

#ifdef BL_USE_POINTLIB
#ifndef BL_CRAY_BUG_DEFARG
template <class T>
void
BaseFab<T>::performCopy (const PointBaseFab<PointDomain,T>& src,
                         const Box&                         srcbox,
                         int                                srccomp,
                         const Box&                         destbox,
                         int                                destcomp,
                         int                                numcomp)
{
    assert(box().contains(destbox));
    assert(destbox.sameSize(srcbox));
    assert(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= nComp());

    destcomp -= srccomp;
    Box bx(domain);
    bx &= srcbox;
    bx &= destbox;
    if (bx.ok())
    {
        PointForAllCXBNN(PointDomain,T,src,bx,srccomp,numcomp)
        {
            this->operator()(ivR, destcomp + nR) = srcR;
        } EndPointForAll
    }
}
#endif /*!BL_CRAY_BUG_DEFARG*/
#endif /*BL_USE_POINTLIB*/

template <class T>
void
BaseFab<T>::performSetVal (T         val,
                           const Box& bx,
                           int        ns,
                           int        num)
{
    assert(domain.contains(bx));
    assert(ns >= 0 && ns + num <= nvar);

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
void
BaseFab<T>::setComplement (T          x,
                           const Box& b,
                           int        ns,
                           int        num)
{
    BoxList b_lst = boxDiff(domain,b);
    for (BoxListIterator bli(b_lst); bli; ++bli)
        performSetVal(x,bli(),ns,num);
}

#endif /*BL_BASEFAB_H*/
