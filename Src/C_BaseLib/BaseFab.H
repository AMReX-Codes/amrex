//BL_COPYRIGHT_NOTICE

#ifndef BL_BASEFAB_H
#define BL_BASEFAB_H

//
// $Id: BaseFab.H,v 1.11 1997-12-09 23:29:59 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
#else
#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#endif

#include <Assert.H>
#include <Boolean.H>
#include <Box.H>
#include <BoxList.H>
#include <Looping.H>
#include <Misc.H>
#include <PArray.H>
#include <REAL.H>

#ifdef BL_USE_POINTLIB
#ifndef BL_CRAY_BUG_DEFARG
#include <PointBaseFab.H>
#include <PointDomain.H>
#endif /*BL_CRAY_BUG_DEFARG*/
#endif /*BL_USE_POINTLIB*/

#ifdef BL_USE_POINTLIB
#include <PointLooping.H>
#endif

template <class T> class BaseFab;

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

  This class does NOT provide a copy constructor or assignment operator.
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
    //@ManDoc: The (virtual) destructor deletes the array memory.
    //
    virtual ~BaseFab ();

    /*@ManDoc: This function resizes a `BaseFab' so it covers the `Box' B
               with N components.  The default action
               is that under resizing, the memory allocated for the `BaseFab'
               only grows and never shrinks.  This function is particularly
               useful when a `BaseFab' is used as a temporary space which
               must be a different size whenever it is used.  Resize'ing a
               temp will often be faster than re-allocating a `BaseFab'
               because memory allocation can often be avoided.
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

    /*@ManDoc: Compute the Lp-norm of this FAB using components
               (comp : comp+numcomp-1).
               p < 0  -> ERROR
               p = 0  -> infinity norm (max norm)
               p = 1  -> sum of ABS(FAB)
    */
    virtual Real norm (int p,
                       int comp=0,
                       int numcomp=1) const;
    //
    //@ManDoc: Same as above except only on given subbox.
    //
    virtual Real norm (const Box& subbox,
                       int        p,
                       int        comp=0,
                       int        numcomp=1) const;
    //
    //@ManDoc:  Compute absolute value for all components of this FAB.
    //
    void abs ();
    //
    //@ManDoc: Same as above except only for components (comp: comp+numcomp-1)
    //
    void abs (int comp,
              int numcomp=1);
    //
    //@ManDoc: Calculate abs() on subbox for given component range.
    //
    void abs (const Box& subbox,
              int        comp = 0,
              int        numcomp=1);
    //
    //@ManDoc: Minimum value of given component.
    //
    T min (int comp = 0) const;
    //
    //@ManDoc: Minimum value of given component in given subbox.
    //
    T min (const Box& subbox,
           int        comp = 0) const;
    //
    //@ManDoc: Maximum value of given component.
    //
    T max (int comp = 0) const;
    //
    //@ManDoc: Maximum value of given component in given subbox.
    //
    T max (const Box& subbox,
           int        comp = 0) const;
    //
    //@ManDoc: Find location of minimum value in given component.
    //
    IntVect minIndex (int comp = 0) const;

    /*@ManDoc: Find location of minimum value in given component in
               given subbox.
    */
    IntVect minIndex (const Box& subbox,
                      int        comp = 0) const;
    //
    //@ManDoc: Find location of maximum value in given component.
    //
    IntVect maxIndex (int comp = 0) const;
 
    /*@ManDoc: Find location of maximum value in given component in given
               subbox.
    */
    IntVect maxIndex (const Box& subbox,
                      int        comp = 0) const;
    //
    //@ManDoc: Fill with a pattern of numbers.
    //
    void patternFill (int mark = 0);

    /*@ManDoc: Copies with index reversal from srcbox region of src
               into destbox region of this FAB.  All components are
               copied and multiplied by corresponding multiplier.
    */
    void copyRev (const Box&        destbox,
                  const BaseFab<T>& src,
                  const Box&         srcbox,
                  int                reversal_index,
                  T*                 multiplier);

#if !defined(BL_CRAY_BUG_DEFARG)

    /*@ManDoc: Copies with a reversal.  The index that is
               reversed is reversal\_index.  Comp is the component copied.
               This is not defined if the macro BL\_CRAY\_BUG\_DEFARG is
               defined.
    */
    friend void copyFABRev (const BaseFab<T>& src,
                            const Box&        srcbox,
                            BaseFab<T>&       dest,
                            const Box&        destbox,
                            int               reversal_index,
                            int               comp,
                            T                 multiplier);
#endif /*!defined(BL_CRAY_BUG_DEFARG)*/
    //
    //@ManDoc: Compute sum of given component of FAB state vector.
    //
    T sum (int comp,
           int numcomp = 1) const;

    /*@ManDoc: Compute sum of given component of FAB state vector in
               given subbox.
    */
    T sum (const Box& subbox,
           int        comp,
           int        numcomp = 1) const;
    //
    //@ManDoc: Most general version, specify subbox and which components.
    //
    BaseFab<T>& invert (T          v,
                        const Box& subbox,
                        int        comp=0,
                        int        numcomp=1);
    //
    //@ManDoc: As above except on entire domain.
    //
    BaseFab<T>& invert (T   v,
                        int comp,
                        int numcomp=1);
    //
    //@ManDoc: As above except on entire domain, all components.
    //
    BaseFab<T>& invert (T v);
    //
    //@ManDoc: Negate BaseFab, most general.
    //
    BaseFab<T>& negate (const Box& subbox,
                        int        comp=0,
                        int        numcomp=1);
    //
    //@ManDoc: As above, except on entire domain.
    //
    BaseFab<T>& negate (int comp,
                        int numcomp=1);
    //
    //@ManDoc: As above, except on entire domain and all components.
    //
    BaseFab<T>& negate ();
    //
    //@ManDoc: Scalar addition (a[i] <- a[i] + r), most general.
    //
    BaseFab<T>& plus (T          r,
                      const Box& b,
                      int        comp=0,
                      int        numcomp=1);
    //
    //@ManDoc: As above, except on entire domain.
    //
    BaseFab<T>& plus (T   r,
                      int comp,
                      int numcomp=1);
    //
    //@ManDoc: As above, except on entire domain and all components.
    //
    BaseFab<T>& plus (T r);

    /*@ManDoc: Addition in place.  This will often be more efficient than
               making new BaseFab for result.
    */
    BaseFab<T>& operator+= (T r);
    //
    //@ManDoc: FAB addition (a[i] <- a[i] + b[i]) in place.
    //
    BaseFab<T>& operator+= (const BaseFab<T>& f);
    //
    //@ManDoc: FAB addition (a[i] <- a[i] + b[i]). The same as += operator.
    //
    BaseFab<T>& plus (const BaseFab<T>& src);

    /*@ManDoc: Add src components (srccomp:srccomp+numcomp-1) to
               this FAB's components (destcomp:destcomp+numcomp-1)
               where the two FABs intersect.
    */
    BaseFab<T>& plus (const BaseFab<T>& src,
                      int                srccomp,
                      int                destcomp,
                      int                numcomp=1);

    /*@ManDoc: Same as above except addition is restricted to intersection
               of subbox and src FAB. NOTE: subbox must be contained in this
               FAB.
    */
    BaseFab<T>& plus (const BaseFab<T>& src,
                      const Box&         subbox,
                      int                srccomp,
                      int                destcomp,
                      int                numcomp=1);

    /*@ManDoc: Add srcbox region of src FAB to destbox region of this FAB.
               The srcbox and destbox must be same size.
    */
    BaseFab<T>& plus (const BaseFab<T>& src,
                      const Box&         srcbox,
                      const Box&         destbox,
                      int                srccomp,
                      int                destcomp,
                      int                numcomp=1);

    /*@ManDoc: Scalar subtraction (a[i] <- a[i] - r).
               Note: use plus(-r) for more general operations.
    */
    BaseFab<T>& operator-= (T r);
    //
    //@ManDoc: FAB subtraction (a[i] <- a[i] - b[i]), in place.
    //
    BaseFab<T>& operator-= (const BaseFab<T>& f);
    //
    //@ManDoc: FAB subtraction (a[i] <- a[i] - b[i]). The same as -= operator.
    //
    BaseFab<T>& minus (const BaseFab<T>& src);

    /*@ManDoc: Subtract src components (srccomp:srccomp+numcomp-1) to
               this FAB's components (destcomp:destcomp+numcomp-1) where
               the two FABs intersect.
    */
    BaseFab<T>& minus (const BaseFab<T>& src,
                       int                srccomp,
                       int                destcomp,
                       int                numcomp=1);

    /*@ManDoc: Same as above except subtraction is restricted to intersection
               of subbox and src FAB.  NOTE: subbox must be contained in
               this FAB.
    */
    BaseFab<T>& minus (const BaseFab<T>& src,
                       const Box&         subbox,
                       int                srccomp,
                       int                destcomp,
                       int                numcomp=1);

    /*@ManDoc: Subtract srcbox region of src FAB from destbox region
               of this FAB. srcbox and destbox must be same size.
    */
    BaseFab<T>& minus (const BaseFab<T>& src,
                       const Box&         srcbox,
                       const Box&         destbox,
                       int                srccomp,
                       int                destcomp,
                       int                numcomp=1);
    //
    //@ManDoc: Scalar multiplication (a[i] <- a[i] * r), in place.
    //
    BaseFab<T>& operator*= (T r);
    //
    //@ManDoc: Scalar multiplication (a[i] <- a[i] * r).  The same as *=.
    //
    BaseFab<T>& mult (T r);

    /*@ManDoc: Scalar multiplication, except control which components
               are multiplied.
    */
    BaseFab<T>& mult (T   r,
                      int comp,
                      int numcomp=1);
    //
    //@ManDoc: As above, except specify sub-box.
    //
    BaseFab<T>& mult (T          r,
                      const Box& b,
                      int        comp=0,
                      int        numcomp=1);
    //
    //@ManDoc: FAB multiplication (a[i] <- a[i] * b[i]), in place.
    //
    BaseFab<T>& operator*= (const BaseFab<T>& f);
    //
    //@ManDoc: As above.
    //
    BaseFab<T>& mult (const BaseFab<T>& src);

    /*@ManDoc: Multiply src components (srccomp:srccomp+numcomp-1) with
               this FAB's components (destcomp:destcomp+numcomp-1) where
               the two FABs intersect.
    */
    BaseFab<T>& mult (const BaseFab<T>& src,
                      int                srccomp,
                      int                destcomp,
                      int                numcomp=1);

    /*@ManDoc: Same as above except multiplication is restricted to
               intersection of subbox and src FAB.  NOTE: subbox must be
               contained in this FAB.
    */
    BaseFab<T>& mult (const BaseFab<T>& src,
                      const Box&         subbox,
                      int                srccomp,
                      int                destcomp,
                      int                numcomp=1);

    /*@ManDoc: Multiply srcbox region of src FAB with destbox region
               of this FAB. The srcbox and destbox must be same size.
    */
    BaseFab<T>& mult (const BaseFab<T>& src,
                      const Box&         srcbox,
                      const Box&         destbox,
                      int                srccomp,
                      int                destcomp,
                      int                numcomp=1);
    //
    //@ManDoc: Scalar division (a[i] <- a[i] / r), in place.
    //
    BaseFab<T>& operator/= (T r);
    //
    //@ManDoc: Scalar division (a[i] <- a[i] / r), in place.
    //
    BaseFab<T>& divide (T r);
    //
    //@ManDoc: As above except specify which components.
    //
    BaseFab<T>& divide (T   r,
                        int comp,
                        int numcomp=1);
    //
    //@ManDoc: As above except specify sub-box.
    //
    BaseFab<T>& divide (T          r,
                        const Box& b,
                        int        comp=0,
                        int        numcomp=1);
    //
    //@ManDoc: FAB division, in place.
    //
    BaseFab<T>& operator/= (const BaseFab<T>& src);
    //
    //@ManDoc: Same as above.
    //
    BaseFab<T>& divide (const BaseFab<T>& src);

    /*@ManDoc: This FAB is numerator, src FAB is denominator
               divide src components (srccomp:srccomp+numcomp-1) into
               this FAB's components (destcomp:destcomp+numcomp-1)
               where the two FABs intersect.
    */
    BaseFab<T>& divide (const BaseFab<T>& src,
                        int                srccomp,
                        int                destcomp,
                        int                numcomp=1);

    /*@ManDoc: Same as above except division is restricted to
               intersection of subbox and src FAB.  NOTE: subbox must be
               contained in this FAB.
    */
    BaseFab<T>& divide (const BaseFab<T>& src,
                        const Box&         subbox,
                        int                srccomp,
                        int                destcomp,
                        int                numcomp=1);
 
    /*@ManDoc: destbox region of this FAB is numerator. srcbox regions of
               src FAB is denominator. srcbox and destbox must be same size.
    */
    BaseFab<T>& divide (const BaseFab<T>& src,
                        const Box&         srcbox,
                        const Box&         destbox,
                        int                srccomp,
                        int                destcomp,
                        int                numcomp=1);

    /*@ManDoc: Linear interpolation / extrapolation.
               Result is (t2-t)/(t2-t1)*f1 + (t-t1)/(t2-t1)*f2
               Data is taken from b1 region of f1, b2 region of f2
               and stored in b region of this FAB.
               Boxes b, b1 and b2 must be the same size.
               Data is taken from component comp1 of f1, comp2 of f2,
               and stored in component comp of this FAB.
               This FAB is returned as a reference for chaining.
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

    /*@ManDoc: Linear combination.  Result is alpha*f1 + beta*f2.
               Data is taken from b1 region of f1, b2 region of f2
               and stored in b region of this FAB.
               Boxes b, b1 and b2 must be the same size.
               Data is taken from component comp1 of f1, comp2 of f2,
               and stored in component comp of this FAB.
               This FAB is returned as a reference for chaining.
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
protected:
    //
    // Allocates memory for the `BaseFab<T>'.
    //
     void define ();
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
#endif /*BL_CRAY_BUG_DEFARG*/
#endif /*BL_USE_POINTLIB*/
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
    BaseFab (const BaseFab<T>& rhs);
    BaseFab<T>& operator= (const BaseFab<T>& rhs);

protected:

    Box  domain;   // My index space.
    int  nvar;     // Number components.
    long numpts;   // Cached number of points in FAB.
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
    assert(!(dptr == 0));
    return dptr[domain.index(p)+n*numpts];
}

template <class T>
inline
T&
BaseFab<T>::operator() (const IntVect& p)
{
    assert(!(dptr == 0));
    return dptr[domain.index(p)];
}

template <class T>
inline
const T&
BaseFab<T>::operator() (const IntVect& p,
                        int            n) const
{
    assert(!(dptr == 0));
    return dptr[domain.index(p)+n*numpts];
}

template <class T>
inline
const T&
BaseFab<T>::operator() (const IntVect& p) const
{
    assert(!(dptr == 0));
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
    domain  += v;
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
void
BaseFab<T>::define ()
{
    assert(nvar > 0);
    assert(numpts > 0);

    if (dptr == 0)
    {
        if ((dptr = new T[nvar*numpts]) == 0)
            BoxLib::OutOfMemory(__FILE__, __LINE__);
    }
}

template <class T>
inline
void
BaseFab<T>::clear ()
{
    if (!(dptr == 0))
        delete [] dptr;
    dptr   = 0;
    domain = Box();
    nvar   = 0;
    numpts = 0;
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
    assert(nvar == src.nvar);
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
    assert(nvar == src.nvar);
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
inline
BaseFab<T>::BaseFab ()
    :
    domain(Box()),
    nvar(0),
    numpts(0),
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
    resize(domain, nvar);
    setVal(val);
    PointForAllCX(PointDomain,T,pbf)
    {
        this->operator()(ivR,nR) = pbfR;
    } EndPointForAll
}
#endif /*BL_CRAY_BUG_DEFARG*/
#endif /*BL_USE_POINTLIB*/

//
// This isn't inlined as it's virtual.
//

template <class T>
BaseFab<T>::~BaseFab ()
{
    delete [] dptr;
}

template <class T>
void
BaseFab<T>::resize (const Box& b,
                    int        n)
{
    long old_size = nvar*numpts;

    nvar   = n;
    domain = b;
    numpts = domain.numPts();

    if (dptr == 0)
    {
        define();
    }
    else
    {
        long new_size = nvar*numpts;

        if (new_size > old_size)
        {
            T* ndptr = new T[new_size];
            if (ndptr == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            
            for (long i = 0; i < old_size; i++)
                ndptr[i] = dptr[i];

            delete [] dptr;

            dptr = ndptr;
        }
    }
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

template <class T>
void
BaseFab<T>::abs ()
{
    ForAllThis(Real)
    {
        thisR = Abs(thisR);
    } EndFor
}

template <class T>
void
BaseFab<T>::abs (int comp,
                 int numcomp)
{
    ForAllThisNN(Real,comp,numcomp)
    {
        thisR = Abs(thisR);
    } EndFor
}

template <class T>
void
BaseFab<T>::abs (const Box& subbox,
                 int        comp,
                 int        numcomp)
{
    ForAllThisBNN(Real,subbox,comp,numcomp)
    {
        thisR = Abs(thisR);
    } EndFor
}

//
// This isn't inlined as it's virtual.
//

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
    assert(comp >= 0 && comp+numcomp <= nComp());
    assert(p >= 0);

    Real* tmp  = 0;
    int tmplen = 0;
    Real nrm   = 0;
    if (p == 0)
    {
        ForAllThisCPencil(T,subbox,comp,numcomp)
        {
            const T* row = &thisR;
            if (tmp == 0)
            {
                tmp = new Real[thisLen];
                if (tmp == 0)
                    BoxLib::OutOfMemory(__FILE__, __LINE__);
                tmplen = thisLen;
                for (int i = 0; i < thisLen; i++)
                {
                    tmp[i] = Abs(Real(row[i]));
		}
            }
            else
            {
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = Max(tmp[i],Real(Abs(row[i])));
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
            nrm = Max(nrm, tmp[i]);
    }
    else if (p == 1)
    {
        ForAllThisCPencil(T,subbox,comp,numcomp)
        {
            const T* row = &thisR;
            if (tmp == 0)
            {
                tmp = new Real[thisLen];
                if (tmp == 0)
                    BoxLib::OutOfMemory(__FILE__, __LINE__);
                tmplen = thisLen;
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = Abs(Real(row[i]));
            }
            else
            {
                for (int i = 0; i < thisLen; i++)
                    tmp[i] += Abs(Real(row[i]));
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
            nrm += tmp[i];
    }
    else
      BoxLib::Error("BaseFab::norm(): only p == 0 or p == 1 are supported");

    delete [] tmp;

    return nrm;
}

#ifdef BL_ARCH_CRAY
#define RESTRICT restrict
//
// Template specialization for Real.
//
Real
BaseFab<Real>::norm (const Box& subbox,
                     int        p,
                     int        comp,
                     int        numcomp) const
{
    assert(comp >= 0 && comp+numcomp <= nComp());
    assert(p >= 0);

    Real* RESTRICT tmp  = 0;
    int tmplen = 0;
    Real nrm   = 0;
    if (p == 0)
    {
        ForAllThisCPencil(Real,subbox,comp,numcomp)
        {
            const Real* RESTRICT row = &thisR;
            if (tmp == 0)
            {
                tmp = new Real[thisLen];
                if (tmp == 0)
                    BoxLib::OutOfMemory(__FILE__, __LINE__);
                tmplen = thisLen;
#pragma _CRI ivdep
                for (int i = 0; i < thisLen; i++)
                {
                    tmp[i] = fabs(row[i]);
		}
            }
            else
            {
#pragma _CRI ivdep
	      for (int i = 0; i < thisLen; i++)
              {
		Real a = fabs(row[i]);
		tmp[i] = tmp[i] > a ? tmp[i] : a ;
	      }
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
        {
	    Real a = tmp[i];
	    nrm = nrm > a ? nrm : a ;
	}
    }
    else if (p == 1)
    {
        ForAllThisCPencil(Real,subbox,comp,numcomp)
        {
            const Real* row = &thisR;
            if (tmp == 0)
            {
                tmp = new Real[thisLen];
                if (tmp == 0)
                    BoxLib::OutOfMemory(__FILE__, __LINE__);
                tmplen = thisLen;
#pragma _CRI ivdep
                for (int i = 0; i < thisLen; i++)
                {
		    tmp[i] = fabs(row[i]);
		}
            }
            else
            {
#pragma _CRI ivdep
	      for (int i = 0; i < thisLen; i++)
              {
		    tmp[i] += fabs( row[i] );
	      }
            }
        } EndForPencil
        nrm = tmp[0];
#pragma _CRI ivdep
        for (int i = 1; i < tmplen; i++)
            nrm += tmp[i];
    }
    else
      BoxLib::Error("BaseFab<Real>::norm(): only p == 0 or p == 1 are supported");

    delete [] tmp;

    return nrm;
}
#endif /*BL_ARCH_CRAY*/

template <class T>
T
BaseFab<T>::min (int comp) const
{
    T *_min_row = 0;
    int _X_len = 0;
    ForAllThisCPencil(T,domain,comp,1)
    {
        const T* _row = &thisR;
        if (_min_row == 0)
        {
            _min_row = new T[thisLen];
            if (_min_row == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            _X_len = thisLen;
            for (int i = 0; i < thisLen; i++)
                _min_row[i] = _row[i];
        }
        else
        {
            for (int i = 0; i < thisLen; i++)
                _min_row[i] = Min(_row[i],_min_row[i]);
        }
    } EndForPencil;

    T _min = _min_row[0];
    for (int i = 1; i < _X_len; i++)
        _min = Min(_min,_min_row[i]);

    delete [] _min_row;

    return _min;
}

template <class T>
T
BaseFab<T>::min (const Box& subbox,
                 int        comp) const
{
    T *_min_row = 0;
    int _X_len = 0;
    ForAllThisCPencil(T,subbox,comp,1)
    {
        const T* _row = &thisR;
        if (_min_row == 0)
        {
            _min_row = new T[thisLen];
            if (_min_row == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            _X_len = thisLen;
            for (int i = 0; i < thisLen; i++)
                _min_row[i] = _row[i];
        }
        else
        {
            for (int i = 0; i < thisLen; i++)
                _min_row[i] = Min(_row[i],_min_row[i]);
        }
    } EndForPencil;

    T _min = _min_row[0];
    for (int i = 1; i < _X_len; i++)
        _min = Min(_min,_min_row[i]);

    delete [] _min_row;

    return _min;
}

template <class T>
T
BaseFab<T>::max (int comp) const
{
    T *_max_row = 0;
    int _X_len = 0;
    ForAllThisCPencil(T,domain,comp,1)
    {
        const T* _row = &thisR;
        if (_max_row== 0)
        {
            _max_row = new T[thisLen];
            if (_max_row == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            _X_len = thisLen;
            for (int i = 0; i < thisLen; i++)
                _max_row[i] = _row[i];
        }
        else
        {
            for (int i = 0; i < thisLen; i++)
                _max_row[i] = Max(_row[i],_max_row[i]);
        }
    } EndForPencil;

    T _max = _max_row[0];
    for (int i = 1; i < _X_len; i++)
        _max = Max(_max,_max_row[i]);

    delete [] _max_row;

    return _max;
}

template <class T>
T
BaseFab<T>::max (const Box& subbox,
                 int        comp) const
{
    T *_max_row = 0;
    int _X_len = 0;
    ForAllThisCPencil(T,subbox,comp,1)
    {
        const T* _row = &thisR;
        if (_max_row == 0)
        {
            _max_row = new T[thisLen];
            if (_max_row == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            _X_len = thisLen;
            for (int i = 0; i < thisLen; i++)
                _max_row[i] = _row[i];
        }
        else
        {
            for (int i = 0; i < thisLen; i++)
                _max_row[i] = Max(_row[i],_max_row[i]);
        }
    } EndForPencil;

    T _max = _max_row[0];
    for (int i = 1; i < _X_len; i++)
        _max = Max(_max,_max_row[i]);

    delete [] _max_row;

    return _max;
}

template <class T>
IntVect
BaseFab<T>::minIndex (int comp) const
{
    IntVect _min_loc(domain.smallEnd());
    T _min_val = (*this).operator()(_min_loc,comp);
    ForAllThisCBNN(T,domain,comp,1)
    {
        if (thisR < _min_val)
        {
            _min_val = thisR;
            D_EXPR(_min_loc[0] = iR,
                   _min_loc[1] = jR,
                   _min_loc[2] = kR);
        }
    } EndFor;

    return _min_loc;
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
            D_EXPR(_min_loc[0] = iR,
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
    IntVect _max_loc(domain.smallEnd());
    T _max_val = (*this).operator()(_max_loc,comp);
    ForAllThisCBNN(T,domain,comp,1)
    {
        if (thisR > _max_val)
        {
            _max_val = thisR;
            D_EXPR(_max_loc[0] = iR,
                   _max_loc[1] = jR,
                   _max_loc[2] = kR);
        }
    } EndFor;

    return _max_loc;
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
            D_EXPR(_max_loc[0] = iR,
                   _max_loc[1] = jR,
                   _max_loc[2] = kR);
        }
    } EndFor;

    return _max_loc;
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::plus (T r)
{
    return operator+=(r);
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& x)
{
    return operator+=(x);
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::operator-= (T r)
{
    return operator+=(-r);
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& x)
{
    return operator-=(x);
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::mult (T r)
{
    return operator*=(r);
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& x)
{
    return operator*=(x);
}

template <class T>
inline
BaseFab<T>&
BaseFab<T>::divide (T r)
{
    return operator/=(r);
}

template <class T>
inline
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
        thisR = D_TERM(iR*100, +jR*10, + kR) + 1000*nR + 10000*mark;
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
void
copyFABRev (const BaseFab<T>& src,
            const Box&        srcbox,
            BaseFab<T>&       dest,
            const Box&        destbox,
            int               reversal_index,
            int               comp,
            T                 multiplier)
{
    ForAllRevXBNYCBNNN(T,dest,destbox,comp,src,srcbox,comp,1,reversal_index)
    {
        destR = multiplier * srcR;
    } EndFor
}

template <class T>
T
BaseFab<T>::sum (int comp,
                 int numcomp) const
{
    T *_sum_row = 0;
    int _sum_len = 0;
    ForAllThisCPencil(T,domain,comp,numcomp)
    {
        const T* _row = &thisR;
        if (_sum_row == 0)
        {
            _sum_row = new T[thisLen];
            if (_sum_row == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            _sum_len = thisLen;
            for (int i = 0; i < thisLen; i++)
                _sum_row[i] = _row[i];
        }
        else
        {
            for (int i = 0; i < thisLen; i++)
                _sum_row[i] += _row[i];
        }
    } EndForPencil;

    T _sum = _sum_row[0];
    for (int i = 1; i < _sum_len; i++)
        _sum += _sum_row[i];
    delete [] _sum_row;

    return _sum;
}

template <class T>
T
BaseFab<T>::sum (const Box& subbox,
                 int        comp,
                 int        numcomp) const
{
    T *_sum_row = 0;
    int _sum_len = 0;
    ForAllThisCPencil(T,subbox,comp,numcomp)
    {
        const T* _row = &thisR;
        if (_sum_row == 0)
        {
            _sum_row = new T[thisLen];
            if (_sum_row == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
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
    ForAllThis(T)
    {
        thisR = - thisR;
    } EndFor
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::negate (int comp,
                    int numcomp)
{
    ForAllThisNN(T,comp,numcomp)
    {
        thisR = - thisR;
    } EndFor
    return *this;
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
    ForAllThis(T)
    {
        thisR = r/thisR;
    } EndFor
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::invert (T   r,
                    int comp,
                    int numcomp)
{
    ForAllThisNN(T,comp,numcomp)
    {
        thisR = r/thisR;
    } EndFor
    return *this;
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
BaseFab<T>::operator += (T r)
{
    ForAllThis(T)
    {
        thisR += r;
    } EndFor
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (T   r,
                  int comp,
                  int numcomp)
{
    ForAllThisNN(T,comp,numcomp)
    {
        thisR += r;
    } EndFor
    return *this;
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
BaseFab<T>::operator += (const BaseFab<T>& x)
{
    ForAllThisXC(T,x)
    {
        thisR += xR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src,
                  int                srccomp,
                  int                destcomp,
                  int                numcomp)
{
    ForAllThisBNNXC(T,domain,destcomp,numcomp,src,srccomp)
    {
        thisR += srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src,
                  const Box&         subbox,
                  int                srccomp,
                  int                destcomp,
                  int                numcomp)
{
    ForAllThisBNNXC(T,subbox,destcomp,numcomp,src,srccomp)
    {
        thisR += srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::plus (const BaseFab<T>& src,
                  const Box&         srcbox,
                  const Box&         destbox,
                  int                srccomp,
                  int                destcomp,
                  int                numcomp)
{
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        thisR += srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator -= (const BaseFab<T>& x)
{
    ForAllThisXC(T,x)
    {
        thisR -= xR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src,
                   int                srccomp,
                   int                destcomp,
                   int                numcomp)
{
    ForAllThisBNNXC(T,domain,destcomp,numcomp,src,srccomp)
    {
        thisR -= srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src,
                   const Box&         subbox,
                   int                srccomp,
                   int                destcomp,
                   int                numcomp)
{
    ForAllThisBNNXC(T,subbox,destcomp,numcomp,src,srccomp)
    {
        thisR -= srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::minus (const BaseFab<T>& src,
                   const Box&         srcbox,
                   const Box&         destbox,
                   int                srccomp,
                   int                destcomp,
                   int                numcomp)
{
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        thisR -= srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator *= (T r)
{
    ForAllThis(T)
    {
        thisR *= r;
    } EndFor
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (T   r,
                  int comp,
                  int numcomp)
{
    ForAllThisNN(T,comp,numcomp)
    {
        thisR *= r;
    } EndFor
    return *this;
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
BaseFab<T>::operator *= (const BaseFab<T> &x)
{
    ForAllThisXC(T,x)
    {
        thisR *= xR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src,
                  int                srccomp,
                  int                destcomp,
                  int                numcomp)
{
    ForAllThisBNNXC(T,domain,destcomp,numcomp,src,srccomp)
    {
        thisR *= srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src,
                  const Box&         subbox,
                  int                srccomp,
                  int                destcomp,
                  int                numcomp)
{
    ForAllThisBNNXC(T,subbox,destcomp,numcomp,src,srccomp)
    {
        thisR *= srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::mult (const BaseFab<T>& src,
                  const Box&         srcbox,
                  const Box&         destbox,
                  int                srccomp,
                  int                destcomp,
                  int                numcomp)
{
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
        thisR *= srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator /= (T r)
{
    ForAllThis(T)
    {
        thisR /= r;
    } EndFor
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (T   r,
                    int comp,
                    int numcomp)
{
    ForAllThisNN(T,comp,numcomp)
    {
        thisR /= r;
    } EndFor
    return *this;
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
BaseFab<T>::operator /= (const BaseFab<T> &x)
{
    ForAllThisXC(T,x)
    {
        thisR /= xR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src,
                    int                srccomp,
                    int                destcomp,
                    int                numcomp)
{
    ForAllThisBNNXC(T,domain,destcomp,numcomp,src,srccomp)
    {
        thisR /= srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src,
                    const Box&         subbox,
                    int                srccomp,
                    int                destcomp,
                    int                numcomp)
{
    ForAllThisBNNXC(T,subbox,destcomp,numcomp,src,srccomp)
    {
        thisR /= srcR;
    } EndForTX
    return *this;
}

template <class T>
BaseFab<T>&
BaseFab<T>::divide (const BaseFab<T>& src,
                    const Box&         srcbox,
                    const Box&         destbox,
                    int                srccomp,
                    int                destcomp,
                    int                numcomp)
{
    ForAllThisBNNXCBN(T,destbox,destcomp,numcomp,src,srcbox,srccomp)
    {
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
    ForAllThisBNNXCBNYCBN(T,b,comp,numcomp,f1,b1,comp1,f2,b2,comp2)
    {
        thisR = (T) (alpha*Real(f1R) + beta*Real(f2R));
    } EndForTX
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
BaseFab<T>::linComb (const BaseFab<T>& f1,
                     const Box&        b1,
                     int               comp1,
                     const BaseFab<T>& f2,
                     const Box&        b2,
                     int               comp2,
                     Real              alpha,
                     Real              beta,
                     const Box&        b,
                     int               comp,
                     int               numcomp)

{
    ForAllThisBNNXCBNYCBN(T,b,comp,numcomp,f1,b1,comp1,f2,b2,comp2)
    {
        thisR = (T) (alpha*Real(f1R) + beta*Real(f2R));
    } EndForTX
    return *this;
}

#endif /*BL_BASEFAB_H*/
