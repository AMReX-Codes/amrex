#ifndef BL_BASEFAB_H
#define BL_BASEFAB_H
#include <winstd.H>

//
// $Id: BaseFab.H,v 1.37 2001-07-24 18:12:34 car Exp $
//

#include <cstdlib>

#include <BLassert.H>
#include <Box.H>
#include <BoxList.H>
#include <CArena.H>
#include <Looping.H>
#include <REAL.H>

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
    
    typedef T value_type;
    
    /*@ManDoc: Construct an invalid `BaseFab'.  The domain is invalid, the
               number of components is zero, and no actually array memory is
               allocated. An invalid `BaseFab' must be resized
               (see `BaseFab::resize') before use.
    */
    BaseFab ();

    BaseFab(const BaseFab&);
    BaseFab& operator=(const BaseFab&);
    //
    //@ManDoc: Make BaseFab with desired domain and number of components.
    //
    explicit BaseFab (const Box& bx,
                      int        n = 1);
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
    //
    //@ManDoc: Returns true if the data for the FAB has been allocated.
    //
    bool isAllocated () const;

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
               (scomp : scomp+ncomp-1).
               p < 0  -> ERROR
               p = 0  -> infinity norm (max norm)
               p = 1  -> sum of ABS(FAB)
    */
    Real norm (int p,
               int scomp = 0,
               int ncomp = 1) const;
    //
    //@ManDoc: Same as above except only on given subbox.
    //
    Real norm (const Box& subbox,
               int        p,
               int        scomp = 0,
               int        ncomp = 1) const;
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

    /*@ManDoc: Compute mask array with value of 1 in cells where
               BaseFab has value less than val, 0 otherwise.
               mask is resized by this function.
               The number of cells marked with 1 returned.
    */
    int maskLT (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
 
    /*@ManDoc: Same as above except mark cells with value less than
               or equal to val.
    */
    int maskLE (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
    //
    //@ManDoc: Same as above except mark cells with value equal to val.
    //
    int maskEQ (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
    //
    //@ManDoc: Same as above except mark cells with value greater than val.
    //
    int maskGT (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;
 
    /*@ManDoc: Same as above except mark cells with value greater than
               or equal to val.
    */
    int maskGE (BaseFab<int>& mask,
                T             val,
                int           comp = 0) const;

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
                  const Box&        srcbox,
                  int               reversal_index,
                  T*                multiplier);
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
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);

    /*@ManDoc: Same as above except addition is restricted to intersection
               of subbox and src FAB. NOTE: subbox must be contained in this
               FAB.
    */
    BaseFab<T>& plus (const BaseFab<T>& src,
                      const Box&        subbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);

    /*@ManDoc: Add srcbox region of src FAB to destbox region of this FAB.
               The srcbox and destbox must be same size.
    */
    BaseFab<T>& plus (const BaseFab<T>& src,
                      const Box&        srcbox,
                      const Box&        destbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);

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
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);

    /*@ManDoc: Same as above except subtraction is restricted to intersection
               of subbox and src FAB.  NOTE: subbox must be contained in
               this FAB.
    */
    BaseFab<T>& minus (const BaseFab<T>& src,
                       const Box&        subbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);

    /*@ManDoc: Subtract srcbox region of src FAB from destbox region
               of this FAB. srcbox and destbox must be same size.
    */
    BaseFab<T>& minus (const BaseFab<T>& src,
                       const Box&        srcbox,
                       const Box&        destbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1);
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
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);

    /*@ManDoc: Same as above except multiplication is restricted to
               intersection of subbox and src FAB.  NOTE: subbox must be
               contained in this FAB.
    */
    BaseFab<T>& mult (const BaseFab<T>& src,
                      const Box&        subbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);

    /*@ManDoc: Multiply srcbox region of src FAB with destbox region
               of this FAB. The srcbox and destbox must be same size.
    */
    BaseFab<T>& mult (const BaseFab<T>& src,
                      const Box&        srcbox,
                      const Box&        destbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp=1);
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
                        int               srccomp,
                        int               destcomp,
                        int               numcomp=1);

    /*@ManDoc: Same as above except division is restricted to
               intersection of subbox and src FAB.  NOTE: subbox must be
               contained in this FAB.
    */
    BaseFab<T>& divide (const BaseFab<T>& src,
                        const Box&        subbox,
                        int               srccomp,
                        int               destcomp,
                        int               numcomp=1);
 
    /*@ManDoc: destbox region of this FAB is numerator. srcbox regions of
               src FAB is denominator. srcbox and destbox must be same size.
    */
    BaseFab<T>& divide (const BaseFab<T>& src,
                        const Box&        srcbox,
                        const Box&        destbox,
                        int               srccomp,
                        int               destcomp,
                        int               numcomp=1);

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
    //
    // This function is called by the `BaseFab' setVal operations.
    //
    void performSetVal (T          x,
                        const Box& bx,
                        int        nstart,
                        int        numcomp);
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
int
BaseFab<T>::nComp () const
{
    return nvar;
}

template <class T>
const Box&
BaseFab<T>::box () const
{
    return domain;
}

template <class T>
const int*
BaseFab<T>::length () const
{
    return domain.length().getVect();
}

template <class T>
const IntVect&
BaseFab<T>::smallEnd () const
{
    return domain.smallEnd();
}

template <class T>
const IntVect&
BaseFab<T>::bigEnd () const
{
    return domain.bigEnd();
}

template <class T>
const int*
BaseFab<T>::loVect () const
{
    return domain.loVect();
}

template <class T>
const int*
BaseFab<T>::hiVect () const
{
    return domain.hiVect();
}

template <class T>
bool
BaseFab<T>::contains (const BaseFab<T>& fab) const
{
    return box().contains(fab.box()) && nvar <= fab.nvar;
}

template <class T>
bool
BaseFab<T>::contains (const Box& bx) const
{
    return box().contains(bx);
}

template <class T>
T*
BaseFab<T>::dataPtr (int n)
{
    BL_ASSERT(!(dptr == 0));
    return &dptr[n*numpts];
}

template <class T>
const T*
BaseFab<T>::dataPtr (int n) const
{
    BL_ASSERT(!(dptr == 0));
    return &dptr[n*numpts];
}

template <class T>
bool
BaseFab<T>::isAllocated () const
{
    return dptr != 0;
}

template <class T>
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
T&
BaseFab<T>::operator() (const IntVect& p)
{
    BL_ASSERT(!(dptr == 0));
    BL_ASSERT(domain.contains(p));

    return dptr[domain.index(p)];
}

template <class T>
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
    const long size    = domain.numPts();

    BL_ASSERT(!(dptr == 0));
    BL_ASSERT(n >= 0 && n + numcomp <= nvar);

    for (int k = 0; k < numcomp; k++)
        data[k] = dptr[loc+(n+k)*size];
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
    BL_ASSERT(!(The_FAB_Arena == 0));

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
    BL_ASSERT(!(The_FAB_Arena == 0));
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
BaseFab<T>::BaseFab ()
    :
    domain(Box()),
    nvar(0),
    numpts(0),
    truesize(0),
    dptr(0)
{}

template <class T>
BaseFab<T>::BaseFab (const BaseFab& fab)
    : domain(fab.box()),
      nvar(fab.nComp())
{
    numpts = domain.numPts();
    define();
    for ( int i = 0; i < numpts*nvar; ++i )
    {
	dptr[i] = fab.dptr[i];
    }
}

template <class T>
BaseFab<T>&
BaseFab<T>::operator= (const BaseFab& fab)
{
    if ( this == &fab ) return *this;
    domain = fab.box();
    nvar   = fab.nComp();
    numpts = domain.numPts();
    define();
    for ( int i = 0; i < numpts*nvar; ++i )
    {
	dptr[i] = fab.dptr[i];
    }
    return *this;
}

template <class T>
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
BaseFab<T>::~BaseFab ()
{
    undefine();
}

template <class T>
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
        const int *_subbox_len = _subbox_.length().getVect(); 
        const int *_bx_lo = (srcbox).loVect(); 
        // const int *_bx_len = (srcbox).length().getVect(); 
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
        const int *_subbox_len = _subbox_.length().getVect(); 
        const int *_bx_lo = (srcbox).loVect(); 
        // const int *_bx_len = (srcbox).length().getVect(); 
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
        const int *_subbox_len = _subbox_.length().getVect(); 
        const int *_bx_lo = (srcbox).loVect(); 
        // const int *_bx_len = (srcbox).length().getVect(); 
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

template <class T>
void
BaseFab<T>::performSetVal (T         val,
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
void
BaseFab<T>::setComplement (T          x,
                           const Box& b,
                           int        ns,
                           int        num)
{
    BoxList b_lst = BoxLib::boxDiff(domain,b);
    for (BoxList::Iterator bli = b_lst.begin(); bli != b_lst.end(); ++bli)
        performSetVal(x, *bli, ns, num);
}

template <class T>
void
BaseFab<T>::abs ()
{
    ForAllThis(Real)
    {
        thisR = std::abs(thisR);
    } EndFor
}

template <class T>
void
BaseFab<T>::abs (int comp,
                   int numcomp)
{
    ForAllThisNN(Real,comp,numcomp)
    {
        thisR = std::abs(thisR);
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
                    tmp[i] = std::abs(Real(row[i]));
            }
            else
            {
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = std::max(tmp[i],Real(std::abs(row[i])));
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
                    tmp[i] = std::abs(Real(row[i]));
            }
            else
            {
                for (int i = 0; i < thisLen; i++)
                    tmp[i] += std::abs(Real(row[i]));
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
            nrm += tmp[i];
    }
    else
    {
      BoxLib::Error("BaseFab::norm(): only p == 0 or p == 1 are supported");
    }

    delete [] tmp;

    return nrm;
}


#ifdef BL_ARCH_CRAY
#define RESTRICT restrict
//
// Template specialization for Real.
//
template<>
Real
BaseFab<Real>::norm (const Box& subbox,
                       int        p,
                       int        comp,
                       int        numcomp) const
{
    BL_ASSERT(comp >= 0 && comp+numcomp <= nComp());
    BL_ASSERT(p >= 0);

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
                tmplen = thisLen;
#pragma _CRI ivdep
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = fabs(row[i]);
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
                tmplen = thisLen;
#pragma _CRI ivdep
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = fabs(row[i]);
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
#endif

template <class T>
T
BaseFab<T>::min (int comp) const
{
    T*  _min_row = 0;
    int _X_len   = 0;
    ForAllThisCPencil(T,domain,comp,1)
    {
        const T* _row = &thisR;
        if (_min_row == 0)
        {
            _min_row = new T[thisLen];
            _X_len = thisLen;
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
    for (int i = 1; i < _X_len; i++)
        _min = std::min(_min,_min_row[i]);

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
            _X_len = thisLen;
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
    for (int i = 1; i < _X_len; i++)
        _min = std::min(_min,_min_row[i]);

    delete [] _min_row;

    return _min;
}

template <class T>
T
BaseFab<T>::max (int comp) const
{
    T* _max_row = 0;
    int _X_len  = 0;
    ForAllThisCPencil(T,domain,comp,1)
    {
        const T* _row = &thisR;
        if (_max_row== 0)
        {
            _max_row = new T[thisLen];
            _X_len = thisLen;
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
    for (int i = 1; i < _X_len; i++)
        _max = std::max(_max,_max_row[i]);

    delete [] _max_row;

    return _max;
}

template <class T>
T
BaseFab<T>::max (const Box& subbox,
                    int        comp) const
{
    T*  _max_row = 0;
    int _X_len   = 0;
    ForAllThisCPencil(T,subbox,comp,1)
    {
        const T* _row = &thisR;
        if (_max_row == 0)
        {
            _max_row = new T[thisLen];
            _X_len = thisLen;
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
    for (int i = 1; i < _X_len; i++)
        _max = std::max(_max,_max_row[i]);

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

#ifdef __GNUC__
typedef BaseFab<int> bfi;
#endif

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
        int ix = D_TERM(_i, +_j*_b_len[0], +_k*_b_len[0]*_b_len[1]);
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
        int ix = D_TERM(_i, +_j*_b_len[0], +_k*_b_len[0]*_b_len[1]);
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
        int ix = D_TERM(_i, +_j*_b_len[0], +_k*_b_len[0]*_b_len[1]);
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
        int ix = D_TERM(_i, +_j*_b_len[0], +_k*_b_len[0]*_b_len[1]);
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
        int ix = D_TERM(_i, +_j*_b_len[0], +_k*_b_len[0]*_b_len[1]);
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
        thisR = D_TERM(iR*100, +jR*10, + kR) + 1000*nR + 10000*mark;
    } EndFor
}

template <class T>
void
BaseFab<T>::copyRev (const Box&         destbox,
                      const BaseFab<T>& src,
                      const Box&         srcbox,
                      int                reversal_index,
                      T*                 multiplier)
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
    T*  _sum_row = 0;
    int _sum_len = 0;
    ForAllThisCPencil(T,domain,comp,numcomp)
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
BaseFab<T>::operator+= (T r)
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
BaseFab<T>::operator+= (const BaseFab<T>& x)
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
BaseFab<T>::operator-= (const BaseFab<T>& x)
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
BaseFab<T>::operator*= (T r)
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
BaseFab<T>::operator*= (const BaseFab<T> &x)
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
BaseFab<T>::operator/= (T r)
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
BaseFab<T>::operator/= (const BaseFab<T> &x)
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
                        const Box&         b1,
                        int                comp1,
                        const BaseFab<T>& f2,
                        const Box&         b2,
                        int                comp2,
                        Real               t1,
                        Real               t2,
                        Real               t,
                        const Box&         b,
                        int                comp,
                        int                numcomp)

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

#endif /*BL_BASEFAB_H*/
