//BL_COPYRIGHT_NOTICE

#ifndef BL_ARITHFAB_H
#define BL_ARITHFAB_H

//
// $Id: ArithFab.H,v 1.3 1999-04-15 23:27:58 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <iostream>
#include <iomanip>
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
#else
#include <iostream.h>
#include <iomanip.h>
#endif

#include <BaseFab.H>
#include <Misc.H>
#include <REAL.H>
#include <Looping.H>

//
//@Man:
//@Memo: A BaseFab<T> Specialized for Numerical Operations
/*@Doc: 

  An ArithFab<T> elaborates the BaseFab concept by adding numerical 
  operations that make sense for arithmetic types T, e.g. integers
  and reals.  For such types we define Fab-wide sums, inverses,
  negations, additions, subtractions, multiplications, divisions,
  linear interpolations, and linear combinations.

  The numerical operations are generally overloaded, allowing a choice
  of a very general operation to a easy-to-use very specific operation.
  Note that many operations offer the same choice of overloaded
  signatures.

  Note that some functions are denoted as scalar operations while other
  functions are denoted as Fab operations.  The difference is that the
  scalar operations take a scalar as one argument, while FAB operations
  utilize two FABS.

  This is NOT a polymorphic class.

  It does NOT provide a copy constructor or assignment operator.
*/

template <class T>
class ArithFab
    :
    public virtual BaseFab<T>
{
public:

    /*@ManDoc: Default constructor uses BaseFab default constructor.
               It doesn't allocate any memory.
    */
    ArithFab ();
    //
    //@ManDoc: Allocate ArithFab for given box and number of components.
    //
    explicit ArithFab (const Box& b,
                       int        n = 1);
    //
    //@ManDoc: The destructor.
    //
    ~ArithFab ();
    //
    //@ManDoc: Fill with a pattern of numbers.
    //
    void patternFill (int mark = 0);

    /*@ManDoc: Copies with index reversal from srcbox region of src
               into destbox region of this FAB.  All components are
               copied and multiplied by corresponding multiplier.
    */
    void copyRev (const Box&         destbox,
                  const ArithFab<T>& src,
                  const Box&         srcbox,
                  int                reversal_index,
                  T*                 multiplier);
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
    ArithFab<T>& invert (T          v,
                         const Box& subbox,
                         int        comp=0,
                         int        numcomp=1);
    //
    //@ManDoc: As above except on entire domain.
    //
    ArithFab<T>& invert (T   v,
                         int comp,
                         int numcomp=1);
    //
    //@ManDoc: As above except on entire domain, all components.
    //
    ArithFab<T>& invert (T v);
    //
    //@ManDoc: Negate ArithFab, most general.
    //
    ArithFab<T>& negate (const Box& subbox,
                         int        comp=0,
                         int        numcomp=1);
    //
    //@ManDoc: As above, except on entire domain.
    //
    ArithFab<T>& negate (int comp,
                         int numcomp=1);
    //
    //@ManDoc: As above, except on entire domain and all components.
    //
    ArithFab<T>& negate ();
    //
    //@ManDoc: Scalar addition (a[i] <- a[i] + r), most general.
    //
    ArithFab<T>& plus (T          r,
                       const Box& b,
                       int        comp=0,
                       int        numcomp=1);
    //
    //@ManDoc: As above, except on entire domain.
    //
    ArithFab<T>& plus (T   r,
                       int comp,
                       int numcomp=1);
    //
    //@ManDoc: As above, except on entire domain and all components.
    //
    ArithFab<T>& plus (T r);

    /*@ManDoc: Addition in place.  This will often be more efficient than
               making new ArithFab for result.
    */
    ArithFab<T>& operator+= (T r);
    //
    //@ManDoc: FAB addition (a[i] <- a[i] + b[i]) in place.
    //
    ArithFab<T>& operator+= (const ArithFab<T>& f);
    //
    //@ManDoc: FAB addition (a[i] <- a[i] + b[i]). The same as += operator.
    //
    ArithFab<T>& plus (const ArithFab<T>& src);

    /*@ManDoc: Add src components (srccomp:srccomp+numcomp-1) to
               this FAB's components (destcomp:destcomp+numcomp-1)
               where the two FABs intersect.
    */
    ArithFab<T>& plus (const ArithFab<T>& src,
                       int                srccomp,
                       int                destcomp,
                       int                numcomp=1);

    /*@ManDoc: Same as above except addition is restricted to intersection
               of subbox and src FAB. NOTE: subbox must be contained in this
               FAB.
    */
    ArithFab<T>& plus (const ArithFab<T>& src,
                       const Box&         subbox,
                       int                srccomp,
                       int                destcomp,
                       int                numcomp=1);

    /*@ManDoc: Add srcbox region of src FAB to destbox region of this FAB.
               The srcbox and destbox must be same size.
    */
    ArithFab<T>& plus (const ArithFab<T>& src,
                       const Box&         srcbox,
                       const Box&         destbox,
                       int                srccomp,
                       int                destcomp,
                       int                numcomp=1);

    /*@ManDoc: Scalar subtraction (a[i] <- a[i] - r).
               Note: use plus(-r) for more general operations.
    */
    ArithFab<T>& operator-= (T r);
    //
    //@ManDoc: FAB subtraction (a[i] <- a[i] - b[i]), in place.
    //
    ArithFab<T>& operator-= (const ArithFab<T>& f);
    //
    //@ManDoc: FAB subtraction (a[i] <- a[i] - b[i]). The same as -= operator.
    //
    ArithFab<T>& minus (const ArithFab<T>& src);

    /*@ManDoc: Subtract src components (srccomp:srccomp+numcomp-1) to
               this FAB's components (destcomp:destcomp+numcomp-1) where
               the two FABs intersect.
    */
    ArithFab<T>& minus (const ArithFab<T>& src,
                        int                srccomp,
                        int                destcomp,
                        int                numcomp=1);

    /*@ManDoc: Same as above except subtraction is restricted to intersection
               of subbox and src FAB.  NOTE: subbox must be contained in
               this FAB.
    */
    ArithFab<T>& minus (const ArithFab<T>& src,
                        const Box&         subbox,
                        int                srccomp,
                        int                destcomp,
                        int                numcomp=1);

    /*@ManDoc: Subtract srcbox region of src FAB from destbox region
               of this FAB. srcbox and destbox must be same size.
    */
    ArithFab<T>& minus (const ArithFab<T>& src,
                        const Box&         srcbox,
                        const Box&         destbox,
                        int                srccomp,
                        int                destcomp,
                        int                numcomp=1);
    //
    //@ManDoc: Scalar multiplication (a[i] <- a[i] * r), in place.
    //
    ArithFab<T>& operator*= (T r);
    //
    //@ManDoc: Scalar multiplication (a[i] <- a[i] * r).  The same as *=.
    //
    ArithFab<T>& mult (T r);

    /*@ManDoc: Scalar multiplication, except control which components
               are multiplied.
    */
    ArithFab<T>& mult (T   r,
                       int comp,
                       int numcomp=1);
    //
    //@ManDoc: As above, except specify sub-box.
    //
    ArithFab<T>& mult (T          r,
                       const Box& b,
                       int        comp=0,
                       int        numcomp=1);
    //
    //@ManDoc: FAB multiplication (a[i] <- a[i] * b[i]), in place.
    //
    ArithFab<T>& operator*= (const ArithFab<T>& f);
    //
    //@ManDoc: As above.
    //
    ArithFab<T>& mult (const ArithFab<T>& src);

    /*@ManDoc: Multiply src components (srccomp:srccomp+numcomp-1) with
               this FAB's components (destcomp:destcomp+numcomp-1) where
               the two FABs intersect.
    */
    ArithFab<T>& mult (const ArithFab<T>& src,
                       int                srccomp,
                       int                destcomp,
                       int                numcomp=1);

    /*@ManDoc: Same as above except multiplication is restricted to
               intersection of subbox and src FAB.  NOTE: subbox must be
               contained in this FAB.
    */
    ArithFab<T>& mult (const ArithFab<T>& src,
                       const Box&         subbox,
                       int                srccomp,
                       int                destcomp,
                       int                numcomp=1);

    /*@ManDoc: Multiply srcbox region of src FAB with destbox region
               of this FAB. The srcbox and destbox must be same size.
    */
    ArithFab<T>& mult (const ArithFab<T>& src,
                       const Box&         srcbox,
                       const Box&         destbox,
                       int                srccomp,
                       int                destcomp,
                       int                numcomp=1);
    //
    //@ManDoc: Scalar division (a[i] <- a[i] / r), in place.
    //
    ArithFab<T>& operator/= (T r);
    //
    //@ManDoc: Scalar division (a[i] <- a[i] / r), in place.
    //
    ArithFab<T>& divide (T r);
    //
    //@ManDoc: As above except specify which components.
    //
    ArithFab<T>& divide (T   r,
                         int comp,
                         int numcomp=1);
    //
    //@ManDoc: As above except specify sub-box.
    //
    ArithFab<T>& divide (T          r,
                         const Box& b,
                         int        comp=0,
                         int        numcomp=1);
    //
    //@ManDoc: FAB division, in place.
    //
    ArithFab<T>& operator/= (const ArithFab<T>& src);
    //
    //@ManDoc: Same as above.
    //
    ArithFab<T>& divide (const ArithFab<T>& src);

    /*@ManDoc: This FAB is numerator, src FAB is denominator
               divide src components (srccomp:srccomp+numcomp-1) into
               this FAB's components (destcomp:destcomp+numcomp-1)
               where the two FABs intersect.
    */
    ArithFab<T>& divide (const ArithFab<T>& src,
                         int                srccomp,
                         int                destcomp,
                         int                numcomp=1);

    /*@ManDoc: Same as above except division is restricted to
               intersection of subbox and src FAB.  NOTE: subbox must be
               contained in this FAB.
    */
    ArithFab<T>& divide (const ArithFab<T>& src,
                         const Box&         subbox,
                         int                srccomp,
                         int                destcomp,
                         int                numcomp=1);
 
    /*@ManDoc: destbox region of this FAB is numerator. srcbox regions of
               src FAB is denominator. srcbox and destbox must be same size.
    */
    ArithFab<T>& divide (const ArithFab<T>& src,
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
    ArithFab<T>& linInterp (const ArithFab<T>& f1,
                            const Box&         b1,
                            int                comp1,
                            const ArithFab<T>& f2,
                            const Box&         b2,
                            int                comp2,
                            Real               t1,
                            Real               t2,
                            Real               t,
                            const Box&         b,
                            int                comp,
                            int                numcomp = 1);

    /*@ManDoc: Linear combination.  Result is alpha*f1 + beta*f2.
               Data is taken from b1 region of f1, b2 region of f2
               and stored in b region of this FAB.
               Boxes b, b1 and b2 must be the same size.
               Data is taken from component comp1 of f1, comp2 of f2,
               and stored in component comp of this FAB.
               This FAB is returned as a reference for chaining.
    */
    ArithFab<T>& linComb (const ArithFab<T>& f1,
                          const Box&         b1,
                          int                comp1,
                          const ArithFab<T>& f2,
                          const Box&         b2,
                          int                comp2,
                          Real               alpha,
                          Real               beta,
                          const Box&         b,
                          int                comp,
                          int                numcomp = 1);

private:
    //
    // Disallowed.
    //
    ArithFab (const ArithFab<T>&);
    ArithFab<T>& operator= (const ArithFab<T>&);
};

//
// Implementation.
//

template <class T>
inline
ArithFab<T>&
ArithFab<T>::plus (T r)
{
    return operator+=(r);
}

template <class T>
inline
ArithFab<T>&
ArithFab<T>::plus (const ArithFab<T>& x)
{
    return operator+=(x);
}

template <class T>
inline
ArithFab<T>&
ArithFab<T>::operator-= (T r)
{
    return operator+=(-r);
}

template <class T>
inline
ArithFab<T>&
ArithFab<T>::minus (const ArithFab<T>& x)
{
    return operator-=(x);
}


template <class T>
inline
ArithFab<T>&
ArithFab<T>::mult (T r)
{
    return operator*=(r);
}

template <class T>
inline
ArithFab<T>&
ArithFab<T>::mult (const ArithFab<T>& x)
{
    return operator*=(x);
}

template <class T>
inline
ArithFab<T>&
ArithFab<T>::divide (T r)
{
    return operator/=(r);
}

template <class T>
inline
ArithFab<T>&
ArithFab<T>::divide (const ArithFab<T>& x)
{
    return operator/=(x);
}

template <class T>
inline
ArithFab<T>::ArithFab ()
{}

template <class T>
inline
ArithFab<T>::ArithFab (const Box &b,
                       int        n)
    :
    BaseFab<T>(b,n)
{}

template <class T>
inline
ArithFab<T>::~ArithFab ()
{}

template <class T>
void
ArithFab<T>::patternFill (int mark)
{
    ForAllThis(T)
    {
        thisR = D_TERM(iR*100, +jR*10, + kR) + 1000*nR + 10000*mark;
    } EndFor
}

template <class T>
void
ArithFab<T>::copyRev (const Box&         destbox,
                      const ArithFab<T>& src,
                      const Box&         srcbox,
                      int                reversal_index,
                      T*                 multiplier)
{
    ArithFab<T>& dest = *this;

    ForAllRevXBNYCBNNN(T,dest,destbox,0,src,srcbox,0,nComp(),reversal_index)
    {
        destR = multiplier[_n]*srcR;
    } EndFor
}

#if !defined(BL_CRAY_BUG_DEFARG)
//
// Yes, a global function.
//
template <class T>
void
copyFABRev (const ArithFab<T>& src,
            const Box&         srcbox,
            ArithFab<T>&       dest,
            const Box&         destbox,
            int                reversal_index,
            int                comp,
            T                  multiplier)
{
    ForAllRevXBNYCBNNN(T,dest,destbox,comp,src,srcbox,comp,1,reversal_index)
    {
        destR = multiplier * srcR;
    } EndFor
}
#endif

template <class T>
T
ArithFab<T>::sum (int comp,
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
ArithFab<T>::sum (const Box& subbox,
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
ArithFab<T>&
ArithFab<T>::negate ()
{
    ForAllThis(T)
    {
        thisR = - thisR;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::negate (int comp,
                     int numcomp)
{
    ForAllThisNN(T,comp,numcomp)
    {
        thisR = - thisR;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::negate (const Box& b,
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
ArithFab<T>&
ArithFab<T>::invert (T r)
{
    ForAllThis(T)
    {
        thisR = r/thisR;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::invert (T   r,
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
ArithFab<T>&
ArithFab<T>::invert (T          r,
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
ArithFab<T>&
ArithFab<T>::operator+= (T r)
{
    ForAllThis(T)
    {
        thisR += r;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::plus (T   r,
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
ArithFab<T>&
ArithFab<T>::plus (T          r,
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
ArithFab<T>&
ArithFab<T>::operator+= (const ArithFab<T>& x)
{
    ForAllThisXC(T,x)
    {
        thisR += xR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::plus (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::plus (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::plus (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::operator-= (const ArithFab<T>& x)
{
    ForAllThisXC(T,x)
    {
        thisR -= xR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::minus (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::minus (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::minus (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::operator*= (T r)
{
    ForAllThis(T)
    {
        thisR *= r;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::mult (T   r,
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
ArithFab<T>&
ArithFab<T>::mult (T          r,
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
ArithFab<T>&
ArithFab<T>::operator*= (const ArithFab<T> &x)
{
    ForAllThisXC(T,x)
    {
        thisR *= xR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::mult (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::mult (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::mult (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::operator/= (T r)
{
    ForAllThis(T)
    {
        thisR /= r;
    } EndFor
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::divide (T   r,
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
ArithFab<T>&
ArithFab<T>::divide (T          r,
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
ArithFab<T>&
ArithFab<T>::operator/= (const ArithFab<T> &x)
{
    ForAllThisXC(T,x)
    {
        thisR /= xR;
    } EndForTX
    return *this;
}

template <class T>
ArithFab<T>&
ArithFab<T>::divide (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::divide (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::divide (const ArithFab<T>& src,
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
ArithFab<T>&
ArithFab<T>::linInterp (const ArithFab<T>& f1,
                        const Box&         b1,
                        int                comp1,
                        const ArithFab<T>& f2,
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
ArithFab<T>&
ArithFab<T>::linComb (const ArithFab<T>& f1,
                      const Box&         b1,
                      int                comp1,
                      const ArithFab<T>& f2,
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

#endif /*BL_ARITHFAB_H*/
