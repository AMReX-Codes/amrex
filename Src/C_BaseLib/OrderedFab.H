//BL_COPYRIGHT_NOTICE

#ifndef BL_ORDEREDFAB_H
#define BL_ORDEREDFAB_H

//
// $Id: OrderedFab.H,v 1.1 1997-12-14 01:05:54 lijewski Exp $
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
#include <Looping.H>

//
//@Man:
//@Memo: A BaseFab Specialized on Ordering between Elements of Type T
/*@Doc: 

  An OrderedFab<T> elaborates the BaseFab<T> concept by adding 
  operations that require an ordering between elements of type T.
  Specifically, an OrderedFab allows the computation of min's, max's,
  minIndex, maxIndex, and the setting of masks.

  The operations are generally overloaded, allowing a choice
  of a very general operation to a easy-to-use very specific operation.
  Note that many operations offer the same choice of overloaded
  signatures.

  This class does NOT provide a copy constructor or assignment operator.
*/

template <class T>
class OrderedFab
    :
    public virtual BaseFab<T>
{
public:

    /*@ManDoc: Default constructor uses BaseFab default constructor;
               i.e. allocates no memory.
    */
    OrderedFab ();
    //
    //@ManDoc: Allocate OrderedFab for given box and number of components.
    //
    explicit OrderedFab (const Box& b,
                         int        n = 1);
    //
    //@ManDoc: The (virtual) destructor.
    //
    virtual ~OrderedFab ();
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
               OrderedFab as value less than val, 0 otherwise.
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

private:
    //
    // These are disallowed.
    //
    OrderedFab (const OrderedFab<T>&);
    OrderedFab<T>& operator= (const OrderedFab<T>&);
};

//
// Implementation.
//

template <class T>
inline
OrderedFab<T>::OrderedFab ()
    : BaseFab<T>()
{}

template <class T>
inline
OrderedFab<T>::OrderedFab (const Box& b,
                           int        n)
    : BaseFab<T>(b,n)
{}

//
// This isn't inlined as it's virtual.
//

template <class T>
OrderedFab<T>::~OrderedFab ()
{}

template <class T>
T
OrderedFab<T>::min (int comp) const
{
    T *_min_row = 0;
    int _X_len = 0;
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
OrderedFab<T>::min (const Box& subbox,
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
OrderedFab<T>::max (int comp) const
{
    T *_max_row = 0;
    int _X_len = 0;
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
OrderedFab<T>::max (const Box& subbox,
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
OrderedFab<T>::minIndex (int comp) const
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
OrderedFab<T>::minIndex (const Box& subbox,
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
OrderedFab<T>::maxIndex (int comp) const
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
OrderedFab<T>::maxIndex (const Box& subbox,
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
OrderedFab<T>::maskLT (BaseFab<int>& mask,
                       T             val,
                       int           comp) const
{
    mask.resize(domain,1);
    mask.setVal(0);

    int *mptr = mask.dataPtr();
    int cnt = 0;

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
OrderedFab<T>::maskLE (BaseFab<int>& mask,
                       T             val,
                       int           comp) const
{
    mask.resize(domain,1);
    mask.setVal(0);

    int *mptr = mask.dataPtr();
    int cnt = 0;

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
OrderedFab<T>::maskEQ (BaseFab<int>& mask,
                       T             val,
                       int           comp) const
{
    mask.resize(domain,1);
    mask.setVal(0);

    int *mptr = mask.dataPtr();
    int cnt = 0;

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
OrderedFab<T>::maskGT (BaseFab<int>& mask,
                       T             val,
                       int           comp) const
{
    mask.resize(domain,1);
    mask.setVal(0);

    int *mptr = mask.dataPtr();
    int cnt = 0;

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
OrderedFab<T>::maskGE(BaseFab<int>& mask,
                      T             val,
                      int           comp) const
{
    mask.resize(domain,1);
    mask.setVal(0);

    int *mptr = mask.dataPtr();
    int cnt = 0;

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

#endif /*BL_ORDEREDFAB_H*/
