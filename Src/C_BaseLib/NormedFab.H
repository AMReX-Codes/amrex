#ifndef BL_NORMEDFAB_H
#define BL_NORMEDFAB_H

//
// $Id: NormedFab.H,v 1.12 2001-07-19 20:02:46 lijewski Exp $
//

#include <algorithm>
#include <iostream>
#include <cmath>

#include <BoxLib.H>
#include <Looping.H>
#include <BaseFab.H>
#include <REAL.H>

//
//@Man:
//@Memo: A BaseFab with Norms
/*@Doc:

  This is NOT a polymorphic class.

  It provides various L-p norms to a BaseFab.

  It does NOT provide a copy constructor or assignment operator.
 */

template <class T>
class NormedFab
    :
    public virtual BaseFab<T>
{
public:
    //
    //@ManDoc: The default constructor.
    //
    NormedFab ();
    //
    //@ManDoc: Construct a NormedFab from a Box.
    //
    explicit NormedFab (const Box& b,
                        int        n = 1);
    //
    //@ManDoc: The destructor.
    //
    ~NormedFab ();

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
private:
    //
    // These are not allowed.
    //
    NormedFab (const NormedFab<T>&);
    NormedFab<T>& operator= (const NormedFab<T>&);
};

template <class T>
NormedFab<T>::NormedFab ()
{}

template <class T>
NormedFab<T>::NormedFab (const Box& b,
                         int        n)
    :
    BaseFab<T>(b,n)
{}

template <class T>
NormedFab<T>::~NormedFab()
{}

template <class T>
void
NormedFab<T>::abs ()
{
    ForAllThis(Real)
    {
        thisR = std::abs(thisR);
    } EndFor
}

template <class T>
void
NormedFab<T>::abs (int comp,
                   int numcomp)
{
    ForAllThisNN(Real,comp,numcomp)
    {
        thisR = std::abs(thisR);
    } EndFor
}

template <class T>
void
NormedFab<T>::abs (const Box& subbox,
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
NormedFab<T>::norm (int p,
                    int comp,
                    int numcomp) const
{
    return norm(domain,p,comp,numcomp);
}

template <class T>
Real
NormedFab<T>::norm (const Box& subbox,
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
      BoxLib::Error("NormedFab::norm(): only p == 0 or p == 1 are supported");
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
NormedFab<Real>::norm (const Box& subbox,
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
      BoxLib::Error("NormedFab<Real>::norm(): only p == 0 or p == 1 are supported");

    delete [] tmp;

    return nrm;
}
#endif

#endif /*BL_NORMEDFAB_H*/
