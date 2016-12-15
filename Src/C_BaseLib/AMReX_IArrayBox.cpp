
#include <winstd.H>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>

#include <IArrayBox.H>
#include <ParmParse.H>
#include <FPC.H>

#include <BLassert.H>
#include <BoxLib.H>
#include <Looping.H>
#include <Utility.H>

namespace
{
    bool initialized = false;
}

IArrayBox::IArrayBox () {}

IArrayBox::IArrayBox (const Box& b,
                      int        n,
		      bool       alloc,
		      bool       shared)
    :
    BaseFab<int>(b,n,alloc,shared)
{
    //
    // For debugging purposes set values to QNAN when possible.
    //
    if ( alloc && do_initval )
	setVal(0);
}

IArrayBox&
IArrayBox::operator= (const int& v)
{
    BaseFab<int>::operator=(v);
    return *this;
}

void
IArrayBox::resize (const Box& b,
                   int        N)
{
    BaseFab<int>::resize(b,N);
    //
    // For debugging purposes set values to QNAN when possible.
    //
    if ( do_initval )
        setVal(0);
}

int
IArrayBox::norm (int p,
                 int comp,
                 int numcomp) const
{
    return norm(domain,p,comp,numcomp);
}

IArrayBox::~IArrayBox () {}

#if !defined(NDEBUG)
bool IArrayBox::do_initval = true;
#else
bool IArrayBox::do_initval = false;
#endif

void
IArrayBox::Initialize ()
{
    if (initialized) return;
//    ParmParse pp("iab");
    BoxLib::ExecOnFinalize(IArrayBox::Finalize);
    initialized = true;
}

void
IArrayBox::Finalize ()
{
    initialized = false;
}

int
IArrayBox::norm (const Box& subbox,
                 int        p,
                 int        comp,
                 int        ncomp) const
{
    BL_ASSERT(p >= 0);
    BL_ASSERT(comp >= 0 && comp+ncomp <= nComp());

    int  nrm    = 0;
    int* tmp    = 0;
    int   tmplen = 0;

    if (p == 0 || p == 1)
    {
        nrm = BaseFab<int>::norm(subbox,p,comp,ncomp);
    }
    else if (p == 2)
    {
        ForAllThisCPencil(int,subbox,comp,ncomp)
        {
            const int* row = &thisR;
            if (tmp == 0)
            {
                tmp    = new int[thisLen];
                tmplen = thisLen;
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = row[i]*row[i];
            }
            else
            {
                for (int i = 0; i < thisLen; i++)
                    tmp[i] += row[i]*row[i];
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
            nrm += tmp[i];
        nrm = std::sqrt(double(nrm));
    }
    else
    {
        int pwr = int(p);
        ForAllThisCPencil(int,subbox,comp,ncomp)
        {
            const int* row = &thisR;
            if (tmp == 0)
            {
                tmp = new int[thisLen];
                tmplen = thisLen;
                for (int i = 0; i < thisLen; i++)
                    tmp[i] = std::pow((double)row[i],pwr);
            }
            else
            {
                for (int i = 0; i < thisLen; i++)
                    tmp[i] += std::pow((double)row[i],pwr);
            }
        } EndForPencil
        nrm = tmp[0];
        for (int i = 1; i < tmplen; i++)
            nrm += tmp[i];
        int invpwr = 1.0/pwr;
        nrm = std::pow((double)nrm,invpwr);
    }

    delete [] tmp;

    return nrm;
}
