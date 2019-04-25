
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>

#include <AMReX_IArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FPC.H>

#include <AMReX_BLassert.H>
#include <AMReX.H>
#include <AMReX_Utility.H>

namespace amrex {

#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
bool IArrayBox::do_initval = true;
#else
bool IArrayBox::do_initval = false;
#endif

namespace
{
    bool initialized = false;
}

void
IArrayBox::Initialize ()
{
    if (initialized) return;
//    ParmParse pp("iab");
    amrex::ExecOnFinalize(IArrayBox::Finalize);
    initialized = true;
}

void
IArrayBox::Finalize ()
{
    initialized = false;
}

IArrayBox::IArrayBox () noexcept {}

IArrayBox::IArrayBox (const Box& b,
                      int        n,
		      bool       alloc,
		      bool       shared)
    :
    BaseFab<int>(b,n,alloc,shared)
{
    // For debugging purposes
    if ( alloc && do_initval ) {
	setVal(std::numeric_limits<int>::max());
    }
}

IArrayBox::IArrayBox (const IArrayBox& rhs, MakeType make_type, int scomp, int ncomp)
    :
    BaseFab<int>(rhs,make_type,scomp,ncomp)
{
}

IArrayBox&
IArrayBox::operator= (int v) noexcept
{
    BaseFab<int>::operator=(v);
    return *this;
}

void
IArrayBox::resize (const Box& b,
                   int        N)
{
    BaseFab<int>::resize(b,N);
    // For debugging purposes
    if ( do_initval ) {
        setVal(std::numeric_limits<int>::max());
    }
}

}
