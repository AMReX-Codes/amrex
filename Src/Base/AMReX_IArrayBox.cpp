
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

IArrayBox::IArrayBox (Arena* ar) noexcept
    : BaseFab<int>(ar)
{}

IArrayBox::IArrayBox (const Box& b, int n, Arena* ar)
    : BaseFab<int>(b,n,ar)
{
#ifndef AMREX_USE_GPU
    // For debugging purposes
    if ( do_initval ) {
	setVal<RunOn::Host>(std::numeric_limits<int>::max());
    }
#endif
}

IArrayBox::IArrayBox (const Box& b, int n, bool alloc, bool shared, Arena* ar)
    : BaseFab<int>(b,n,alloc,shared,ar)
{
#ifndef AMREX_USE_GPU
    // For debugging purposes
    if ( alloc && do_initval ) {
	setVal<RunOn::Host>(std::numeric_limits<int>::max());
    }
#endif
}

IArrayBox::IArrayBox (const IArrayBox& rhs, MakeType make_type, int scomp, int ncomp)
    : BaseFab<int>(rhs,make_type,scomp,ncomp)
{
}

void
IArrayBox::resize (const Box& b, int N)
{
    BaseFab<int>::resize(b,N);
    // For debugging purposes
    if ( do_initval ) {
#if defined(AMREX_USE_GPU)
        bool run_on_device = Gpu::inLaunchRegion() and
            (arena() == The_Arena() ||
             arena() == The_Device_Arena() ||
             arena() == The_Managed_Arena());
        if (run_on_device) {
            setVal<RunOn::Device>(std::numeric_limits<int>::max());
            Gpu::streamSynchronize();
        } else
#endif
        {
            setVal<RunOn::Host>(std::numeric_limits<int>::max());
        }
    }
}

}
