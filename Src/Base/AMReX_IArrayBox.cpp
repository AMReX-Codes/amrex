#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_VectorIO.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>

namespace amrex {

#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
bool IArrayBox::do_initval = true;
#else
bool IArrayBox::do_initval = false;
#endif

std::unique_ptr<IFABio> IArrayBox::ifabio;

namespace
{
    bool initialized = false;
}

void
IArrayBox::Initialize ()
{
    if (initialized) return;
//    ParmParse pp("iab");

    ifabio = std::make_unique<IFABio>();

    amrex::ExecOnFinalize(IArrayBox::Finalize);
    initialized = true;
}

void
IArrayBox::Finalize ()
{
    ifabio.reset();
    initialized = false;
}

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
IArrayBox::resize (const Box& b, int N, Arena* ar)
{
    BaseFab<int>::resize(b,N,ar);
    // For debugging purposes
    if ( do_initval ) {
#if defined(AMREX_USE_GPU)
        if (Gpu::inLaunchRegion() && arena()->isDeviceAccessible()) {
            setVal<RunOn::Device>(std::numeric_limits<int>::max());
            Gpu::streamSynchronize();
        } else
#endif
        {
            setVal<RunOn::Host>(std::numeric_limits<int>::max());
        }
    }
}

std::unique_ptr<IntDescriptor>
IArrayBox::getDataDescriptor ()
{
    return std::make_unique<IntDescriptor>(FPC::NativeIntDescriptor());
}

std::string
IArrayBox::getClassName ()
{
    return std::string("amrex::IArrayBox");
}

IFABio const&
IArrayBox::getFABio ()
{
    return *ifabio;
}

void
IArrayBox::readFrom (std::istream& is)
{
    std::string type;
    is >> type;
    if (type != "IFAB") {
        amrex::Error(std::string("IArrayBox::readFrom: IFAB is expected, but instead we have ")
                     +type);
    }

    IntDescriptor data_descriptor;
    is >> data_descriptor;

    Box tmp_box;
    int tmp_ncomp;
    is >> tmp_box;
    is >> tmp_ncomp;
    is.ignore(99999, '\n');

    if (this->box() != tmp_box || this->nComp() != tmp_ncomp) {
        this->resize(tmp_box, tmp_ncomp);
    }

#ifdef AMREX_USE_GPU
    if (this->arena()->isManaged() || this->arena()->isDevice()) {
        IArrayBox hostfab(this->box(), this->nComp(), The_Pinned_Arena());
        ifabio->read(is, hostfab, data_descriptor);
        Gpu::htod_memcpy_async(this->dataPtr(), hostfab.dataPtr(),
                               hostfab.size()*sizeof(IArrayBox::value_type));
        Gpu::streamSynchronize();
    } else
#endif
    {
        ifabio->read(is, *this, data_descriptor);
    }
}

void
IFABio::write_header (std::ostream& os, const IArrayBox& fab, int nvar)
{
    AMREX_ASSERT(nvar <= fab.nComp());
    os <<"IFAB " << FPC::NativeIntDescriptor();
    os << fab.box() << ' ' << nvar << '\n';
}

void
IFABio::read (std::istream& is, IArrayBox& fab, IntDescriptor const& data_descriptor)
{
    readIntData(fab.dataPtr(), fab.size(), is, data_descriptor);
}

}
