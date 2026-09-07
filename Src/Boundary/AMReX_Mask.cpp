
#include <AMReX_Mask.H>
#include <AMReX_Utility.H>
#include <cstdlib>

namespace amrex {

Mask::Mask (Arena* ar) noexcept
    : BaseFab<int>(ar) {}

Mask::Mask (const Box& bx, int nc, Arena* ar)
    : BaseFab<int>(bx,nc,ar) {}

Mask::Mask (const Box& bx, int nc, bool alloc, bool shared, Arena* ar)
    : BaseFab<int>(bx,nc,alloc,shared,ar) {}

Mask::Mask (Mask const& rhs, MakeType make_type, int scomp, int ncomp)
    :
    BaseFab<int>(rhs,make_type,scomp,ncomp) {}

Mask::Mask (std::istream& is)
{
    readFrom(is);
}

std::ostream&
operator<< (std::ostream& os,
            const Mask&   m)
{
    int ncomp = m.nComp();

    os << "(Mask: " << m.box() << " " << ncomp << "\n";

    const Mask* mp = &m;
#ifdef AMREX_USE_GPU
    Mask hostmask(The_Pinned_Arena());
    if (m.arena()->isManaged() || m.arena()->isDevice()) {
        hostmask.resize(m.box(), ncomp);
        Gpu::dtoh_memcpy_async(hostmask.dataPtr(), m.dataPtr(),
                               static_cast<std::size_t>(hostmask.size())*sizeof(int));
        Gpu::streamSynchronize();
        mp = &hostmask;
    }
#endif

    IntVect sm = m.box().smallEnd();
    IntVect bg = m.box().bigEnd();
    for (IntVect p = sm; p <= bg; m.box().next(p))
    {
        os << p;
        for (int k = 0; k < ncomp; k++) {
            os << "  " << (*mp)(p,k);
        }
        os << "\n";
    }
    os << ")\n";

    BL_ASSERT(os.good());

    return os;
}

std::istream&
operator>> (std::istream& is,
            Mask&         m)
{
    is.ignore(BL_IGNORE_MAX,':');
    Box b;
    int ncomp;
    is >> b >> ncomp;
    AMREX_ASSERT(ncomp >= 0 && ncomp < std::numeric_limits<int>::max());
    is.ignore(BL_IGNORE_MAX, '\n');
    m.resize(b,ncomp);

    Mask* mp = &m;
#ifdef AMREX_USE_GPU
    Mask hostmask(The_Pinned_Arena());
    if (m.arena()->isManaged() || m.arena()->isDevice()) {
        hostmask.resize(b,ncomp);
        mp = &hostmask;
    }
#endif

    IntVect sm = b.smallEnd();
    IntVect bg = b.bigEnd();
    IntVect q;
    for (IntVect p = sm; p <= bg; b.next(p))
    {
        is >> q;
        BL_ASSERT( p == q);
        for( int k=0; k<ncomp; k++ ) { is >> (*mp)(p,k); }
        is.ignore(BL_IGNORE_MAX, '\n');
    }
    is.ignore(BL_IGNORE_MAX,'\n');

#ifdef AMREX_USE_GPU
    if (mp != &m) {
        Gpu::htod_memcpy_async(m.dataPtr(), hostmask.dataPtr(),
                               static_cast<std::size_t>(hostmask.size())*sizeof(int));
        Gpu::streamSynchronize();
    }
#endif

    BL_ASSERT(is.good());
    return is;
}

void
Mask::writeOn (std::ostream& os) const
{
    os << "(Mask: " << domain << " " << nvar << "\n";
    const int* ptr = dataPtr();
    auto len = static_cast<std::size_t>(domain.numPts()) * static_cast<std::size_t>(nvar);
#ifdef AMREX_USE_GPU
    Mask hostmask(The_Pinned_Arena());
    if (this->arena()->isManaged() || this->arena()->isDevice()) {
        hostmask.resize(domain, nvar);
        Gpu::dtoh_memcpy_async(hostmask.dataPtr(), ptr, len*sizeof(int));
        Gpu::streamSynchronize();
        ptr = hostmask.dataPtr();
    }
#endif
    os.write(reinterpret_cast<char const*>(ptr),
             static_cast<std::streamsize>(len*sizeof(int)));
    os << ")\n";
}

void
Mask::readFrom (std::istream& is)
{
    is.ignore(BL_IGNORE_MAX,':');
    Box b;
    int ncomp;
    is >> b >> ncomp;
    AMREX_ASSERT(ncomp >= 0 && ncomp < std::numeric_limits<int>::max());
    is.ignore(BL_IGNORE_MAX, '\n');
    resize(b,ncomp);
    int *ptr = dataPtr();
    auto len = static_cast<std::size_t>(domain.numPts()) * static_cast<std::size_t>(nvar);
#ifdef AMREX_USE_GPU
    Mask hostmask(The_Pinned_Arena());
    if (this->arena()->isManaged() || this->arena()->isDevice()) {
        hostmask.resize(domain, nvar);
        ptr = hostmask.dataPtr();
    }
#endif
    is.read(reinterpret_cast<char*>(ptr),
            static_cast<std::streamsize>(len*sizeof(int)));
#ifdef AMREX_USE_GPU
    if (ptr != dataPtr()) {
        Gpu::htod_memcpy_async(dataPtr(), ptr, len*sizeof(int));
        Gpu::streamSynchronize();
    }
#endif
    is.ignore(BL_IGNORE_MAX, '\n');
}

}
