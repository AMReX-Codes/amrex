
#include <cstdlib>
#include <AMReX_Mask.H>
#include <AMReX_Utility.H>

namespace amrex {

Mask::Mask () noexcept
    :
    BaseFab<int>() {}

Mask::Mask (const Box& bx,
            int        nc,
	    bool       alloc,
	    bool       shared)
    :
    BaseFab<int>(bx,nc,alloc,shared) {}

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

    IntVect sm = m.box().smallEnd();
    IntVect bg = m.box().bigEnd();
    for (IntVect p = sm; p <= bg; m.box().next(p))
    {
        os << p;
        for (int k = 0; k < ncomp; k++)
            os << "  " << m(p,k);
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
    is.ignore(BL_IGNORE_MAX, '\n');
    m.resize(b,ncomp);
    IntVect sm = b.smallEnd();
    IntVect bg = b.bigEnd();
    IntVect q;
    for (IntVect p = sm; p <= bg; b.next(p))
    {
        is >> q;
        BL_ASSERT( p == q);
        for( int k=0; k<ncomp; k++ ) is >> m(p,k);
        is.ignore(BL_IGNORE_MAX, '\n');
    }
    is.ignore(BL_IGNORE_MAX,'\n');
    BL_ASSERT(is.good());
    return is;
}

void
Mask::writeOn (std::ostream& os) const
{
    os << "(Mask: " << domain << " " << nvar << "\n";
    const int* ptr = dataPtr();
    int len = domain.numPts();
    os.write( (char*) ptr, len*sizeof(int) );
    os << ")\n";
}

void
Mask::readFrom (std::istream& is)
{
    is.ignore(BL_IGNORE_MAX,':');
    Box b;
    int ncomp;
    is >> b >> ncomp;
    is.ignore(BL_IGNORE_MAX, '\n');
    resize(b,ncomp);
    int *ptr = dataPtr();
    int len = domain.numPts();
    is.read( (char*) ptr, len*sizeof(int) );
    is.ignore(BL_IGNORE_MAX, '\n');
}

}
