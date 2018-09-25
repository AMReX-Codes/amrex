
#include <cstdlib>
#include <AMReX_Mask.H>
#include <AMReX_Utility.H>

namespace amrex {

Mask::Mask ()
    :
    BaseFab<int>() {}

Mask::Mask (const Box& bx,
            int        nc,
	    bool       alloc,
	    bool       shared)
    :
    BaseFab<int>(bx,nc,alloc,shared) {}

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

Mask&
Mask::And (const Mask& src)
{
    return this->And(src,domain,domain,0,0,nvar);
}

Mask&
Mask::And (const Mask& src,
           int         srccomp,
           int         destcomp,
           int         numcomp)
{
    return this->And(src,domain,domain,srccomp,destcomp,numcomp);
}

Mask&
Mask::And (const Mask& src,
           const Box&  subbox,
           int         srccomp,
           int         destcomp,
           int         numcomp)
{
    return this->And(src,subbox,subbox,srccomp,destcomp,numcomp);
}

Mask&
Mask::And (const Mask& src,
           const Box&  srcbox,
           const Box&  destbox,
           int         srccomp,
           int         destcomp,
           int         numcomp)
{
    const auto dp0 = this->stridedPtr(destbox,destcomp);
    const auto sp0 =   src.stridedPtr( srcbox, srccomp);
    const auto len = destbox.length3d();

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len[2]; ++k) {
            for     (int j = 0; j < len[1]; ++j) {
                value_type       * AMREX_RESTRICT dp = dp0(0,j,k,n);
                value_type const * AMREX_RESTRICT sp = sp0(0,j,k,n);
                for (int i = 0; i < len[0]; ++i) {
                    dp[i] = dp[i] ? sp[0] : 0;
                }
            }
        }
    }

    return *this;
}

Mask&
Mask::Or (const Mask& src)
{
    return this->Or(src,domain,domain,0,0,nvar);
}

Mask&
Mask::Or (const Mask& src,
          int         srccomp,
          int         destcomp,
          int         numcomp)
{
    return this->Or(src,domain,domain,srccomp,destcomp,numcomp);
}

Mask&
Mask::Or (const Mask& src,
          const Box&  subbox,
          int         srccomp,
          int         destcomp,
          int         numcomp)
{
    return this->Or(src,subbox,subbox,srccomp,destcomp,numcomp);
}

Mask&
Mask::Or (const Mask& src,
          const Box&  srcbox,
          const Box&  destbox,
          int         srccomp,
          int         destcomp,
          int         numcomp)
{
    const auto dp0 = this->stridedPtr(destbox,destcomp);
    const auto sp0 =   src.stridedPtr( srcbox, srccomp);
    const auto len = destbox.length3d();

    for (int n = 0; n < numcomp; ++n) {
        for         (int k = 0; k < len[2]; ++k) {
            for     (int j = 0; j < len[1]; ++j) {
                value_type       * AMREX_RESTRICT dp = dp0(0,j,k,n);
                value_type const * AMREX_RESTRICT sp = sp0(0,j,k,n);
                for (int i = 0; i < len[0]; ++i) {
                    dp[i] = dp[i] ? 1: sp[0];
                }
            }
        }
    }

    return *this;
}

}
