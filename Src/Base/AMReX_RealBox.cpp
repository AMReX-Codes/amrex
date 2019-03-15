
#include <iostream>
#include <string>

#include <AMReX_RealBox.H>

namespace amrex {

RealBox::RealBox (const Box&  bx,
                  const Real* dx,
                  const Real* base) noexcept
{
    const int* blo = bx.loVect();
    const int* bhi = bx.hiVect();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        xlo[i] = base[i] + dx[i]*blo[i];
        int shft = (bx.type(i) == IndexType::CELL ? 1 : 0);
        xhi[i] = base[i] + dx[i]*(bhi[i]+ shft);
    }   
}

RealBox::RealBox (const std::array<Real,AMREX_SPACEDIM>& a_lo,
                  const std::array<Real,AMREX_SPACEDIM>& a_hi) noexcept
{
    AMREX_D_EXPR(xlo[0] = a_lo[0] , xlo[1] = a_lo[1] , xlo[2] = a_lo[2]);
    AMREX_D_EXPR(xhi[0] = a_hi[0] , xhi[1] = a_hi[1] , xhi[2] = a_hi[2]);
}

std::ostream&
operator << (std::ostream &os, const RealBox& b)
{
    os << "(RealBox ";
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        os << b.lo(i) << ' ' << b.hi(i) << ' ';
    os << ')';
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

std::istream&
operator >> (std::istream &is, RealBox& b)
{
    is.ignore(BL_IGNORE_MAX,'(');

    std::string s;

    is >> s;

    if (s != "RealBox")
    {
        amrex::ErrorStream() << "unexpected token in RealBox: " << s << '\n';
        amrex::Abort();
    }

    Real lo[AMREX_SPACEDIM];
    Real hi[AMREX_SPACEDIM];
#ifdef BL_USE_FLOAT
    double dlotemp, dhitemp;
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        is >> dlotemp >> dhitemp;
        lo[i] = dlotemp;
        hi[i] = dhitemp;
    }
#else
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        is >> lo[i] >> hi[i];
#endif

    is.ignore(BL_IGNORE_MAX, ')');

    b = RealBox(lo,hi);

    return is;
}

}
