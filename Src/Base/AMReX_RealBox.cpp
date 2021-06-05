
#include <AMReX_RealBox.H>
#include <AMReX_Algorithm.H>

#include <iostream>
#include <string>

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
        lo[i] = static_cast<Real>(dlotemp);
        hi[i] = static_cast<Real>(dhitemp);
    }
#else
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        is >> lo[i] >> hi[i];
#endif

    is.ignore(BL_IGNORE_MAX, ')');

    b = RealBox(lo,hi);

    return is;
}

bool AlmostEqual (const RealBox& box1,
                  const RealBox& box2,
                  Real eps /* = 0.0 */) noexcept
{
    bool almost_equal = true;
    for(int i = 0; i < AMREX_SPACEDIM && almost_equal; ++i)
    {
        almost_equal = almost_equal &&
            (std::abs(box1.lo(i) - box2.lo(i)) <= eps || amrex::almostEqual(box1.lo(i),box2.lo(i)));
        almost_equal = almost_equal &&
            (std::abs(box1.hi(i) - box2.hi(i)) <= eps || amrex::almostEqual(box1.hi(i),box2.hi(i)));
    }
    return almost_equal;
}

}
