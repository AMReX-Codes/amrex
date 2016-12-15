
#include <iostream>
#include <string>

#include <RealBox.H>
//
// The definition of lone static data member.
//
Real RealBox::eps = 1.0e-8;

RealBox::RealBox (const Box&  bx,
                  const Real* dx,
                  const Real* base)
{
    const int* lo = bx.loVect();
    const int* hi = bx.hiVect();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        xlo[i] = base[i] + dx[i]*lo[i];
        int shft = (bx.type(i) == IndexType::CELL ? 1 : 0);
        xhi[i] = base[i] + dx[i]*(hi[i]+ shft);
    }   
}

RealBox::RealBox ()
{
    D_TERM(xlo[0] , = xlo[1] , = xlo[2] ) = 0.;
    D_TERM(xhi[0] , = xhi[1] , = xhi[2] ) = -1.;
}

RealBox::RealBox (const Real* lo,
                  const Real* hi)
{
    D_EXPR(xlo[0] = lo[0] , xlo[1] = lo[1] , xlo[2] = lo[2]);
    D_EXPR(xhi[0] = hi[0] , xhi[1] = hi[1] , xhi[2] = hi[2]);
}

RealBox::RealBox (D_DECL(Real x0, Real y0, Real z0),
                  D_DECL(Real x1, Real y1, Real z1))
{
    D_EXPR(xlo[0] = x0 , xlo[1] = y0 , xlo[2] = z0);
    D_EXPR(xhi[0] = x1 , xhi[1] = y1 , xhi[2] = z1);
}

bool
RealBox::contains (const RealBox& rb) const
{
    return contains(rb.xlo) && contains(rb.xhi);
}

bool
RealBox::ok () const
{
    return (length(0) > eps)
#if (BL_SPACEDIM > 1)
        && (length(1) > eps)
#endif   
#if (BL_SPACEDIM > 2)
        && (length(2) > eps)
#endif
   ;
}

bool
RealBox::contains (const Real* point) const
{
    return  D_TERM((xlo[0]-eps < point[0]) && (point[0] < xhi[0]+eps),
                   && (xlo[1]-eps < point[1]) && (point[1] < xhi[1]+eps),
                   && (xlo[2]-eps < point[2]) && (point[2] < xhi[2]+eps));
}

std::ostream&
operator << (std::ostream &os, const RealBox& b)
{
    os << "(RealBox ";
    for (int i = 0; i < BL_SPACEDIM; i++)
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
        std::cerr << "unexpected token in RealBox: " << s << '\n';
        BoxLib::Abort();
    }

    Real lo[BL_SPACEDIM];
    Real hi[BL_SPACEDIM];
#ifdef BL_USE_FLOAT
    double dlotemp, dhitemp;
    for (int i = 0; i < BL_SPACEDIM; i++) {
        is >> dlotemp >> dhitemp;
        lo[i] = dlotemp;
        hi[i] = dhitemp;
    }
#else
    for (int i = 0; i < BL_SPACEDIM; i++)
        is >> lo[i] >> hi[i];
#endif

    is.ignore(BL_IGNORE_MAX, ')');

    b = RealBox(lo,hi);

    return is;
}
