//
// $Id: RealBox.cpp,v 1.9 2001-08-01 21:50:50 lijewski Exp $
//
#include <iostream>

#include <RealBox.H>
//
// The definition of lone static data member.
//
Real RealBox::eps = 1.0e-6;

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
    computeBoxLen();
}

void
RealBox::setEpsilon (Real epsilon)
{
    eps = epsilon;
}

Real
RealBox::epsilon ()
{
    return eps;
}

const Real*
RealBox::lo () const
{
    return xlo;
}

const Real*
RealBox::hi () const
{
    return xhi;
}

const Real*
RealBox::length () const
{
    return len;
}

Real
RealBox::lo (int dir) const
{
    return xlo[dir];
}

Real
RealBox::hi (int dir) const
{
    return xhi[dir];
}

Real
RealBox::length (int dir) const
{
    return len[dir];
}

void
RealBox::computeBoxLen ()
{
    D_EXPR(len[0] = xhi[0]-xlo[0],
           len[1] = xhi[1]-xlo[1],
           len[2] = xhi[2]-xlo[2]);
}

RealBox::RealBox ()
{
    D_TERM(xlo[0] , = xlo[1] , = xlo[2] ) = 0.;
    D_TERM(xhi[0] , = xhi[1] , = xhi[2] ) = -1.;
    computeBoxLen();
}

RealBox::RealBox (const Real* lo,
                  const Real* hi)
{
    D_EXPR(xlo[0] = lo[0] , xlo[1] = lo[1] , xlo[2] = lo[2]);
    D_EXPR(xhi[0] = hi[0] , xhi[1] = hi[1] , xhi[2] = hi[2]);
    computeBoxLen() ;
}

RealBox::RealBox (D_DECL(Real x0, Real y0, Real z0),
                  D_DECL(Real x1, Real y1, Real z1))
{
    D_EXPR(xlo[0] = x0 , xlo[1] = y0 , xlo[2] = z0);
    D_EXPR(xhi[0] = x1 , xhi[1] = y1 , xhi[2] = z1);
    computeBoxLen() ;
}

void
RealBox::setLo (const Real* lo)
{
    D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]);
    computeBoxLen();
}

void
RealBox::setLo (const Array<Real> &lo)
{
    D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]);
    computeBoxLen();
}

void
RealBox::setHi (const Real* hi)
{
    D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]);
    computeBoxLen();
}

void
RealBox::setHi (const Array<Real>& hi)
{
    D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]);
    computeBoxLen();
}

void
RealBox::setLo (int  indx,
                Real lo)
{
   BL_ASSERT(indx >= 0 && indx < BL_SPACEDIM);
   xlo[indx] = lo;
   computeBoxLen();
}

void
RealBox::setHi (int  indx,
                Real hi)
{
    BL_ASSERT(indx >= 0 && indx < BL_SPACEDIM);
    xhi[indx] = hi;
    computeBoxLen();
}

bool
RealBox::contains (const RealBox& rb) const
{
    return contains(rb.xlo) && contains(rb.xhi);
}

bool
RealBox::ok () const
{
    return (len[0] > eps)
#if (BL_SPACEDIM > 1)
        && (len[1] > eps)
#endif   
#if (BL_SPACEDIM > 2)
        && (len[2] > eps)
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
    for (int i = 0; i < BL_SPACEDIM; i++)
        is >> lo[i] >> hi[i];

    is.ignore(BL_IGNORE_MAX, ')');

    b = RealBox(lo,hi);

    return is;
}
