//BL_COPYRIGHT_NOTICE

//
// $Id: RealBox.cpp,v 1.1 1997-12-10 19:07:43 lijewski Exp $
//

#include <aString.H>
#include <Misc.H>
#include <Utility.H>
#include <RealBox.H>

Real RealBox::eps = 1.0e-6;

RealBox::RealBox()
{
    D_TERM(xlo[0] , = xlo[1] , = xlo[2] ) = 0.;
    D_TERM(xhi[0] , = xhi[1] , = xhi[2] ) = -1.;
    computeBoxLen();
}

RealBox::RealBox(const Real* lo, const Real* hi)
{
    D_EXPR(xlo[0] = lo[0] , xlo[1] = lo[1] , xlo[2] = lo[2]);
    D_EXPR(xhi[0] = hi[0] , xhi[1] = hi[1] , xhi[2] = hi[2]);
    computeBoxLen() ;
}

RealBox::RealBox(const Box& bx, const Real* dx, const Real* base)
{
    const int* lo = bx.loVect();
    const int* hi = bx.hiVect();
    int i;
    for (i = 0; i < BL_SPACEDIM; i++) {
	xlo[i] = base[i] + dx[i]*lo[i];
	int shft = (bx.type(i)==IndexType::CELL ? 1 : 0);
	xhi[i] = base[i] + dx[i]*(hi[i]+ shft);
    }   
    computeBoxLen() ;
}

RealBox::RealBox(D_DECL(Real x0, Real y0, Real z0),
		 D_DECL(Real x1, Real y1, Real z1))
{
    D_EXPR(xlo[0] = x0 , xlo[1] = y0 , xlo[2] = z0);
    D_EXPR(xhi[0] = x1 , xhi[1] = y1 , xhi[2] = z1);
    computeBoxLen() ;
}

void RealBox::setLo(const Real* lo)
{
    D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]);
    computeBoxLen();
}

void RealBox::setLo(const Array<Real> &lo)
{
    D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]);
    computeBoxLen();
}

void RealBox::setHi(const Real* hi)
{
    D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]);
    computeBoxLen();
}

void RealBox::setHi(const Array<Real>& hi)
{
    D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]);
    computeBoxLen();
}

void
RealBox::setLo(int indx, Real lo)
{
   assert( indx >= 0 && indx < BL_SPACEDIM);
   xlo[indx] = lo;
   computeBoxLen();
}

void
RealBox::setHi(int indx, Real hi)
{
    assert( indx >= 0 && indx < BL_SPACEDIM);
    xhi[indx] = hi;
    computeBoxLen();
}

int
RealBox::contains(const Real* point)
{
    return  (xlo[0]-eps < point[0]) && (point[0] < xhi[0]+eps)
#if (BL_SPACEDIM > 1)   
        && (xlo[1]-eps < point[1]) && (point[1] < xhi[1]+eps)
#endif
#if (BL_SPACEDIM > 2)   
        && (xlo[2]-eps < point[2]) && (point[2] < xhi[2]+eps)
#endif
   ;
}

int
RealBox::contains(const RealBox& rb)
{
    return (contains(rb.xlo) && contains(rb.xhi));
}

int
RealBox::ok() const
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

ostream&
operator << (ostream &os, const RealBox& b)
{
	os << "(RealBox ";
        int i;
	for (i = 0; i < BL_SPACEDIM; i++)
	    os << b.xlo[i] << ' ' << b.xhi[i] << ' ';
	os << ')' << flush;
    return os;
}

istream&
operator >> (istream &is, RealBox& b)
{
	is.ignore(BL_IGNORE_MAX,'(');
	aString s;
	is >> s;
	if (s != "RealBox")
        {
  	    cerr << "unexpected token in RealBox: " << s << '\n';
            BoxLib::Abort();
	}
        int i;
	for (i = 0; i < BL_SPACEDIM; i++)
	    is >> b.xlo[i] >> b.xhi[i];
	is.ignore(BL_IGNORE_MAX, ')');
	b.computeBoxLen();
    return is;
}
