//BL_COPYRIGHT_NOTICE

//
// $Id: RealBox.cpp,v 1.2 1997-12-11 05:01:10 lijewski Exp $
//

#include <aString.H>
#include <Misc.H>
#include <Utility.H>
#include <RealBox.H>

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

bool
RealBox::contains (const Real* point)
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

ostream&
operator << (ostream &os, const RealBox& b)
{
    os << "(RealBox ";
    for (int i = 0; i < BL_SPACEDIM; i++)
        os << b.xlo[i] << ' ' << b.xhi[i] << ' ';
    os << ')';
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
    for (int i = 0; i < BL_SPACEDIM; i++)
        is >> b.xlo[i] >> b.xhi[i];
    is.ignore(BL_IGNORE_MAX, ')');
    b.computeBoxLen();
    return is;
}
