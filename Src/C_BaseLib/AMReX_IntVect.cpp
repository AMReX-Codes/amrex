
#include <winstd.H>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>

#include <BLassert.H>
#include <BoxLib.H>
#include <IntVect.H>
#include <IndexType.H>

const IntVect&
IntVect::TheUnitVector ()
{
    static const IntVect Unit(D_DECL(1,1,1));
    return Unit;
}

const IntVect&
IntVect::TheZeroVector ()
{
    static const IntVect Zero(D_DECL(0,0,0));
    return Zero;
}

const IntVect&
IntVect::TheDimensionVector (int d)
{
    switch (d) {
    case (0) :
    {
	static const IntVect xdim(D_DECL(1,0,0));
	return xdim;
    }
    case (1) :
    {
	static const IntVect ydim(D_DECL(0,1,0));
	return ydim;
    }
    default:
    {
	static const IntVect zdim(D_DECL(0,0,1));
	return zdim;
    }
    };
}

const IntVect&
IntVect::TheNodeVector ()
{
    static const IntVect Node(D_DECL(IndexType::NODE,IndexType::NODE,IndexType::NODE));
    return Node;
}

const IntVect&
IntVect::TheCellVector ()
{
    static const IntVect Cell(D_DECL(IndexType::CELL,IndexType::CELL,IndexType::CELL));
    return Cell;
}

const IntVect&
IntVect::TheMaxVector ()
{
    static const IntVect mx(D_DECL(std::numeric_limits<int>::max(),
                                   std::numeric_limits<int>::max(),
                                   std::numeric_limits<int>::max()));
    return mx;
}

const IntVect&
IntVect::TheMinVector ()
{
    static const IntVect mn(D_DECL(std::numeric_limits<int>::min(),
                                   std::numeric_limits<int>::min(),
                                   std::numeric_limits<int>::min()));
    return mn;
}

IntVect::IntVect (const int *a)
{
    D_EXPR(vect[0] = a[0], vect[1] = a[1], vect[2] = a[2]);
}

IntVect::IntVect (const Array<int> &a)
{
    BL_ASSERT(a.size() == BL_SPACEDIM);
    D_EXPR(vect[0] = a[0], vect[1] = a[1], vect[2] = a[2]);
}

bool
IntVect::lexLT (const IntVect &s) const
{
#if (BL_SPACEDIM == 1)
    return vect[0] < s[0];
#elif (BL_SPACEDIM == 2)
    return (vect[1] < s[1]) || ((vect[1] == s[1]) && (vect[0] < s[0]));
#else
    return (vect[2] < s[2]) || ((vect[2] == s[2]) &&
	(  (vect[1] < s[1]) || ((vect[1] == s[1]) && (vect[0] < s[0])) ));
#endif
}

bool
IntVect::lexGT (const IntVect& s) const
{
#if (BL_SPACEDIM == 1)
    return vect[0] > s[0];
#elif (BL_SPACEDIM == 2)
    return (vect[1] > s[1]) || ((vect[1] == s[1]) && (vect[0] > s[0]));
#else
    return (vect[2] > s[2]) || ((vect[2] == s[2]) &&
	(  (vect[1] > s[1]) || ((vect[1] == s[1]) && (vect[0] > s[0])) ));
#endif
}

IntVect
BoxLib::min (const IntVect& p1,
	     const IntVect& p2)
{
    IntVect p(p1);
    p.min(p2);
    return p;
}

IntVect
BoxLib::max (const IntVect& p1,
	     const IntVect& p2)
{
    IntVect p(p1);
    p.max(p2);
    return p;
}

IntVect
BoxLib::BASISV (int dir)
{
    BL_ASSERT(dir >= 0 && dir < BL_SPACEDIM);
    IntVect tmp;
    tmp[dir] = 1;
    return tmp;
}

IntVect
BoxLib::scale (const IntVect& p, int s)
{
    return IntVect(D_DECL(s * p[0], s * p[1], s * p[2]));
}

IntVect
BoxLib::reflect (const IntVect& a,
		 int            ref_ix,
		 int            idir)
{
    BL_ASSERT(idir >= 0 && idir < BL_SPACEDIM);
    IntVect b(a);
    b[idir] = -b[idir] + 2*ref_ix;
    return b;
}

IntVect
BoxLib::diagShift (const IntVect& p, int s)
{
    return IntVect(D_DECL(p[0] + s, p[1] + s, p[2] + s));
}

IntVect
BoxLib::coarsen (const IntVect& p,
		 int            s)
{
    BL_ASSERT(s > 0);
    IntVect v = p;
    v.coarsen(IntVect(D_DECL(s,s,s)));
    return v;
}

IntVect
BoxLib::coarsen (const IntVect& p1,
		 const IntVect& p2)
{
    IntVect v = p1;
    v.coarsen(p2);
    return v;
}

IntVect&
IntVect::coarsen (int s)
{
    BL_ASSERT(s > 0);
    return this->coarsen(IntVect(D_DECL(s,s,s)));
}

IntVect&
IntVect::coarsen (const IntVect& p)
{
    BL_ASSERT(p > IntVect::TheZeroVector());
    for (int i = 0; i <BL_SPACEDIM; ++i)
    {
        const int s = p.vect[i];
        vect[i] = ((vect[i]<0) ? -abs(vect[i]+1)/s-1 : vect[i]/s);
    }
    return *this;
}

//
// Returns IntVect which is the componentwise integer projection
// of IntVect p1 by IntVect p2.
//

std::ostream&
operator<< (std::ostream&  os,
            const IntVect& p)
{
    os << D_TERM( '(' << p[0] , <<
                  ',' << p[1] , <<
                  ',' << p[2])  << ')';
    if (os.fail())
        BoxLib::Error("operator<<(ostream&,IntVect&) failed");
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000
#include <Utility.H>

std::istream&
operator>> (std::istream& is,
            IntVect&      iv)
{
    is >> std::ws;
    char c;
    is >> c;

    if (c == '(')
    {
        D_EXPR(is >> iv[0],
               is.ignore(BL_IGNORE_MAX, ',') >> iv[1],
               is.ignore(BL_IGNORE_MAX, ',') >> iv[2]);
        is.ignore(BL_IGNORE_MAX, ')');
    }
    else
    {
        BoxLib::Error("operator>>(istream&,IntVect&): expected \'(\'");
    }

    if (is.fail())
        BoxLib::Error("operator>>(istream&,IntVect&) failed");

    return is;
}
