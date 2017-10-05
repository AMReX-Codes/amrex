
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>

#include <AMReX_BLassert.H>
#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_IndexType.H>
#include <AMReX_Utility.H>


namespace amrex {


const IntVect&
IntVect::TheUnitVector ()
{
    static const IntVect Unit(AMREX_D_DECL(1,1,1));
    return Unit;
}

const IntVect&
IntVect::TheZeroVector ()
{
    static const IntVect Zero(AMREX_D_DECL(0,0,0));
    return Zero;
}

//
// Static object initialization.
//
int IntVect::InitStatics()
{
  IntVect* pz = const_cast<IntVect*>( &IntVect::Zero );
  *pz = IntVect(AMREX_D_DECL(0,0,0));

  IntVect* pu = const_cast<IntVect*>( &IntVect::Unit );
  *pu = IntVect(AMREX_D_DECL(1,1,1));

  // No danger of IntVect::Zero and Unit not having been allocated, as ARM section
  // 3.4 says "The initialization of nonlocal static objects in a translation unit
  // is done before the first use of any function...defined in that translation
  // unit."
  //
  // Had to go through the const_cast stuff because it's nice to be able to declare
  // IntVect::Zero and IntVect::Unit as const.

  return 0; // arbitrary
}

const IntVect IntVect::Zero;
const IntVect IntVect::Unit;
static int s_dummyForIntVectCpp( IntVect::InitStatics() );

const IntVect&
IntVect::TheDimensionVector (int d)
{
    switch (d) {
    case (0) :
    {
	static const IntVect xdim(AMREX_D_DECL(1,0,0));
	return xdim;
    }
    case (1) :
    {
	static const IntVect ydim(AMREX_D_DECL(0,1,0));
	return ydim;
    }
    default:
    {
	static const IntVect zdim(AMREX_D_DECL(0,0,1));
	return zdim;
    }
    };
}

const IntVect&
IntVect::TheNodeVector ()
{
    static const IntVect Node(AMREX_D_DECL(IndexType::NODE,IndexType::NODE,IndexType::NODE));
    return Node;
}

const IntVect&
IntVect::TheCellVector ()
{
    static const IntVect Cell(AMREX_D_DECL(IndexType::CELL,IndexType::CELL,IndexType::CELL));
    return Cell;
}

const IntVect&
IntVect::TheMaxVector ()
{
    static const IntVect mx(AMREX_D_DECL(std::numeric_limits<int>::max(),
                                   std::numeric_limits<int>::max(),
                                   std::numeric_limits<int>::max()));
    return mx;
}

const IntVect&
IntVect::TheMinVector ()
{
    static const IntVect mn(AMREX_D_DECL(std::numeric_limits<int>::min(),
                                   std::numeric_limits<int>::min(),
                                   std::numeric_limits<int>::min()));
    return mn;
}

IntVect::IntVect (const int *a)
{
    AMREX_D_EXPR(vect[0] = a[0], vect[1] = a[1], vect[2] = a[2]);
}

IntVect::IntVect (const Vector<int> &a)
{
    BL_ASSERT(a.size() == BL_SPACEDIM);
    AMREX_D_EXPR(vect[0] = a[0], vect[1] = a[1], vect[2] = a[2]);
}

IntVect
min (const IntVect& p1,
	     const IntVect& p2)
{
    IntVect p(p1);
    p.min(p2);
    return p;
}

IntVect
max (const IntVect& p1,
     const IntVect& p2)
{
    IntVect p(p1);
    p.max(p2);
    return p;
}

IntVect
BASISV (int dir)
{
    BL_ASSERT(dir >= 0 && dir < BL_SPACEDIM);
    IntVect tmp;
    tmp[dir] = 1;
    return tmp;
}

IntVect
scale (const IntVect& p, int s)
{
    return IntVect(AMREX_D_DECL(s * p[0], s * p[1], s * p[2]));
}

IntVect
reflect (const IntVect& a,
         int            ref_ix,
         int            idir)
{
    BL_ASSERT(idir >= 0 && idir < BL_SPACEDIM);
    IntVect b(a);
    b[idir] = -b[idir] + 2*ref_ix;
    return b;
}

IntVect
diagShift (const IntVect& p, int s)
{
    return IntVect(AMREX_D_DECL(p[0] + s, p[1] + s, p[2] + s));
}

IntVect
coarsen (const IntVect& p,
         int            s)
{
    BL_ASSERT(s > 0);
    IntVect v = p;
    v.coarsen(IntVect(AMREX_D_DECL(s,s,s)));
    return v;
}

IntVect
coarsen (const IntVect& p1,
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
    return this->coarsen(IntVect(AMREX_D_DECL(s,s,s)));
}

IntVect&
IntVect::coarsen (const IntVect& p)
{
    BL_ASSERT(p.allGT(IntVect::TheZeroVector()));
    if (p != 1) {
        for (int i = 0; i <BL_SPACEDIM; ++i)
        {
            const int s = p.vect[i];
            vect[i] = ((vect[i]<0) ? -abs(vect[i]+1)/s-1 : vect[i]/s);
        }
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
    os << AMREX_D_TERM( '(' << p[0] , <<
                  ',' << p[1] , <<
                  ',' << p[2])  << ')';
    if (os.fail())
        amrex::Error("operator<<(ostream&,IntVect&) failed");
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

std::istream&
operator>> (std::istream& is,
            IntVect&      iv)
{
    is >> std::ws;
    char c;
    is >> c;

    if (c == '(')
    {
        AMREX_D_EXPR(is >> iv[0],
               is.ignore(BL_IGNORE_MAX, ',') >> iv[1],
               is.ignore(BL_IGNORE_MAX, ',') >> iv[2]);
        is.ignore(BL_IGNORE_MAX, ')');
    }
    else
    {
        amrex::Error("operator>>(istream&,IntVect&): expected \'(\'");
    }

    if (is.fail())
        amrex::Error("operator>>(istream&,IntVect&) failed");

    return is;
}

}
