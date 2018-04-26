
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
