
#include <iostream>

#include <AMReX_IntVect.H>
#include <AMReX_BLassert.H>
#include <AMReX.H>
#include <AMReX_Utility.H>
#include <AMReX_IndexType.H>

namespace amrex {

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
