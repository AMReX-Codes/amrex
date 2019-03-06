
#include <iostream>

#include <AMReX_IntVect.H>
#include <AMReX_BLassert.H>
#include <AMReX.H>
#include <AMReX_Utility.H>
#include <AMReX_IndexType.H>

namespace amrex {

const IntVect IntVect::Zero = IntVect::TheZeroVector();
const IntVect IntVect::Unit = IntVect::TheUnitVector();

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

#define BL_IGNORE_MAX 100000

std::istream&
operator>> (std::istream& is,
            IntVect&      iv)
{
    is >> std::ws;
    char c;
    is >> c;

    AMREX_D_TERM(iv[0]=0;, iv[1]=0;, iv[2]=0);

    if (c == '(')
    {
        is >> iv[0];
#if (AMREX_SPACEDIM >= 2)
        is >> std::ws;
        int ic = is.peek();
        if (ic == static_cast<int>(',')) {
            is.ignore(BL_IGNORE_MAX, ',');
            is >> iv[1];
#if (AMREX_SPACEDIM == 3)
            is >> std::ws;
            ic = is.peek();
            if (ic == static_cast<int>(',')) {
                is.ignore(BL_IGNORE_MAX, ',');
                is >> iv[2];
            }
#endif
        }
#endif
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
