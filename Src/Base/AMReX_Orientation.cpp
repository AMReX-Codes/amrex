
#include <iostream>

#include <AMReX.H>
#include <AMReX_Orientation.H>

namespace amrex {

std::ostream&
operator<< (std::ostream&      os,
            const Orientation& o)
{
    os << '('<< int(o) << ')' ;
    if (os.fail())
        amrex::Error("operator<<(ostream&,Orientation&) failed");
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

std::istream&
operator>> (std::istream& is,
            Orientation&  o)
{
    char c;
    is >> c;

    if (c == '(')
    {
        is >> o.val;
        is.ignore(BL_IGNORE_MAX, ')');
    }
    else
    {
        amrex::Error("operator>>(istream&,Orientation&): expected \'(\'");
    }

    if (is.fail())
        amrex::Error("operator>>(ostream&,Orientation&) failed");

    return is;
}

}
