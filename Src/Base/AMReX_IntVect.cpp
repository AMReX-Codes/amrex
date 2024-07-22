
#include <AMReX_IntVect.H>
#include <AMReX_BLassert.H>
#include <AMReX.H>
#include <AMReX_Utility.H>
#include <AMReX_IndexType.H>

#include <iostream>

namespace amrex::detail {

std::ostream&
int_vector_write (std::ostream& os, const int* iv, int dim)
{
    os << '(' << iv[0];
    for (int i=1; i<dim; ++i) {
        os << ',' << iv[i];
    }
    os << ')';
    if (os.fail()) {
        amrex::Error("operator<<(ostream&,IntVect&) failed");
    }
    return os;
}

#define BL_IGNORE_MAX 100000

std::istream&
int_vector_read (std::istream& is, int* iv, int dim)
{
    is >> std::ws;
    char c;
    is >> c;

    for (int i=0; i<dim; ++i) {
        iv[i] = 0;
    }

    if (c == '(')
    {
        is >> iv[0];
        for (int i=1; i<dim; ++i) {
            is >> std::ws;
            int ic = is.peek();
            if (ic == static_cast<int>(',')) {
                is.ignore(BL_IGNORE_MAX, ',');
                is >> iv[i];
                continue;
            }
            break;
        }
        is.ignore(BL_IGNORE_MAX, ')');
    }
    else
    {
        amrex::Error("operator>>(istream&,IntVect&): expected \'(\'");
    }

    if (is.fail()) {
        amrex::Error("operator>>(istream&,IntVect&) failed");
    }

    return is;
}

}
