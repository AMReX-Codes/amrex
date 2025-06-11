#include <AMReX_RealVect.H>
#include <iostream>

namespace amrex::detail {

std::ostream&
real_vector_write (std::ostream& os, const Real* p, int dim)
{
    os << '(' << p[0];
    for (int i=1; i<dim; ++i) {
        os << ',' << p[i];
    }
    os << ')';
    if (os.fail()) {
        amrex::Error("operator<<(ostream&,RealVect&) failed");
    }
    return os;
}

#define BL_IGNORE_MAX 100000

std::istream&
real_vector_read (std::istream& is, Real* p, int dim)
{
    is >> std::ws;
    char c;
    is >> c;

    for (int i=0; i<dim; ++i) {
        p[i] = 0;
    }

    if (c == '(')
    {
        is >> p[0];
        for (int i=1; i<dim; ++i) {
            is >> std::ws;
            int ic = is.peek();
            if (ic == static_cast<int>(',')) {
                is.ignore(BL_IGNORE_MAX, ',');
                is >> p[i];
                continue;
            }
            break;
        }
        is.ignore(BL_IGNORE_MAX, ')');
    }
    else
    {
        amrex::Error("operator>>(istream&,RealVect&): expected \'(\'");
    }

    if (is.fail()) {
        amrex::Error("operator>>(istream&,RealVect&) failed");
    }

    return is;
}

} //namespace amrex::detail
